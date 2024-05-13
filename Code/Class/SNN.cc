#include "SNN.h"

using namespace std;

SNN::SNN(int _NL0, int _NL1,
         float _alpha,
         float _CFI0, float _CFI1, float _CF01,
         float _L1inhibitfactor,
         float _K, float _K1, float _K2,
         float _IE_Pot_const, double _IPSP_dt_dilation,
         double _MaxDelay,

         double _tau_m, double _tau_s, double _tau_r, double _tau_plus, double _tau_minus,
         double _a_plus, double _a_minus,

         double _taud_plus, double _taud_minus,
         double _d_plus, double _d_minus,

         int _N_InputStreams,
         float _Threshold0, float _Threshold1, float _sparsity, bool _split_layer0) :
                                                 // Initializations
                                                 alpha(_alpha),

                                                 CFI0(_CFI0),
                                                 CFI1(_CFI1),
                                                 CF01(_CF01),

                                                 L1inhibitfactor(_L1inhibitfactor),

                                                 K(_K),
                                                 K1(_K1),
                                                 K2(_K2),

                                                 IE_Pot_const(_IE_Pot_const),

                                                 IPSP_dt_dilation(_IPSP_dt_dilation),

                                                 MaxDelay(_MaxDelay),

                                                 tau_m(_tau_m), // membrane time constant (the potential will decrease ~ exp(-t/tau_m))
                                                 tau_s(_tau_s), // synaptic time constant
                                                 tau_r(_tau_r), // refractory time constant
                                                 tau_plus(_tau_plus),
                                                 tau_minus(_tau_minus),

                                                 a_plus(_a_plus),
                                                 a_minus(_a_minus),

                                                 taud_plus(_taud_plus),
                                                 taud_minus(_taud_minus),

                                                 d_plus(_d_plus),
                                                 d_minus(_d_minus),


                                                 N_InputStreams(_N_InputStreams),
                                                 sparsity(_sparsity),
                                                 split_layer0(_split_layer0)

{
    Threshold[0] = _Threshold0;
    Threshold[1] = _Threshold1;
    
    N_neuronsL[0] = _NL0;
    N_neuronsL[1] = _NL1;
    N_neurons = N_neuronsL[0] + N_neuronsL[1];
    N_streams = N_InputStreams + _NL0;
    tmax = tau_s * tau_m / (tau_m - tau_s) * (log(tau_m) - log(tau_s));
    MaxDeltaT = 7. * tau_m;

    fire_granularity = tau_s / 5.;
    fire_precision = min(Threshold[0], Threshold[1]) *2.5 / 100.;
    myRNG = new TRandom3(static_cast<unsigned int>(std::time(0)));
    largenumber = 999999999.;
    epsilon = 1. / largenumber;

    Weight = new float *[N_neurons];         // Weight of synapse-neuron strength
    Weight_initial = new float *[N_neurons];
    Delay = new double *[N_neurons];          // Delay in incoming signals
    Delay_initial = new double *[N_neurons];
    check_LTD = new bool *[N_neurons];       // checks to generate LTD after neuron discharge
    Void_weight = new bool *[N_neurons];     // These may be used to model disconnections
    
    EnableIPSP = new bool *[N_neurons];

    for (int in = 0; in < N_neurons; in++)
    {
        Weight[in] = new float[N_streams];
        Weight_initial[in] = new float[N_streams];
        check_LTD[in] = new bool[N_streams];
        Void_weight[in] = new bool[N_streams];
        Delay[in] = new double[N_streams];
        Delay_initial[in] = new double[N_streams];
        EnableIPSP[in] = new bool[N_neurons];
    }

    History_time = new vector<double>[N_neurons]; // Time of signal events per each 1neuron
    History_type = new vector<int>[N_neurons];   // Type of signal
    History_ID = new vector<int>[N_neurons];     // ID of generating signal stream or neuron
    History_class = new vector<int>[N_neurons];  // Class of the signal

    Fire_time = new vector<double>[N_neurons];    // Times of firing of each neuron
    Neuron_layer = new int[N_neurons];
    sumweight = new float[N_neurons]; // summed weights of streams for each neurons for the purpose of normalization
    sumdelays = new double[N_neurons]; // summed delays of streams for each neurons for the purpose of normalization


    Delta_delay = MaxDelay/10.;
    Mean_delay = MaxDelay/2.;

    Init_neurons();
    Init_connection_map();
    Init_weights();
    Init_delays_uniform();
}

SNN::~SNN()
{
}

void SNN::Reset_weights(){
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_InputStreams; is++)
            Weight[in][is] = Weight_initial[in][is];
    }
}

double SNN::bisectionMethod(double a, double b, int in, double epsilon, std::function<float(int, double, bool)> func)
{
    float fa = func(in, a, false);
    float fb = func(in, b, false);
    double c = 0;

    if (fa * fb > 0) return a; 


    int maxIterations = 10; // Choose an appropriate maximum number of iterations

    for (int i = 0; i < maxIterations; i++)
    {
        c = (a + b) / 2;
        if (fa * fb > 0) return a;

        float fc = func(in, c, false);

        // Check if the root is found within the specified tolerance
        if (std::abs(fc) < epsilon)
        {
            return c;
        }

        // Update the interval based on the sign of the function at the midpoint
        if (fa * fc < 0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    return c;
}

// Initialize neuron potentials
// ----------------------------
void SNN::Init_neurons()
{
    for (int in = 0; in < N_neurons; in++)
    {
        // Set first event in history of this neuron
        History_time[in].clear();
        History_type[in].clear();
        History_ID[in].clear();
        History_class[in].clear();

        History_time[in].push_back(0);
        History_type[in].push_back(SPIKE);
        History_ID[in].push_back(0);
        History_class[in].push_back(-2);

        if (in < N_neuronsL[0])
            Neuron_layer[in] = 0;
        else
            Neuron_layer[in] = 1;
    }
    return;
}

// Set same synapse weights
// --------------------------
void SNN::Set_weights()
{
    for (int in = 0; in < N_neurons; in++)
    {
        sumweight[in] = 0;
        for (int is = 0; is < N_streams; is++)
        {
            check_LTD[in][is] = true; // flags used to see if we need to create a LTD signal after a neuron discharge
            if (Void_weight[in][is])
                Weight[in][is] = -1;
            else
            {
                Weight[in][is] = 1;
                sumweight[in] += Weight[in][is];
            }
        }
    }

    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
            if (sumweight[in] > 0 && !Void_weight[in][is])
                Weight[in][is] = Weight[in][is] / sumweight[in];
            Weight_initial[in][is] = Weight[in][is];
        }
    }

    return;
}

void SNN::Init_delays_man(){
    //we assume to read a specific file
    string filename = string(getenv("SNN_PATH"))+"/Code/Data/delays.txt";

    // Open the file
    ifstream file(filename);

    if (!file.is_open()) throw runtime_error("Error: Unable to open the file " + filename); 

    //set them by default to MaxDelay for connections involving the input stream
    for (int in = 0; in < N_neuronsL[0]; in++){
        for (int is = 0; is < N_InputStreams; is++){
            if (!(file >> Delay[in][is] && file >> Delay_initial[in][is])) {
                file.close(); // Close the file before throwing exception
                throw runtime_error("Error: Unable to read delay value at position (" + to_string(in) + ", " + to_string(is) + ").");    
            }
        }
    }
    for (int in = N_neuronsL[0]; in < N_neurons; in++){
        for (int is = 0; is < N_InputStreams; is++){
            Delay[in][is] = 0;
            Delay_initial[in][is] = 0;  

        }      
    }
    // Close the file
    file.close();

    //set them by default to 0 for connections involving layers
    for (int in = 0; in < N_neurons; in++){
        for (int is = N_InputStreams; is < N_streams; is++){
            Delay[in][is] = 0;
            Delay_initial[in][is] = 0;
        }
    }
    return;
}

void SNN::Init_delays_gauss()
{
    //gaussian delays for all synapses connecting to input streams
    for (int in = 0; in < N_neurons; in++){
        for (int is = 0; is < N_InputStreams; is++){
            Delay[in][is] = myRNG->Gaus(0.5*MaxDelay, sparsity/sqrt(N_InputStreams)*MaxDelay);
            Delay_initial[in][is] = Delay[in][is];
        }
    }
    //0 delay for everything else
    for (int in = 0; in < N_neurons; in++){
        for (int is = N_InputStreams; is < N_streams; is++){
            Delay[in][is] = 0;
            Delay_initial[in][is] = 0;
        }
    }

    return;
}

void SNN::Init_weights()
{
    for (int in = 0; in < N_neurons; in++)
    {
        sumweight[in] = 0;
        for (int is = 0; is < N_streams; is++)
        {
            check_LTD[in][is] = true; // flags used to see if we need to create a LTD signal after a neuron discharge
            if (Void_weight[in][is])
                Weight[in][is] = -1;
            else
            {
                Weight[in][is] = myRNG->Gaus(1, sparsity/sqrt(N_InputStreams));
                if(Weight[in][is]<0) Weight[in][is]=0;
                sumweight[in] += Weight[in][is];
            }
        }
    }

    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
            if (sumweight[in] > 0 && !Void_weight[in][is])
                Weight[in][is] = Weight[in][is] / sumweight[in];
            Weight_initial[in][is] = Weight[in][is];    
        }
    }
    return;
}

void SNN::Init_delays_uniform()
{
   for (int in = 0; in < N_neurons; in++){
        sumdelays[in] = 0.;
        for (int is = 0; is < N_InputStreams; is++){
            Delay[in][is] = Delta_delay * (2*myRNG->Uniform()-1.) + Mean_delay;
            Delay_initial[in][is] = Delay[in][is];
            sumdelays[in]+=Delay[in][is];
        }
    }
    //0 delay for everything else
    for (int in = 0; in < N_neurons; in++){
        for (int is = N_InputStreams; is < N_streams; is++){
            Delay[in][is] = 0.;
            Delay_initial[in][is] = 0.;
            sumdelays[in]+=Delay[in][is];
        }
        //cout << in << " " << sumdelays[in] << endl;
    }

    return;
}

void SNN::Init_weights_uniform()
{
    for (int in = 0; in < N_neurons; in++)
    {
        sumweight[in] = 0;
        for (int is = 0; is < N_streams; is++)
        {
            check_LTD[in][is] = true; // flags used to see if we need to create a LTD signal after a neuron discharge
            if (Void_weight[in][is])
                Weight[in][is] = -1;
            else
            {
                Weight[in][is] = myRNG->Uniform();
                sumweight[in] += Weight[in][is];
            }
        }
    }

    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
            if (sumweight[in] > 0 && !Void_weight[in][is])
                Weight[in][is] = Weight[in][is] / sumweight[in];
            Weight_initial[in][is] = Weight[in][is];    
        }
    }
    return;
}

// Initialize connection map
// -------------------------
void SNN::Init_connection_map()
{

    // Setting L0 input connections
    for (int in = 0; in < N_neuronsL[0]; in++)
    {
        // input connections Tracking layers -> L0
        for (int is = 0; is < N_InputStreams; is++)
        {
            Void_weight[in][is] = false;
            if (myRNG->Uniform() > CFI0)
                Void_weight[in][is] = true;
        }
        // input connections L0 -> L0
        for (int is = N_InputStreams; is < N_streams; is++)
            Void_weight[in][is] = true;
    }

    // Setting L1 input connections
    for (int in = N_neuronsL[0]; in < N_neurons; in++)
    {
        // input connections Tracking layers -> L1
        for (int is = 0; is < N_InputStreams; is++)
        {
            Void_weight[in][is] = false;
            if (myRNG->Uniform() > CFI1)
                Void_weight[in][is] = true;
        }
        // input connections L0 -> L1
        for (int is = N_InputStreams; is < N_streams; is++)
        {
            Void_weight[in][is] = false;
            if (myRNG->Uniform() > CF01)
                Void_weight[in][is] = true;
        }
    }

    //initialize IPSP matrix
    for (int in=0; in < N_neurons; in++){
        for (int jn=0; jn < N_neurons; jn++){
            EnableIPSP[in][jn] = false;
        }
    }

    //enabFle IPSP in layer 1
    for(int in=N_neuronsL[0]; in<N_neurons; in++){
        for(int jn=N_neuronsL[0]; jn<N_neurons; jn++)
            EnableIPSP[in][jn] = true;
    }

    //enable IPSP in layer 0
    for (int in=0; in < N_neuronsL[0]; in++){
        //in this case we create 2 sublayers out of layer 0
        if(split_layer0){
            for (int jn=0; jn < N_neuronsL[0]; jn++){
                if((jn-N_neuronsL[0]/2)*(in-N_neuronsL[0]/2)>0)
                    EnableIPSP[in][jn] = true;
            }   
        }

        //otherwise we enable IPSP across all the layer
        else{
            for (int jn=0; jn < N_neuronsL[0]; jn++)
                EnableIPSP[in][jn] = true;
        }
    }
    
    return;
}

//Function to insert new spikes in the correct temporal position
void SNN::insert_spike(int id_neuron, double spike_time, int type, int id, int spike_class){
    // Find the position where the new value should be inserted
    auto it = lower_bound(History_time[id_neuron].begin(), History_time[id_neuron].end(), spike_time);
    int position = distance(History_time[id_neuron].begin(), it);
    // Insert the new value at the found position
    History_time[id_neuron].insert(it, spike_time);
    History_type[id_neuron].insert(History_type[id_neuron].begin()+position, type);
    History_ID[id_neuron].insert(History_ID[id_neuron].begin()+position, id);
    History_class[id_neuron].insert(History_class[id_neuron].begin()+position, type);
}

// Model Excitatory Post-Synaptic Potential
// We take this as parametrized in T. Masquelier et al., "Competitive STDP-Based Spike Pattern Learning", DOI: 10.1162/neco.2008.06-08-804
// ---------------------------------------------------------------------------------------------------------------------------------------
float SNN::EPS_potential(double delta_t)
{
    float psp = 0.;
    if (delta_t >= 0. && delta_t < MaxDeltaT)
        psp = K * (exp(-delta_t / tau_m) - exp(-delta_t / tau_s));
    return psp;
}

// Model membrane potential after spike
// Also modeled as in paper cited above, like IPSP and other signals below
// -----------------------------------------------------------------------
float SNN::Spike_potential(double delta_t, int ilayer)
{
    float sp = 0.;
    if (delta_t >= 0. && delta_t < MaxDeltaT)
        sp = Threshold[ilayer] * (K1 * exp(-delta_t / tau_m) - K2 * (exp(-delta_t / tau_m) - exp(-delta_t / tau_s)));
    return sp;
}

// Model Inhibitory Post-Synaptic Potential (IPSP)
// -----------------------------------------------
float SNN::Inhibitory_potential(double delta_t, int ilayer)
{
    // In order to dilate the inhibition to larger times (in the attempt at obtaining higher selectivity),
    // we kludge it by multiplying delta_t and maxdeltat by a factor
    delta_t = delta_t * IPSP_dt_dilation;
    float ip = 0.;
    float thisalpha = alpha;
    if (ilayer > 0)
        thisalpha = L1inhibitfactor * alpha; // Different inhibition in L1
    if (delta_t >= 0. && delta_t < MaxDeltaT)
        ip = -thisalpha * Threshold[ilayer] * EPS_potential(delta_t);
    return ip;
}

// Compute collective effect of excitatory, post-spike, and inhibitory potentials on a neuron
// ------------------------------------------------------------------------------------------
double SNN::Neuron_firetime(int in, double t)
{
    double t0 = History_time[in][0];
    double delta_t = t - t0;
    if (delta_t < tau_r)
        return largenumber;

    int ilayer = Neuron_layer[in];
    // now we will scan the interval in between the last EPSP and this time looking for an activation according to the defined granularity
    double last_EPSP = -1;

    for (int ih = History_type[in].size() - 1; ih > 1; ih--)
    {
        // longer approach: add "&& History_ID[in][ih] < N_InputStreams" to rescan from the last InputStream spike
        if (History_type[in][ih] == EPSP && !Void_weight[in][History_ID[in][ih]] && History_time[in][ih]<t)
        {
            last_EPSP = History_time[in][ih];
            break;
        }
    }
    // in that case it's impossible that the neuron is firing
    if (last_EPSP < 0)
        return largenumber;

    // now I want to scan the potential from the last EPSP to time t at fire_granularity steps
    double t_neg;
    if(last_EPSP - t0 > tau_r) t_neg = last_EPSP;
    else(t_neg) = t0+tau_r;

    double P_neg = Neuron_Potential(in, t_neg, true);
    if (P_neg > Threshold[ilayer]) return t_neg;
    
    double time = t_neg + fire_granularity;
    float P_t;
    bool fire = false;
    while (time < t)
    {
        P_t = Neuron_Potential(in, time, false);
        // If I have a value below the threshold I save it
        if (P_t < Threshold[ilayer])
        {
            P_neg = P_t;
            t_neg = time;
            time += fire_granularity;
        }
        // if I find a value higher than the treshold a neuron is gonna fire!
        else
        {
            fire = true;
            break;
        }
    }
    // Let's check at the time t
    if (!fire)
    {
        time = t;
        P_t = Neuron_Potential(in, time, false);
        // if still the potential is below the threshold we know that the neuron is not activating
        if (P_t < Threshold[ilayer]) return largenumber;
    }

    // if we are here the potential has reached the threshold at some point!
    // we need to determine when the neuron has fired given a certain confidence
    return bisectionMethod(t_neg, time, in, fire_precision,
                           [this](int in, double time, bool delete_history)
                           {
                               return Neuron_Potential(in, time, delete_history) - Threshold[Neuron_layer[in]];
                           });
}

//Handle the activation of a neuron
void SNN::Activate_Neuron(int in, double t){
    // Reset history of this neuron
    //TODO: loop pack in History_time to clear just the spikes before the activation
    auto it = upper_bound(History_time[in].begin(), History_time[in].end(), t);
    int position = distance(History_time[in].begin(), it);
    // Erase elements before the position found by binary search
    History_time[in].erase(History_time[in].begin(), it);
    History_type[in].erase(History_type[in].begin(), History_type[in].begin() + position);
    History_ID[in].erase(History_ID[in].begin(), History_ID[in].begin() + position);

    /*
    History_time[in].clear();
    History_type[in].clear();
    History_ID[in].clear();
    */

    insert_spike(in, t, SPIKE, 0, NOCLASS);
    
    return;
}

// Compute collective effect of excitatory, post-spike, and inhibitory potentials on a neuron
// ------------------------------------------------------------------------------------------
float SNN::Neuron_Potential(int in, double t, bool delete_history)
{
    int ilayer = Neuron_layer[in];
    float P0 = 0.;
    double t0 = History_time[in][0];
    double delta_t = t - t0;
    if (t0 > 0. && delta_t >= 0. && delta_t < MaxDeltaT)
    {
        int ilayer = Neuron_layer[in];
        P0 = Spike_potential(delta_t, ilayer); // the first event in the history sequence is a spike
    }
    float P = P0;

    // Now we extrapolate the effect of all past spikes and inhibitions to time t, to compute the potential when EPSP arrives
    int len = History_time[in].size();
    if (len > 1)
    {
        for (int ih = 1; ih < len; ih++)
        {
            delta_t = t - History_time[in][ih];
            if (History_type[in][ih] == EPSP)
            { // EPSP
                if (delta_t < MaxDeltaT && (History_time[in][ih] - t0) > tau_r)
                {
                    if (!Void_weight[in][History_ID[in][ih]]) // for type 1 or 3 signals, ID is the stream
                        P += Weight[in][History_ID[in][ih]] * EPS_potential(delta_t);
                }
                else if (delete_history)
                {
                    History_time[in].erase(History_time[in].begin() + ih, History_time[in].begin() + ih + 1);
                    History_type[in].erase(History_type[in].begin() + ih, History_type[in].begin() + ih + 1);
                    History_ID[in].erase(History_ID[in].begin() + ih, History_ID[in].begin() + ih + 1);
                    History_class[in].erase(History_class[in].begin() + ih, History_class[in].begin() + ih + 1);
                    len = len - 1;
                }
            }
            else if (History_type[in][ih] == IPSP)
            { // IPSP
                if (delta_t < MaxDeltaT && EnableIPSP[in][History_ID[in][ih] - N_InputStreams])
                {
                    int ilayer = Neuron_layer[in];
                    P += Inhibitory_potential(delta_t, ilayer);
                }
                else if (delete_history)
                {
                    // get rid of irrelevant events
                    History_time[in].erase(History_time[in].begin() + ih, History_time[in].begin() + ih + 1);
                    History_type[in].erase(History_type[in].begin() + ih, History_type[in].begin() + ih + 1);
                    History_ID[in].erase(History_ID[in].begin() + ih, History_ID[in].begin() + ih + 1);
                    History_class[in].erase(History_class[in].begin() + ih, History_class[in].begin() + ih + 1);
                    
                    len = len - 1;
                }
            }
        }
    }
    return P;
}

// Model Inhibitory-Excitatory signal (IE) as combination of two EPSP-like shapes, a negative one followed by a positive one
// This is a crude model, loosely inspired by shapes in Fig.2 of F. Sandin and M. Nilsson, "Synaptic Delays for Insect-Inspired
// Feature Detection in Dynamic Neuromorphic Processors", doi.org/10.3389/fnins.2020.00150
// ----------------------------------------------------------------------------------------------------------------------------
float SNN::IE_potential(double delta_t, int in, int is)
{
    float sp = 0.;
    if (delta_t >= 0. && delta_t < Delay[in][is])
    {
        sp = -IE_Pot_const * EPS_potential(delta_t);
    }
    else if (delta_t >= Delay[in][is] && delta_t < MaxDeltaT + Delay[in][is])
    {
        delta_t = delta_t - Delay[in][is];
        sp = IE_Pot_const * EPS_potential(delta_t); // So for zero delay, this is an EPSP
    }
    return sp;
}

//LTP rule for weights
void SNN::LTP(int in, double fire_time, bool nearest_spike_approx, SNN &old)
{
    for (int is = 0; is < N_streams; is++)
    {            
        if (Void_weight[in][is])
            continue;
        check_LTD[in][is] = true;
        // Use nearest-spike approximation: search for closest pre-spike
        bool no_prespikes = true;
        int isp = History_time[in].size() - 1;
        float delta_weight = 0;
        
        do
        {
            double delta_t = History_time[in][isp] - fire_time;
            if (History_ID[in][isp] == is && History_type[in][isp]==EPSP && delta_t <= 0) 
            {
                Weight[in][is] += a_plus * exp(delta_t / tau_plus);
                
                if (Weight[in][is] > 1.)
                    Weight[in][is] = 1.;

                no_prespikes = false;
                delta_weight += Weight[in][is] - old.Weight[in][is];

                // in this approximation we're interested only in the first spike
                if (nearest_spike_approx)
                    break;
            }
            isp--;
        } while (isp >= 0 && History_time[in][isp] > fire_time - 7. * tau_plus);
        
        do
        {   
            double delta_t = History_time[in][isp] - fire_time;
            if (History_ID[in][isp] == is && History_type[in][isp]==EPSP && delta_t <= 0) 
            {
                if(is < N_InputStreams){
                    Delay[in][is]  += d_plus * exp(delta_t / taud_plus);
                
                    if (Delay[in][is] > MaxDelay)
                        Delay[in][is] = MaxDelay;
                }          
                no_prespikes = false;

                // in this approximation we're interested only in the first spike
                if (nearest_spike_approx)
                    break;
            }
            isp--;
        } while (isp >= 0 && History_time[in][isp] > fire_time - 7. * taud_plus);

        if (!no_prespikes)
            Renorm(in, old);
    }
    return;
}

void SNN::LTD(int in, int is, double spike_time, bool nearest_spike_approx, SNN &old)
{
    return;
    if(!check_LTD[in][is]) return;
    if (Fire_time[in].size() == 0){
        //cout <<"skip" << endl; 
        return;
    }
    if (Void_weight[in][is]){
        //cout << "void skip" <<endl;
        return;
    }
    
    double delta_t = spike_time - Fire_time[in].back();
    // if nearest_spike_approx we prevent to compute future LTD until the next activation of the neuron
    if (nearest_spike_approx)
        check_LTD[in][is] = false;

    //cout << "Go " << delta_t << " " << spike_time<< " " << Fire_time[in].back() << endl;
    if (delta_t >= 0 && delta_t < 7. * tau_minus)
    {
        Weight[in][is] -= a_minus * exp(-delta_t / tau_minus);
                    
        if (Weight[in][is] < 0.)
            Weight[in][is] = 0.;
     
    }
    if (delta_t >= 0 && delta_t < 7. * taud_minus)
    {
        if(is<N_InputStreams){
            Delay[in][is] -= d_minus * exp(-delta_t / taud_minus);
            
            if (Delay[in][is] < 0.)
                Delay[in][is] = 0.; 
        }                    
    }
    Renorm(in, old);
    if(delta_t >= 0 && delta_t < 7. * max(taud_minus, tau_minus)){
        Fire_time[in].clear();
    }
    return;
}

//different approach: compute ltd bewtween two consecutive spikes
void SNN::Compute_LTD(int in, double fire_time, bool nearest_spike_approx, SNN &old)
{
    if (Fire_time[in].size() == 0){
        //cout <<"skip" << endl; 
        return;
    }

    bool no_prespikes = true;
    int isp = 1;
    double previous_firetime = Fire_time[in].back();

    for (int is = 0; is < N_streams; is++)
    {
        if (Void_weight[in][is])
            continue;
        
        do
        {
            double delta_t = History_time[in][isp] - previous_firetime;
            if (History_ID[in][isp] == is && History_type[in][isp]==EPSP && delta_t >= 0) 
            {
                Weight[in][is] -= a_minus * exp(-delta_t / tau_minus);
                
                if (Weight[in][is] < 0.)
                    Weight[in][is] = 0.;

                no_prespikes = false;
                
                // in this approximation we're interested only in the first spike
                if (nearest_spike_approx)
                    break;
            }
            isp++;
        } while (isp < History_time[in].size() && History_time[in][isp] < previous_firetime + 7. * tau_minus && History_time[in][isp]<fire_time);
        
        do
        {
            double delta_t = History_time[in][isp] - previous_firetime;
            if (History_ID[in][isp] == is && History_type[in][isp]==EPSP && delta_t >= 0) 
            {
                if(is<N_InputStreams){
                    Delay[in][is] -= d_minus * exp(-delta_t / taud_minus);
                    
                    if (Delay[in][is] < 0.)
                        Delay[in][is] = 0.; 
                }
                no_prespikes = false;
                
                // in this approximation we're interested only in the first spike
                if (nearest_spike_approx)
                    break;
            }
            isp++;
        } while (isp < History_time[in].size() && History_time[in][isp] < previous_firetime + 7. * taud_minus && History_time[in][isp] < fire_time);
        if(!no_prespikes)
            Renorm(in, old);        
    }
    return;
}

void SNN::Renorm(int in, SNN &old) {
    float weight_sum = 0.0;
    double delay_factor = 0.;

    // Calculate the sum of weights for the 'in' neuron
    for (int is = 0; is < N_streams; is++) {
        if(!Void_weight[in][is]) weight_sum += Weight[in][is];
        delay_factor+=Delay[in][is];
    }
    //cout << "delay factor " << in << "  " << delay_factor << endl;
    delay_factor /= sumdelays[in];
    // Check if the sum is greater than 0 to avoid division by zero
    if (weight_sum > 0.0) {
        // Normalize the weights for the 'in' neuron
        for (int is = 0; is < N_streams; is++) {
            if(!Void_weight[in][is]){
                Weight[in][is] = Weight[in][is]/weight_sum;
                old.Weight[in][is] = Weight[in][is];
           }
           Delay[in][is]/=delay_factor;
           //TODO: maybe I should recheck the boundary conditions for the delays.
        }
    }

    return;    
}


void SNN::Renorm_Opt(int in, float delta_weight, SNN &old)
{
    float norm_factor = 1. + delta_weight;
    for (int is = 0; is < N_streams; is++)
    {
        if (!Void_weight[in][is])
        {
            Weight[in][is] /= norm_factor;
            old.Weight[in][is] = Weight[in][is];
        }
    }
    return;
}

void SNN::PrintWeights(){
    for(int in = 0; in<N_neurons; in++){
        cout << "Neuron " << in << endl;
        for (int is = 0; is < N_streams; is++)
        {
            if (!Void_weight[in][is]){
            cout << Weight[in][is] << ", ";
            }
        }
        cout << endl <<endl;
    }
}

void SNN::PrintSNN(){
    // Print all the variables
    cout << "alpha = " << alpha << endl;
    cout << "CFI0 = " << CFI0 << endl;
    cout << "CFI1 = " << CFI1 << endl;
    cout << "CF01 = " << CF01 << endl;
    cout << "L1inhibitfactor = " << L1inhibitfactor << endl;
    cout << "K = " << K << endl;
    cout << "K1 = " << K1 << endl;
    cout << "K2 = " << K2 << endl;
    cout << "IE_Pot_const = " << IE_Pot_const << endl;
    cout << "IPSP_dt_dilation = " << IPSP_dt_dilation << endl;
    cout << "MaxDelay = " << MaxDelay << endl;
    cout << "tau_m = " << tau_m << endl;
    cout << "tau_s = " << tau_s << endl;
    cout << "tau_r = " << tau_r << endl;
    cout << "tau_plus = " << tau_plus << endl;
    cout << "tau_minus = " << tau_minus << endl;
    cout << "a_plus = " << a_plus << endl;
    cout << "a_minus = " << a_minus << endl;
    cout << "taud_plus = " << taud_plus << endl;
    cout << "taud_minus = " << taud_minus << endl;
    cout << "d_plus = " << d_plus << endl;
    cout << "d_minus = " << d_minus << endl;
    cout << "N_neurons = " << N_neurons << endl;
    cout << "N_streams = " << N_streams << endl;
    cout << "Threshold[0] = " << Threshold[0] << endl;
    cout << "Threshold[1] = " << Threshold[1] << endl;
    cout << "tmax = " << tmax << endl;
    cout << "MaxDeltaT = " << MaxDeltaT << endl;
    cout << "N_neuronsL[0] = " << N_neuronsL[0] << endl;
    cout << "N_neuronsL[1] = " << N_neuronsL[1] << endl;
    cout << "N_InputStreams = " << N_InputStreams << endl;
    cout << "largenumber = " << largenumber << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "sparsity = "<< sparsity << endl;
    cout << "split layer0 = " << split_layer0 << endl;
    cout << "-------------------------------------" << endl;
}