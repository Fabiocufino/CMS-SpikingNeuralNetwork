#include "SNN.h"
#include <iostream>
#include <cstring>

using json = nlohmann::json;
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
    tmax = tau_s * tau_m / (tau_m - tau_s) * log(tau_m/tau_s);

    MaxDeltaT = 7. * tau_m;

    fire_granularity = tau_s / 5.;
    fire_precision = min(Threshold[0], Threshold[1]) *2.5 / 100.;
    myRNG = new TRandom3(static_cast<unsigned int>(time(0)));
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
    History_ev_class = new vector<pair <int, int>>[N_neurons];  // Class of the signal

    Fire_time = new vector<double>[N_neurons];    // Times of firing of each neuron
    Neuron_layer = new int[N_neurons];
    sumweight = new float[N_neurons]; // summed weights of streams for each neurons for the purpose of normalization
    sumdelays = new double[N_neurons]; // summed delays of streams for each neurons for the purpose of normalization


    Delta_delay = MaxDelay/10.;
    Mean_delay = MaxDelay/2.;

    Init_neurons(0);
    Init_connection_map();
    Init_weights();
    Init_delays_uniform();
}
SNN::SNN() : SNN(1, 1, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0f, 0.0f, 0.0f, false) {}
SNN::~SNN() {
    // Release memory for dynamically allocated arrays
    for (int i = 0; i < N_neurons; ++i) {
        delete[] Weight[i];
        delete[] Weight_initial[i];
        delete[] check_LTD[i];
        delete[] Void_weight[i];
        delete[] Delay[i];
        delete[] Delay_initial[i];
        delete[] EnableIPSP[i];
    }

    delete[] Weight;
    delete[] Weight_initial;
    delete[] check_LTD;
    delete[] Void_weight;
    delete[] Delay;
    delete[] Delay_initial;
    delete[] EnableIPSP;
    
    // Clear vectors
    delete[] History_time;
    delete[] History_type;
    delete[] History_ID;
    delete[] History_ev_class;
    delete[] Fire_time;

    delete[] Neuron_layer;
    delete[] sumweight;
    delete[] sumdelays;

    // Delete random number generator
    delete myRNG;
}

void SNN::Reset_weights(){
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_InputStreams; is++)
            Weight[in][is] = Weight_initial[in][is];
    }
}

double SNN::bisectionMethod(double a, double b, int in, double epsilon, function<float(int, double, bool)> func)
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
        if (abs(fc) < epsilon)
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
void SNN::Init_neurons(int ievent)
{

    for (int in = 0; in < N_neurons; in++)
    {
        // Set first event in history of this neuron
        History_time[in].clear();
        History_type[in].clear();
        History_ID[in].clear();
        History_ev_class[in].clear();

        History_time[in].push_back(0);
        History_type[in].push_back(SPIKE);
        History_ID[in].push_back(0);
        History_ev_class[in].push_back({ievent, NOCLASS});

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
void SNN::insert_spike(int id_neuron, double spike_time, int type, int id, int spike_class, int ievent){
    // Find the position where the new value should be inserted
    auto it = lower_bound(History_time[id_neuron].begin(), History_time[id_neuron].end(), spike_time);
    int position = distance(History_time[id_neuron].begin(), it);
    // Insert the new value at the found position
    History_time[id_neuron].insert(it, spike_time);
    History_type[id_neuron].insert(History_type[id_neuron].begin()+position, type);
    History_ID[id_neuron].insert(History_ID[id_neuron].begin()+position, id);
    History_ev_class[id_neuron].insert(History_ev_class[id_neuron].begin()+position, {ievent, spike_class});
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

//function to extract the unique values in a subset of a vector
vector<pair<int, int>> uniquePairsInRange(const vector<pair<int, int>>& my_pairs, int start_idx, int end_idx) {
    if (start_idx < 0 || end_idx >= my_pairs.size() || start_idx > end_idx) {
        cout << "start_idx: " << start_idx << " end_idx: " << end_idx << " my_pairs.size(): " << my_pairs.size() << endl;
        throw std::out_of_range("Index out of bounds in uniquePairsInRange");
        //print start_idx, end_idx e my_pairs.size()
        
    }
    vector<pair<int, int>> range_subset(my_pairs.begin() + start_idx, my_pairs.begin() + end_idx + 1);

    // Sort the subrange
    sort(range_subset.begin(), range_subset.end());

    // Remove consecutive duplicates
    auto it = unique(range_subset.begin(), range_subset.end(), [](const auto& a, const auto& b) {
        return a.first == b.first && a.second == b.second;
    });

    // Resize the vector to remove the duplicates
    range_subset.resize(distance(range_subset.begin(), it));

    return range_subset;
}

//Inspect the history of a neuron before the activation
vector<pair <int, int>> SNN::Inspect_History(int in, double fire_time, double window){
    //find the position of the fire_time and the first spike in the window
    int start_pos = distance(History_time[in].begin(), lower_bound(History_time[in].begin(), History_time[in].end(), fire_time - window));
    int end_pos   = distance(History_time[in].begin(), upper_bound(History_time[in].begin(), History_time[in].end(), fire_time));
    if (end_pos >= History_ev_class[in].size()) {
        end_pos = History_ev_class[in].size() - 1;
    }
    return (uniquePairsInRange(History_ev_class[in], start_pos, end_pos));
    
}

//Handle the activation of a neuron
void SNN::Activate_Neuron(int in, double t) {
    // Reset history of this neuron before time t
    auto it = upper_bound(History_time[in].begin(), History_time[in].end(), t);
    if (it != History_time[in].begin()) {
        int position = distance(History_time[in].begin(), it);
       
        // Erase elements before the position found by binary search
        History_time[in].erase(History_time[in].begin(), it);
        History_type[in].erase(History_type[in].begin(), History_type[in].begin() + position);
        History_ID[in].erase(History_ID[in].begin(), History_ID[in].begin() + position);
        History_ev_class[in].erase(History_ev_class[in].begin(), History_ev_class[in].begin() + position);
    }
    
    // Insert the spike at time t
    if (!History_ev_class[in].empty()) {
        insert_spike(in, t, SPIKE, 0, NOCLASS, History_ev_class[in].front().first);
    } else {
        // Handle the case when the history is empty
        insert_spike(in, t, SPIKE, 0, NOCLASS, 0); // Assuming default value for last event class
    }
    return;
    /*
    History_time[in].clear();
    History_type[in].clear();
    History_ID[in].clear();
    */
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
        double last_spike = History_time[in].front();
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
                else if (delete_history && (History_time[in][ih] - last_spike) > 7. * max(tau_minus, taud_minus) && delta_t > 7. * max(tau_plus, taud_plus))
                {
                    History_time[in].erase(History_time[in].begin() + ih, History_time[in].begin() + ih + 1);
                    History_type[in].erase(History_type[in].begin() + ih, History_type[in].begin() + ih + 1);
                    History_ID[in].erase(History_ID[in].begin() + ih, History_ID[in].begin() + ih + 1);
                    History_ev_class[in].erase(History_ev_class[in].begin() + ih, History_ev_class[in].begin() + ih + 1);
                    
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
                    History_ev_class[in].erase(History_ev_class[in].begin() + ih, History_ev_class[in].begin() + ih + 1);
                    
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
void SNN::LTP_weights(int in, double fire_time, bool nearest_spike_approx, SNN &old)
{
    for (int is = 0; is < N_streams; is++)
    {            
        if (Void_weight[in][is])
            continue;
        
        check_LTD[in][is] = true;
        bool no_prespikes = true;
        int isp = History_time[in].size() - 1;

        while (isp >= 0 && History_time[in][isp] > fire_time - 7. * tau_plus)
        {
            double delta_t = History_time[in][isp] - fire_time;
            if (History_ID[in][isp] == is && History_type[in][isp] == EPSP && delta_t <= 0) 
            {
                Weight[in][is] += a_plus * exp(delta_t / tau_plus);
                if (Weight[in][is] > 1.)
                    Weight[in][is] = 1.;

                no_prespikes = false;
    
                if (nearest_spike_approx)
                    break;
            }
            isp--;
        }

        // Decomment if you want to renormalize the weights
        //if (!no_prespikes)
        //    Renorm_weights(in, old);
    }
    return;
}

void SNN::LTP_delays(int in, double fire_time, bool nearest_spike_approx, SNN &old)
{
    for (int is = 0; is < N_streams; is++)
    {            
        if (Void_weight[in][is])
            continue;
        
        bool no_prespikes = true;
        int isp = History_time[in].size() - 1;
        
        while (isp >= 0 && History_time[in][isp] > fire_time - 7. * taud_plus)
        {
            double delta_t = History_time[in][isp] - fire_time;
            if (History_ID[in][isp] == is && History_type[in][isp] == EPSP && delta_t < 0) 
            {
                if (is < N_InputStreams)
                {
                    if(delta_t < -tmax){
                        Delay[in][is] += d_plus * exp(delta_t / taud_plus);
                        if (Delay[in][is] > MaxDelay)
                            Delay[in][is] = MaxDelay;
                    }
                    else{
                        Delay[in][is] -= d_minus * exp(delta_t / taud_minus); 
                        if (Delay[in][is] < 0.)
                            Delay[in][is] = 0.;
                    } 
                }
                no_prespikes = false;

                if (nearest_spike_approx)
                    break;
            }
            isp--;
        }

        if (!no_prespikes)
            Renorm_delays(in, old);
    }
    return;
}

void SNN::LTD_weights(int in, double fire_time, bool nearest_spike_approx, SNN &old)
{
    if (Fire_time[in].empty()) return;

    int isp = 1;
    double previous_firetime = Fire_time[in].back();

    for (int is = 0; is < N_streams; is++)
    {
        if (Void_weight[in][is])
            continue;
        
        bool no_prespikes = true; 
        isp = 1; 

        while (isp < History_time[in].size() && History_time[in][isp] < previous_firetime + 7 * tau_minus && History_time[in][isp] < fire_time)
        {
            double delta_t = History_time[in][isp] - previous_firetime;
            if (History_ID[in][isp] == is && History_type[in][isp] == EPSP && delta_t >= 0) 
            {
                Weight[in][is] -= a_minus * exp(-delta_t / tau_minus);
                
                if (Weight[in][is] < 0.)
                    Weight[in][is] = 0.;

                no_prespikes = false;
                
                if (nearest_spike_approx)
                    break;
            }
            isp++;
        }

        // Decomment if you want to renormalize the weights
        //if (!no_prespikes)
        //    Renorm_weights(in, old);
    }
    return;
}

void SNN::LTD_delays(int in, double fire_time, bool nearest_spike_approx, SNN &old)
{
    if (Fire_time[in].empty()) return;

    int isp = 1;
    double previous_firetime = Fire_time[in].back();

    for (int is = 0; is < N_streams; is++)
    {
        if (Void_weight[in][is])
            continue;
        
        bool no_prespikes = true;  
        isp = 1;

        while (isp < History_time[in].size() && History_time[in][isp] < previous_firetime + 7 * taud_minus && History_time[in][isp] < fire_time)
        { 
            double delta_t = History_time[in][isp] - previous_firetime;
            if (History_ID[in][isp] == is && History_type[in][isp] == EPSP && delta_t > 0) 
            {
                if (is < N_InputStreams)
                {
                    Delay[in][is] -= d_minus * exp(- delta_t / taud_minus); 
                    if (Delay[in][is] < 0.)
                        Delay[in][is] = 0.;
                }
                no_prespikes = false;
                
                if (nearest_spike_approx)
                    break;
            }
            isp++;
        }
        
        if (!no_prespikes)
            Renorm_delays(in, old);
    }
    return;
}

void SNN::Renorm_weights(int in, SNN &old) {
    if (in < 0 || in >= N_neurons) {
        std::cerr << "Error: Index 'in' out of bounds in Renorm. Aborting operation." << std::endl;
        return;
    }

    if (N_streams <= 0) {
        std::cerr << "Error: N_streams is non-positive in Renorm. Aborting operation." << std::endl;
        return;
    }

    if (!Void_weight || !Weight || !old.Weight) {
        std::cerr << "Error: Null pointer encountered in Renorm. Aborting operation." << std::endl;
        return;
    }

    float weight_sum = 0.0;
    // Calculate the sum of weights for the 'in' neuron
    for (int is = 0; is < N_streams; is++) {
        if (is < N_neurons && !Void_weight[in][is]) {
            weight_sum += Weight[in][is];
        }
    }

    // Check if the sum is greater than 0 to avoid division by zero
    if (weight_sum > 0.0) {
        // Normalize the weights for the 'in' neuron
        for (int is = 0; is < N_streams; is++) {
            if (is < N_neurons && !Void_weight[in][is]) {
                Weight[in][is] =  Weight[in][is]/weight_sum;
                old.Weight[in][is] = Weight[in][is];
            }
        }
    }
}

void SNN::Renorm_delays(int in, SNN &old) {
    if (in < 0 || in >= N_neurons) {
        std::cerr << "Error: Index 'in' out of bounds in Renorm. Aborting operation." << std::endl;
        return;
    }

    if (N_streams <= 0) {
        std::cerr << "Error: N_streams is non-positive in Renorm. Aborting operation." << std::endl;
        return;
    }

    if (!Delay || !sumdelays) {
        std::cerr << "Error: Null pointer encountered in Renorm. Aborting operation." << std::endl;
        return;
    }

    double delay_factor = 0.0;

    for (int is = 0; is < N_streams; is++) {
        if (is < N_neurons && !Void_weight[in][is]) {
            delay_factor += Delay[in][is];
        }
    }

    if (sumdelays[in] == 0.0) {
        std::cerr << "Error: sumdelays[" << in << "] is zero. Aborting operation to avoid division by zero." << std::endl;
        return;
    }

    delay_factor /= sumdelays[in];

    // Check if the sum is greater than 0 to avoid division by zero
    if (delay_factor > 0.0) {
        // Normalize the weights for the 'in' neuron
        for (int is = 0; is < N_streams; is++) {
            if (is < N_neurons) {
                Delay[in][is] /= delay_factor;
            }
        }
    }
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

void SNN::PrintDelays(){
    for(int in = 0; in<N_neurons; in++){
        cout << "Neuron " << in << endl;
        for (int is = 0; is < N_streams; is++)
        {
            if (!Void_weight[in][is]){
            cout << Delay[in][is] << ", ";
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

void SNN::copy_from(const SNN& other) {
    // Deallocate previously allocated memory
    for (int i = 0; i < N_neurons; ++i) {
        delete[] Weight[i];
        delete[] Weight_initial[i];
        delete[] check_LTD[i];
        delete[] Void_weight[i];
        delete[] Delay[i];
        delete[] Delay_initial[i];
        delete[] EnableIPSP[i];
    }
    delete[] Weight;
    delete[] Weight_initial;
    delete[] check_LTD;
    delete[] Void_weight;
    delete[] Delay;
    delete[] Delay_initial;
    delete[] EnableIPSP;

    // Deallocate and reallocate vector arrays
    delete[] History_time;
    delete[] History_type;
    delete[] History_ID;
    delete[] History_ev_class;
    delete[] Fire_time;

    delete[] sumdelays;
    delete[] Neuron_layer;
    delete[] sumweight;

    // Copy scalar values
    alpha = other.alpha;
    CFI0 = other.CFI0;
    CFI1 = other.CFI1;
    CF01 = other.CF01;
    L1inhibitfactor = other.L1inhibitfactor;
    K = other.K;
    K1 = other.K1;
    K2 = other.K2;
    IE_Pot_const = other.IE_Pot_const;
    IPSP_dt_dilation = other.IPSP_dt_dilation;
    MaxDelay = other.MaxDelay;
    tau_m = other.tau_m;
    tau_s = other.tau_s;
    tau_r = other.tau_r;
    tau_plus = other.tau_plus;
    tau_minus = other.tau_minus;
    taud_plus = other.taud_plus;
    taud_minus = other.taud_minus;
    a_plus = other.a_plus;
    a_minus = other.a_minus;
    d_plus = other.d_plus;
    d_minus = other.d_minus;
    sparsity = other.sparsity;
    split_layer0 = other.split_layer0;
    N_InputStreams = other.N_InputStreams;
    N_streams = other.N_streams;
    fire_granularity = other.fire_granularity;
    fire_precision = other.fire_precision;
    Delta_delay = other.Delta_delay;
    Mean_delay = other.Mean_delay;
    N_neurons = other.N_neurons;
    Threshold[0] = other.Threshold[0];
    Threshold[1] = other.Threshold[1];
    tmax = other.tmax;
    MaxDeltaT = other.MaxDeltaT;

    // Allocate memory for new arrays
    Weight = new float*[N_neurons];
    Weight_initial = new float*[N_neurons];
    check_LTD = new bool*[N_neurons];
    Void_weight = new bool*[N_neurons];
    Delay = new double*[N_neurons];
    Delay_initial = new double*[N_neurons];
    EnableIPSP = new bool*[N_neurons];

    for (int i = 0; i < N_neurons; ++i) {
        Weight[i] = new float[N_streams];
        Weight_initial[i] = new float[N_streams];
        check_LTD[i] = new bool[N_streams];
        Void_weight[i] = new bool[N_streams];
        Delay[i] = new double[N_streams];
        Delay_initial[i] = new double[N_streams];
        EnableIPSP[i] = new bool[N_neurons];

        // Manually copy elements in a for loop
        for (int j = 0; j < N_streams; ++j) {
            Weight[i][j] = other.Weight[i][j];
            Weight_initial[i][j] = other.Weight_initial[i][j];
            check_LTD[i][j] = other.check_LTD[i][j];
            Void_weight[i][j] = other.Void_weight[i][j];
            Delay[i][j] = other.Delay[i][j];
            Delay_initial[i][j] = other.Delay_initial[i][j];
        }
        for (int j = 0; j< N_neurons; ++j)
            EnableIPSP[i][j] = other.EnableIPSP[i][j];
    }


    History_time = new std::vector<double>[N_neurons];
    History_type = new std::vector<int>[N_neurons];
    History_ID = new std::vector<int>[N_neurons];
    History_ev_class = new std::vector<std::pair<int, int>>[N_neurons];
    Fire_time = new std::vector<double>[N_neurons];

    for (int i = 0; i < N_neurons; ++i) {
        History_time[i] = other.History_time[i];
        History_type[i] = other.History_type[i];
        History_ID[i] = other.History_ID[i];
        History_ev_class[i] = other.History_ev_class[i];
        Fire_time[i] = other.Fire_time[i];
    }

    // Copy the Neuron_layer array
    Neuron_layer = new int[N_neurons];
    for (int i = 0; i < N_neurons; ++i) {
        Neuron_layer[i] = other.Neuron_layer[i];
    }

    // Copy the sumweight array
    sumweight = new float[N_neurons];
    for (int i = 0; i < N_neurons; ++i) {
        sumweight[i] = other.sumweight[i];
    }

    // Copy the sumdelays array
    sumdelays = new double[N_neurons];
    for (int i = 0; i < N_neurons; ++i) {
        sumdelays[i] = other.sumdelays[i];
    }

    // Copy the N_neuronsL array
    N_neuronsL[0] = other.N_neuronsL[0];
    N_neuronsL[1] = other.N_neuronsL[1];

    // Copy the random number generator
    if (myRNG) {
        delete myRNG;
    }
    myRNG = new TRandom3(*other.myRNG);

    // Copy the remaining scalar values
    largenumber = other.largenumber;
    epsilon = other.epsilon;
}

void SNN::dumpToJson(const string& filename) {
    json j;

    j["alpha"] = alpha;
    j["CFI0"] = CFI0;
    j["CFI1"] = CFI1;
    j["CF01"] = CF01;
    j["L1inhibitfactor"] = L1inhibitfactor;
    j["K"] = K;
    j["K1"] = K1;
    j["K2"] = K2;
    j["IE_Pot_const"] = IE_Pot_const;
    j["IPSP_dt_dilation"] = IPSP_dt_dilation;
    j["MaxDelay"] = MaxDelay;
    j["tau_m"] = tau_m;
    j["tau_s"] = tau_s;
    j["tau_r"] = tau_r;
    j["tau_plus"] = tau_plus;
    j["tau_minus"] = tau_minus;
    j["taud_plus"] = taud_plus;
    j["taud_minus"] = taud_minus;
    j["a_plus"] = a_plus;
    j["a_minus"] = a_minus;
    j["d_plus"] = d_plus;
    j["d_minus"] = d_minus;
    j["sparsity"] = sparsity;
    j["split_layer0"] = split_layer0;
    j["N_InputStreams"] = N_InputStreams;
    j["N_streams"] = N_streams;
    j["fire_granularity"] = fire_granularity;
    j["fire_precision"] = fire_precision;
    j["Delta_delay"] = Delta_delay;
    j["Mean_delay"] = Mean_delay;
    j["N_neurons"] = N_neurons;
    j["Threshold"] = {Threshold[0], Threshold[1]};
    j["tmax"] = tmax;
    j["MaxDeltaT"] = MaxDeltaT;

    json weights;
    json weight_initials;
    json check_ltds;
    json void_weights;
    json delays;
    json delay_initials;
    json enable_ipsp;

    for (int i = 0; i < N_neurons; ++i) {
        json weight_row;
        json weight_initial_row;
        json check_ltd_row;
        json void_weight_row;
        json delay_row;
        json delay_initial_row;
        json enable_ipsp_row;

        for (int j = 0; j < N_streams; ++j) {
            weight_row.push_back(Weight[i][j]);
            weight_initial_row.push_back(Weight_initial[i][j]);
            check_ltd_row.push_back(check_LTD[i][j]);
            void_weight_row.push_back(Void_weight[i][j]);
            delay_row.push_back(Delay[i][j]);
            delay_initial_row.push_back(Delay_initial[i][j]);
        }

        for (int j = 0; j < N_neurons; j++)
            enable_ipsp_row.push_back(EnableIPSP[i][j]);

        weights[to_string(i)] = weight_row;
        weight_initials[to_string(i)] = weight_initial_row;
        check_ltds[to_string(i)] = check_ltd_row;
        void_weights[to_string(i)] = void_weight_row;
        delays[to_string(i)] = delay_row;
        delay_initials[to_string(i)] = delay_initial_row;
        enable_ipsp[to_string(i)] = enable_ipsp_row;
    }

    j["Weight"] = weights;
    j["Weight_initial"] = weight_initials;
    j["check_LTD"] = check_ltds;
    j["Void_weight"] = void_weights;
    j["Delay"] = delays;
    j["Delay_initial"] = delay_initials;
    j["EnableIPSP"] = enable_ipsp;

    j["History_time"] = json::array();
    j["History_type"] = json::array();
    j["History_ID"] = json::array();
    j["History_ev_class"] = json::array();
    j["Fire_time"] = json::array();

    for (int i = 0; i < N_neurons; ++i) {
        j["History_time"].push_back(History_time[i]);
        j["History_type"].push_back(History_type[i]);
        j["History_ID"].push_back(History_ID[i]);
        j["History_ev_class"].push_back(History_ev_class[i]);
        j["Fire_time"].push_back(Fire_time[i]);
    }

    j["Neuron_layer"] = json::array();
    for (int i = 0; i < N_neurons; ++i) {
        j["Neuron_layer"].push_back(Neuron_layer[i]);
    }

    j["sumweight"] = json::array();
    for (int i = 0; i < N_neurons; ++i) {
        j["sumweight"].push_back(sumweight[i]);
    }

    j["sumdelays"] = json::array();
    for (int i = 0; i < N_neurons; ++i) {
        j["sumdelays"].push_back(sumdelays[i]);
    }

    j["N_neuronsL"] = {N_neuronsL[0], N_neuronsL[1]};
    j["myRNG"] = myRNG->GetSeed();
    j["largenumber"] = largenumber;
    j["epsilon"] = epsilon;

    ofstream file(filename);
    file << j.dump(4);
}

void SNN::loadFromJson(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open file");
    }

    json j;
    file >> j;

    alpha = j["alpha"].get<float>();
    CFI0 = j["CFI0"].get<float>();
    CFI1 = j["CFI1"].get<float>();
    CF01 = j["CF01"].get<float>();
    L1inhibitfactor = j["L1inhibitfactor"].get<float>();
    K = j["K"].get<float>();
    K1 = j["K1"].get<float>();
    K2 = j["K2"].get<float>();
    IE_Pot_const = j["IE_Pot_const"].get<float>();
    IPSP_dt_dilation = j["IPSP_dt_dilation"].get<double>();
    MaxDelay = j["MaxDelay"].get<double>();
    tau_m = j["tau_m"].get<double>();
    tau_s = j["tau_s"].get<double>();
    tau_r = j["tau_r"].get<double>();
    tau_plus = j["tau_plus"].get<double>();
    tau_minus = j["tau_minus"].get<double>();
    taud_plus = j["taud_plus"].get<double>();
    taud_minus = j["taud_minus"].get<double>();
    a_plus = j["a_plus"].get<double>();
    a_minus = j["a_minus"].get<double>();
    d_plus = j["d_plus"].get<double>();
    d_minus = j["d_minus"].get<double>();
    sparsity = j["sparsity"].get<float>();
    split_layer0 = j["split_layer0"].get<bool>();
    N_InputStreams = j["N_InputStreams"].get<int>();
    N_streams = j["N_streams"].get<int>();
    fire_granularity = j["fire_granularity"].get<double>();
    fire_precision = j["fire_precision"].get<float>();
    Delta_delay = j["Delta_delay"].get<double>();
    Mean_delay = j["Mean_delay"].get<double>();
    N_neurons = j["N_neurons"].get<int>();
    Threshold[0] = j["Threshold"][0].get<float>();
    Threshold[1] = j["Threshold"][1].get<float>();
    tmax = j["tmax"].get<double>();
    MaxDeltaT = j["MaxDeltaT"].get<double>();

    // Allocate memory for arrays
    Weight = new float*[N_neurons];
    Weight_initial = new float*[N_neurons];
    check_LTD = new bool*[N_neurons];
    Void_weight = new bool*[N_neurons];
    Delay = new double*[N_neurons];
    Delay_initial = new double*[N_neurons];
    EnableIPSP = new bool*[N_neurons];

    for (int i = 0; i < N_neurons; ++i) {
        Weight[i] = new float[N_streams];
        Weight_initial[i] = new float[N_streams];
        check_LTD[i] = new bool[N_streams];
        Void_weight[i] = new bool[N_streams];
        Delay[i] = new double[N_streams];
        Delay_initial[i] = new double[N_streams];
        EnableIPSP[i] = new bool[N_neurons];
    }

    json weights            = j["Weight"];
    json weight_initials    = j["Weight_initial"];
    json check_ltds         = j["check_LTD"];
    json void_weights       = j["Void_weight"];
    json delays             = j["Delay"];
    json delay_initials     = j["Delay_initial"];
    json enable_ipsp        = j["EnableIPSP"];

    // Copy data for Weight_initial matrix
    for (json::iterator it = weights.begin(); it != weights.end(); ++it) {
        int neuron_index = stoi(it.key());
        json weights_row = it.value();
        for (size_t j = 0; j < weights_row.size(); ++j) {
            Weight[neuron_index][j] = weights_row[j].get<float>();
        }
    }

    // Copy data for Weight_initial matrix
    for (json::iterator it = weight_initials.begin(); it != weight_initials.end(); ++it) {
        int neuron_index = stoi(it.key());
        json weight_initial_row = it.value();
        for (size_t j = 0; j < weight_initial_row.size(); ++j) {
            Weight_initial[neuron_index][j] = weight_initial_row[j].get<float>();
        }
    }

    // Copy data for check_LTD matrix
    for (json::iterator it = check_ltds.begin(); it != check_ltds.end(); ++it) {
        int neuron_index = stoi(it.key());
        json check_ltd_row = it.value();
        for (size_t j = 0; j < check_ltd_row.size(); ++j) {
            check_LTD[neuron_index][j] = check_ltd_row[j].get<bool>();
        }
    }

    // Copy data for Void_weight matrix
    for (json::iterator it = void_weights.begin(); it != void_weights.end(); ++it) {
        int neuron_index = stoi(it.key());
        json void_weight_row = it.value();
        for (size_t j = 0; j < void_weight_row.size(); ++j) {
            Void_weight[neuron_index][j] = void_weight_row[j].get<bool>();
        }
    }

    // Copy data for Delay matrix
    for (json::iterator it = delays.begin(); it != delays.end(); ++it) {
        int neuron_index = stoi(it.key());
        json delay_row = it.value();
        for (size_t j = 0; j < delay_row.size(); ++j) {
            Delay[neuron_index][j] = delay_row[j].get<double>();
        }
    }

    // Copy data for Delay_initial matrix
    for (json::iterator it = delay_initials.begin(); it != delay_initials.end(); ++it) {
        int neuron_index = stoi(it.key());
        json delay_initial_row = it.value();
        for (size_t j = 0; j < delay_initial_row.size(); ++j) {
            Delay_initial[neuron_index][j] = delay_initial_row[j].get<double>();
        }
    }

    // Copy data for EnableIPSP matrix
    for (json::iterator it = enable_ipsp.begin(); it != enable_ipsp.end(); ++it) {
        int neuron_index = stoi(it.key());
        json enable_ipsp_row = it.value();
        for (size_t j = 0; j < enable_ipsp_row.size(); ++j) {
            EnableIPSP[neuron_index][j] = enable_ipsp_row[j].get<bool>();
        }
    }

    History_time = new vector<double>[N_neurons];
    History_type = new vector<int>[N_neurons];
    History_ID = new vector<int>[N_neurons];
    History_ev_class = new vector<pair<int, int>>[N_neurons];
    Fire_time = new vector<double>[N_neurons];

    for (int i = 0; i < N_neurons; ++i) {
        History_time[i] = j["History_time"][i].get<vector<double>>();
        History_type[i] = j["History_type"][i].get<vector<int>>();
        History_ID[i] = j["History_ID"][i].get<vector<int>>();
        History_ev_class[i] = j["History_ev_class"][i].get<vector<pair<int, int>>>();
        Fire_time[i] = j["Fire_time"][i].get<vector<double>>();
    }

    Neuron_layer = new int[N_neurons];
    for (int i = 0; i < N_neurons; ++i) {
        Neuron_layer[i] = j["Neuron_layer"][i].get<int>();
    }

    sumweight = new float[N_neurons];
    for (int i = 0; i < N_neurons; ++i) {
        sumweight[i] = j["sumweight"][i].get<float>();
    }

    sumdelays = new double[N_neurons];
    for (int i = 0; i < N_neurons; ++i) {
        sumdelays[i] = j["sumdelays"][i].get<double>();
    }

    N_neuronsL[0] = j["N_neuronsL"][0].get<int>();
    N_neuronsL[1] = j["N_neuronsL"][1].get<int>();

    if (myRNG) {
        delete myRNG;
    }
    myRNG = new TRandom3(j["myRNG"].get<unsigned int>());

    largenumber = j["largenumber"].get<double>();
    epsilon = j["epsilon"].get<double>();

    cout << "File loaded succesfully" << endl;
}

