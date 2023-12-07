#include "SNN.h"

using namespace std;

SNN::SNN(int _NL0, int _NL1,
         float _alpha,
         float _CFI0, float _CFI1, float _CF01,
         float _L1inhibitfactor,
         float _K, float _K1, float _K2,
         float _IE_Pot_const, float _IPSP_dt_dilation,
         float _MaxDelay,

         float _tau_m, float _tau_s, float _tau_r, float _tau_plus, float _tau_minus,
         float _a_plus, float _a_minus,

         int _N_InputStreams,
         float _Threshold0, float _Threshold1) :
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

                                                 N_InputStreams(_N_InputStreams)

{
    Threshold[0] = _Threshold0;
    Threshold[1] = _Threshold1;

    N_neuronsL[0] = _NL0;
    N_neuronsL[1] = _NL1;
    N_neurons = N_neuronsL[0] + N_neuronsL[1];
    N_streams = N_InputStreams + _NL0;
    tmax = tau_s * tau_m / (tau_m - tau_s) * (log(tau_m) - log(tau_s));
    MaxDeltaT = 7. * tau_m;

    fire_granularity = tau_s / 5;
    fire_precision = Threshold[0] / 100;
    myRNG = new TRandom3(23);
    largenumber = 999999999.;
    epsilon = 1. / largenumber;

    Weight = new float *[N_neurons];         // Weight of synapse-neuron strength
    check_LTD = new bool *[N_neurons];       // checks to generate LTD after neuron discharge
    Void_weight = new bool *[N_neurons];     // These may be used to model disconnections
    Delay = new float *[N_neurons];          // Delay in incoming signals

    for (int in = 0; in < N_neurons; in++)
    {
        Weight[in] = new float[N_streams];
        check_LTD[in] = new bool[N_streams];
        Void_weight[in] = new bool[N_streams];
        Delay[in] = new float[N_streams];
    }

    History_time = new vector<float>[N_neurons]; // Time of signal events per each 1neuron
    History_type = new vector<int>[N_neurons];   // Type of signal
    History_ID = new vector<int>[N_neurons];     // ID of generating signal stream or neuron
    Fire_time = new vector<float>[N_neurons];    // Times of firing of each neuron
    Neuron_layer = new int[N_neurons];
    sumweight = new float[N_neurons]; // summed weights of streams for each neurons for the purpose of normalization

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
    cout << "-------------------------------------" << endl;

    cout << "Init_neurons" << endl;
    Init_neurons();
    cout << "Init_connection_map" << endl;
    Init_connection_map();
    cout << "Init_Weights" << endl;
    Init_weights();
    cout << "Done" << endl;
}

SNN::~SNN()
{
}

void SNN::Init_delays(){
    return;
}
void SNN::Reset_weights(){
    return;
}

float SNN::bisectionMethod(float a, float b, int in, float epsilon, std::function<float(int, float, bool)> func)
{
    float fa = func(in, a, false);
    float fb = func(in, b, false);
    float c = 0;
    // cout << "Initial values: fa = " << fa << ", fb = " << fb << endl;

    if (fa * fb > 0)
    {
        cerr << "Error: The function values at the endpoints have the same sign. Interval: [" << a << ", " << b << "]\n";
        cerr << "fa: " << fa << " fb: " << fb << "in: " << in << endl;
        return largenumber; // Indicate failure
    }

    int maxIterations = 10; // Choose an appropriate maximum number of iterations

    for (int i = 0; i < maxIterations; ++i)
    {
        c = (a + b) / 2;
        if (fa * fb > 0)
        {
            cout << "---------- Bisection problem -------------" << endl;
            cout << "fa: " << fa << " fb: " << fb << endl;
            return largenumber;
        }

        float fc = func(in, c, false);

        // cout << "Iteration " << i << ": Interval [" << a << ", " << b << "], Root estimate: " << c << ", Function value: " << fc << endl;

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

    return c; // Indicate failure
}

// Initialize neuron potentials
// ----------------------------

void SNN::Init_neurons()
{
    cout << N_neurons << endl;
    for (int in = 0; in < N_neurons; in++)
    {
        // Set first event in history of this neuron
        History_time[in].clear();
        History_type[in].clear();
        History_ID[in].clear();

        History_time[in].push_back(0);
        History_type[in].push_back(0);
        History_ID[in].push_back(0);
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
        }
    }

    return;
}

/*

void SNN::Init_delays()
{
    // Define delays for IE signals
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
            Delay[in][is] = 0.;
            if (learnDelays || updateDelays)
                Delay[in][is] = MaxDelay / 2.;
            //            if (is<N_bin_r) { // no IE delay for neuron-originated spikes into L1
            //                Delay[in][is] = myRNG->Uniform(MaxDelay);
            //            }
        }
    }
    return;
}
*/

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
        }
    }
    PrintWeights();
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
    cout << endl;
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
            cout << Void_weight[in][is];
        cout << endl;
    }

    return;
}

// Model Excitatory Post-Synaptic Potential
// We take this as parametrized in T. Masquelier et al., "Competitive STDP-Based Spike Pattern Learning", DOI: 10.1162/neco.2008.06-08-804
// ---------------------------------------------------------------------------------------------------------------------------------------
float SNN::EPS_potential(float delta_t)
{
    float psp = 0.;
    if (delta_t >= 0. && delta_t < MaxDeltaT)
        psp = K * (exp(-delta_t / tau_m) - exp(-delta_t / tau_s));
    return psp;
}

// Model membrane potential after spike
// Also modeled as in paper cited above, like IPSP and other signals below
// -----------------------------------------------------------------------
float SNN::Spike_potential(float delta_t, int ilayer)
{
    float sp = 0.;
    if (delta_t >= 0. && delta_t < MaxDeltaT)
        sp = Threshold[ilayer] * (K1 * exp(-delta_t / tau_m) - K2 * (exp(-delta_t / tau_m) - exp(-delta_t / tau_s)));
    return sp;
}

// Model Inhibitory Post-Synaptic Potential (IPSP)
// -----------------------------------------------
float SNN::Inhibitory_potential(float delta_t, int ilayer)
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
float SNN::Neuron_firetime(int in, float t)
{
    float t0 = History_time[in][0];
    float delta_t = t - t0;
    if (delta_t < tau_r)
        return largenumber;

    int ilayer = Neuron_layer[in];
    // now we will scan the interval in between the last EPSP and this time looking for an activation according to the defined granularity
    float last_EPSP = -1;

    for (int ih = History_type[in].size() - 1; ih > 1; ih--)
    {
        // longer approach: add "&& History_ID[in][ih] < N_InputStreams" to rescan from the last InputStream spike
        if (History_type[in][ih] == 1 && !Void_weight[in][History_ID[in][ih]])
        {
            last_EPSP = History_time[in][ih];
            break;
        }
    }
    // in that case it's impossible that the neuron is firing
    if (last_EPSP < 0)
        return largenumber;

    // now I want to scan the potential from the last EPSP to time t at fire_granularity steps
    float t_neg = last_EPSP;
    float time = last_EPSP + fire_granularity;
    float P_neg = Neuron_Potential(in, t_neg, true);
    if (P_neg > Threshold[ilayer] && t_neg - t0 > tau_r)
        return t_neg;

    float P_t = 0;
    bool fire = false;
    while (time < t)
    {
        P_t = Neuron_Potential(in, time, false);
        // If I a value below the threshold I save it
        if (P_t < Threshold[ilayer] || t_neg - t0 < tau_r)
        {
            P_neg = P_t;
            t_neg = time;
        }
        // if I find a value higher than the treshold a neuron is gonna fire!
        else
        {
            fire = true;
            break;
        }
        time += fire_granularity;
    }
    // Let's check at the time t
    if (!fire)
    {
        time = t;
        P_t = Neuron_Potential(in, time, false);
        // if still the potential is below the threshold we know that the neuron is not activating
        if (P_t < Threshold[ilayer])
            return largenumber;
    }
    if (P_neg > Threshold[ilayer])
    {
        cout << "Weird" << endl;
        cout << fire << ", " << time - t_neg << endl;
    }

    // if we are here the potential has reached the threshold at some point!
    // we need to determine when the neuron has fired given a certain confidence
    return bisectionMethod(t_neg, time, in, fire_precision,
                           [this](int in, float time, bool delete_history)
                           {
                               return Neuron_Potential(in, time, delete_history) - Threshold[Neuron_layer[in]];
                           });
}

// Compute collective effect of excitatory, post-spike, and inhibitory potentials on a neuron
// ------------------------------------------------------------------------------------------
float SNN::Neuron_Potential(int in, float t, bool delete_history)
{
    int ilayer = Neuron_layer[in];
    float P0 = 0.;
    float t0 = History_time[in][0];
    float delta_t = t - t0;
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
            if (History_type[in][ih] == 1)
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
                    len = len - 1;
                }
            }
            else if (History_type[in][ih] == 2)
            { // IPSP
                if (delta_t < MaxDeltaT)
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
                    len = len - 1;
                }
            }
            else if (History_type[in][ih] == 3)
            { // IE
                if (delta_t < MaxDeltaT && (History_time[in][ih] - t0) > tau_r)
                {
                    if (!Void_weight[in][History_ID[in][ih]]) // for type 1 or 3 signals, ID is the stream
                        P += IE_potential(delta_t, in, History_ID[in][ih]);
                }
                else if (delete_history)
                {
                    // get rid of irrelevant events
                    History_time[in].erase(History_time[in].begin() + ih, History_time[in].begin() + ih + 1);
                    History_type[in].erase(History_type[in].begin() + ih, History_type[in].begin() + ih + 1);
                    History_ID[in].erase(History_ID[in].begin() + ih, History_ID[in].begin() + ih + 1);
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
float SNN::IE_potential(float delta_t, int in, int is)
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

void SNN::LTP(int in, int this_spike, float fire_time, bool nearest_spike_approx, SNN &old)
{
    for (int is = 0; is < N_streams; is++)
        {            
        if (Void_weight[in][is])
            break;
        check_LTD[in][is] = true;
        // Use nearest-spike approximation: search for closest pre-spike
        bool no_prespikes = true;
        int isp = this_spike - 1;
        float delta_weight = 0;
        do
        {
            if (History_ID[in][isp] == is && History_type[in][isp]==1)
            {
                float delta_t = History_time[in][isp] - fire_time;
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

        if (!no_prespikes)
            Renorm(in, delta_weight, old);
        }
    return;
}

void SNN::LTD(int in, int is, float spike_time, bool nearest_spike_approx, SNN &old)
{
    if(!check_LTD[in][is]) return;
    if (Fire_time[in].size() == 0)
        return;
    if (Void_weight[in][is])
        return;
    float delta_t = spike_time - Fire_time[in].back();
    // if nearest_spike_approx we prevent to compute future LTD until the next activation of the neuron
    if (nearest_spike_approx)
        check_LTD[in][is] = false;

    if (delta_t >= 0 && delta_t < 7. * tau_minus)
    {
        Weight[in][is] -= a_minus * exp(-delta_t / tau_minus);
        if (Weight[in][is] < 0.)
            Weight[in][is] = 0.;
        Renorm(in, Weight[in][is] - old.Weight[in][is], old);
    }
    return;
}

void SNN::Renorm(int in, float delta_weight, SNN &old)
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
            if (!Void_weight[in][is])
            cout << Weight[in][is] << ", ";
        }
        cout << endl <<endl;
    }
}