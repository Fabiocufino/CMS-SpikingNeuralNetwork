#include "SNN.h"
#include "Snnt_constants.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TF1.h"

using namespace std;


SNN::SNN(int NL0, int NL1): //Initializations
                            int N_ev(100000), int N_ep(10), char *rootInput("100k_100br.root"), bool batch(false),  float _tau_m(1e-09), double _tau_s(0.25e-09),
                  double _tau_plus(1.68e-09), double _tau_minus(3.37e-09), double _a_plus(0.03125), floa _a_minus(0.02656), double CFI0(1), double CFI1(1), double CF01(1), double a(0.25),
                  double Thresh0 (15), double Thresh1 (10), double _MaxFactor (0.2), double l1if (1.), double k (1.), double k1(2.), double k2 (4.),
                  double IEPC (2.5), double ipspdf (1.0), double _MaxDelay (0.1e-9),
                  int N_cl(6),
                  int TrainingCode(4), bool ReadPars(false), long int _NROOT(100000)
{
    Threshold[0] = 0.1;
    Threshold[1] = 0.1;


    N_neuronsL[0] = NL0;
    N_neuronsL[1] = NL1;
    //fire_granularity = tau_s/5;
    //fire_precision = 1e-11;
    
    //Print all the variables
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
    cout << "MaxFactor = " << MaxFactor << endl;
    cout << "Threshold[0] = " << Threshold[0] << endl;
    cout << "Threshold[1] = " << Threshold[1] << endl;
    cout << "tmax = " << tmax << endl;
    cout << "MaxDeltaT = " << MaxDeltaT << endl;
    cout << "N_neuronsL[0] = " << N_neuronsL[0] << endl;
    cout << "N_neuronsL[1] = " << N_neuronsL[1] << endl;
    cout << "N_InputStreams = " << N_InputStreams << endl;
    cout << "NROOT = " << NROOT << endl;
    cout << "MaxClasses = " << MaxClasses << endl;
    cout << "BGR = " << BGR << endl;
    cout << "SIG = " << SIG << endl;
    cout << "MaxNeurons = " << MaxNeurons << endl;
    cout << "largenumber = " << largenumber << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "MaxEvents = " << MaxEvents << endl;
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
//intiialisation

void SNN::Init_neurons()
{
    for (int in = 0; in < N_neurons; in++)
    {
        // Set first event in history of this neuron
        History_time[in].push_back(0.);
        History_type[in].push_back(0);
        History_ID[in].push_back(0);
        if (in < N_neuronsL[0])
        {
            Neuron_layer[in] = 0;
        }
        else
        {
            Neuron_layer[in] = 1;
        }
    }
    return;
}

// Initialize synapse weights
// --------------------------


void SNN::Init_weights()
{
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
            check_LTD[in][is] = true; // flags used to see if we need to create a LTD signal after a neuron discharge
            Weight[in][is] = myRNG->Uniform();
           /// Weight_initial[infor (int is=0; is< Nsynapses (or however this is called); is++) {
 
 
            sumweight[in]+=Weight[in][is];
            
             
   
        }
    }
    
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
          if(sumweight[in]>0)
             {
           Weight[in][is]=Weight[in][is]/sumweight[in];
           checksumweight[in]+=Weight[in][is];
           Weight_initial[in][is] = Weight[in][is];
           OldWeight[in][is]=Weight[in][is];//this will be used for the renorm
             }
        }
    
    cout<<"ini weightsum "<<checksumweight[in]<<endl;
    }
    
    return;

}
//Renormalization function
    
     
void SNN::Renorm(int in) {
    double weight_sum = 0.0;
    double checks=0;
    // Calculate the sum of weights for the 'in' neuron
    for (int is = 0; is < N_streams; is++) {
        weight_sum += Weight[in][is];
        //cout<<"wieght sum "<<weight_sum<<endl;
    }
    
    
    if (weight_sum > 0.0) {
        // Normalize the weights for the 'in' neuron
        for (int is = 0; is < N_streams; is++) {
            Weight[in][is] =Weight[in][is]/weight_sum;
           // cout<<"wiegh ltp "<<Weight[in][is]<<endl;
            checks+=Weight[in][is];
        }
        
    // cout<<"sum "<<checks<<endl;   
    }
return;    
}

 
            
 
// Reset synapse weights
// ---------------------
void SNN::Reset_weights()
{
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
            Weight[in][is] = Weight_initial[in][is];
        }
    }
    return;
}

// Initialize synaptic delays
// --------------------------
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
            if (myRNG->Uniform() > ConnectedFraction_Input_L0)
                Void_weight[in][is] = true;
        }
        // input connections L0 -> L0
        for (int is = N_InputStreams; is < N_streams; is++)
        {
            Void_weight[in][is] = false;
        }
    }

    // Setting L1 input connections
    for (int in = N_neuronsL[0]; in < N_neurons; in++)
    {
        // input connections Tracking layers -> L1
        for (int is = 0; is < N_InputStreams; is++)
        {
            Void_weight[in][is] = false;
            if (myRNG->Uniform() > ConnectedFraction_Input_L1)
                Void_weight[in][is] = true;
        }
        // input connections L0 -> L1
        for (int is = N_InputStreams; is < N_streams; is++)
        {
            Void_weight[in][is] = false;
            if (myRNG->Uniform() > ConnectedFraction_L0_L1)
                Void_weight[in][is] = true;
        }
    }
    return;
}

// clera hits vector
void Reset_hits()
{
    hit_pos.clear();
    return;
}

// Start binning r-z plane ---------------
int GetBinR(float r_hit)
{
    if (r_hit < 0)
        r_hit = 0;
    if (r_hit > max_R)
        r_hit = max_R - epsilon;

    return (int)(r_hit / r_bin_length);
}

int GetBinZ(float z)
{
    double tmp = (z + z_range / 2.);
    if (tmp < 0)
        tmp = 0;
    else if (tmp > z_range)
        tmp = z_range - epsilon;

    return (int)(tmp / z_bin_length);
}

int GetStreamID(int r, int z)
{
    return r + z * N_bin_r; // Visually streams are sorted by z and then by r
    // return r*N_bin_z+z;  // Visually streams are sorted by r and then by z
}
// End binning r-z plane -----------------

// Transforms hits into spikes streams by scanning the event
void Encode(double t_in)
{
    // sort by phi angle
    sort(hit_pos.begin(), hit_pos.end(), [](const Hit &h1, const Hit &h2)
         { return h1.phi < h2.phi; });

    for (auto &&row : hit_pos)
    {
        double time = t_in + row.phi / omega;
        int itl = GetStreamID(GetBinR(row.r), GetBinZ(row.z));

        PreSpike_Time.push_back(time);
        PreSpike_Stream.push_back(itl);
        PreSpike_Signal.push_back(row.id - 1); // 0,1,2 -> -1,0,1 respectively NoHit, Backgroung, Signal
    }

    // rescan from [0, delta]
    for (auto &&row : hit_pos)
    {
        if (row.phi > delta)
            break;
        double time = t_in + (row.phi + M_PI * 2.) / omega;
        // uncomment when implemented:
        int itl = GetStreamID(GetBinR(row.r), GetBinZ(row.z));

        PreSpike_Time.push_back(time);
        PreSpike_Stream.push_back(itl);
        PreSpike_Signal.push_back(row.id - 1); // 0,1,2 -> -1,0,1 respectively NoHit, Backgroung, Signal
    }
}

// Model Excitatory Post-Synaptic Potential
// We take this as parametrized in T. Masquelier et al., "Competitive STDP-Based Spike Pattern Learning", DOI: 10.1162/neco.2008.06-08-804
// ---------------------------------------------------------------------------------------------------------------------------------------
float SNN::EPS_potential(double delta_t)
{
    double psp = 0.;
    if (delta_t >= 0. && delta_t < MaxDeltaT)
        psp = K * (exp(-delta_t / tau_m) - exp(-delta_t / tau_s));
    return psp;
}

// Model membrane potential after spike
// Also modeled as in paper cited above, like IPSP and other signals below
// -----------------------------------------------------------------------
float SNN::Spike_potential(double delta_t, int ilayer)
{
    double sp = 0.;
    if (delta_t >= 0. && delta_t < MaxDeltaT)
        sp = Threshold[ilayer] * (K1 * exp(-delta_t / tau_m) - K2 * (exp(-delta_t / tau_m) - exp(-delta_t / tau_s)));
    return sp;
}

// Model Inhibitory-Excitatory signal (IE) as combination of two EPSP-like shapes, a negative one followed by a positive one
// This is a crude model, loosely inspired by shapes in Fig.2 of F. Sandin and M. Nilsson, "Synaptic Delays for Insect-Inspired
// Feature Detection in Dynamic Neuromorphic Processors", doi.org/10.3389/fnins.2020.00150
// ----------------------------------------------------------------------------------------------------------------------------
float SNN::IE_potential(double delta_t, int in, int is)
{
    double sp = 0.;
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

// Model Inhibitory Post-Synaptic Potential (IPSP)
// -----------------------------------------------
float SNN::Inhibitory_potential(double delta_t, int ilayer)
{
    // In order to dilate the inhibition to larger times (in the attempt at obtaining higher selectivity),
    // we kludge it by multiplying delta_t and maxdeltat by a factor
    delta_t = delta_t * IPSP_dt_dilation;
    double ip = 0.;
    double thisalpha = alpha;
    if (ilayer > 0)
        thisalpha = L1inhibitfactor * alpha; // Different inhibition in L1
    if (delta_t >= 0. && delta_t < MaxDeltaT)
        ip = -thisalpha * Threshold[ilayer] * EPS_potential(delta_t);
    return ip;
}

// Model Spike-Timing-Dependent Plasticity (LTP - Long-Term Potentiation) - modify weights
// based on how close a previous synapse fired before the spike)
// ---------------------------------------------------------------------------------------
void SNN::LTP(int in, int is, int this_spike, double fire_time)
{
    if (Void_weight[in][is])
        return;
    // Use nearest-spike approximation: search for closest pre-spike
    bool no_prespikes = true;
    int isp = this_spike - 1;
    do
    {
        if (PreSpike_Stream[isp] == is)
        {
            double delta_t = PreSpike_Time[isp] - fire_time;
            Weight[in][is] += a_plus * exp(delta_t / tau_plus);
            if (Weight[in][is] > 1.)
                Weight[in][is] = 1.;
            no_prespikes = false;
        }
        isp--;
    } while (isp >= 0 && PreSpike_Time[isp] > fire_time - 7. * tau_plus && no_prespikes);

    // Also modify delays (experimental)
    if (learnDelays)
    {
        double dtbig = MaxDelay / NevPerEpoch;
        if (Delay[in][is] >= dtbig)
            Delay[in][is] -= dtbig;
        double dt = dtbig / N_neuronsL[Neuron_layer[in]];
        int inmin = 0;
        int inmax = N_neuronsL[0];
        if (Neuron_layer[in] == 1)
        {
            inmin = N_neuronsL[0];
            inmax = N_neurons;
        }
        for (int in2 = inmin; in2 < inmax; in2++)
        {
            if (in2 != in && Delay[in2][is] <= MaxDelay - dt)
                Delay[in2][is] += dt;
        }
    }
    return;
}

// Model Spike-Timing-Dependent Plasticity (LTD - Long-Term Depression)
// --------------------------------------------------------------------
void SNN::LTD(int in, int is, double spike_time)
{
    // Use nearest-spike approximation: search for closest neuron spike to this input spike
    if (Fire_time[in].size() == 0)
        return;
    if (Void_weight[in][is])
        return;
    double delta_t = spike_time - Fire_time[in].back();
    check_LTD[in][is] = false;
    if (delta_t >= 0 && delta_t < 7. * tau_minus)
    {
        Weight[in][is] -= a_minus * exp(-delta_t / tau_minus);
        if (Weight[in][is] < 0.)
            Weight[in][is] = 0.;
        //cout<<"LTD "<<Weight[in][is]<<endl;    
    }
    // Also modify delay (experimental)
    if (learnDelays)
    {
        double dtbig = MaxDelay / NevPerEpoch;
        if (Delay[in][is] < MaxDelay - dtbig)
            Delay[in][is] += dtbig;
        double dt = dtbig / (N_neurons - 1);
        for (int in2 = 0; in2 < N_neurons; in2++)
        {
            if (in2 != in && Neuron_layer[in] == Neuron_layer[in2] && Delay[in2][is] > dt)
                Delay[in2][is] -= dt;
        }
    }
    return;
}

// Compute collective effect of excitatory, post-spike, and inhibitory potentials on a neuron
// ------------------------------------------------------------------------------------------
float SNN::Neuron_firetime(int in, double t)
{
    int ilayer = Neuron_layer[in];
    double P0 = 0.;
    double t0 = History_time[in][0];
    double delta_t = t - t0;
    if (t0 > 0. && delta_t >= 0. && delta_t < MaxDeltaT)
    {
        int ilayer = Neuron_layer[in];
        P0 = Spike_potential(delta_t, ilayer); // the first event in the history sequence is a spike
    }
    double P = P0;

    // Now we extrapolate the effect of all past spikes and inhibitions to time t, to compute the potential when EPSP arrives
    int len = History_time[in].size();
    if (len > 1)
    {
        for (int ih = 1; ih < len - 1; ih++)
        {
            delta_t = t - History_time[in][ih];
            if (History_type[in][ih] == 1)
            { // EPSP
                if (delta_t > MaxDeltaT)
                { // Get rid of irrelevant events
                    History_time[in].erase(History_time[in].begin() + ih, History_time[in].begin() + ih + 1);
                    History_type[in].erase(History_type[in].begin() + ih, History_type[in].begin() + ih + 1);
                    History_ID[in].erase(History_ID[in].begin() + ih, History_ID[in].begin() + ih + 1);
                    len = len - 1;
                }
                else
                {
                    if (!Void_weight[in][History_ID[in][ih]]) // for type 1 or 3 signals, ID is the stream
                        P += Weight[in][History_ID[in][ih]] * EPS_potential(delta_t);
                }
            }
            else if (History_type[in][ih] == 2)
            { // IPSP
                if (delta_t > MaxDeltaT)
                { // get rid of irrelevant events
                    History_time[in].erase(History_time[in].begin() + ih, History_time[in].begin() + ih + 1);
                    History_type[in].erase(History_type[in].begin() + ih, History_type[in].begin() + ih + 1);
                    History_ID[in].erase(History_ID[in].begin() + ih, History_ID[in].begin() + ih + 1);
                    len = len - 1;
                }
                else
                {
                    int ilayer = Neuron_layer[in];
                    P += Inhibitory_potential(delta_t, ilayer);
                }
            }
            else if (History_type[in][ih] == 3)
            { // IE
                if (delta_t > MaxDeltaT)
                { // get rid of irrelevant events
                    History_time[in].erase(History_time[in].begin() + ih, History_time[in].begin() + ih + 1);
                    History_type[in].erase(History_type[in].begin() + ih, History_type[in].begin() + ih + 1);
                    History_ID[in].erase(History_ID[in].begin() + ih, History_ID[in].begin() + ih + 1);
                    len = len - 1;
                }
                else
                {
                    if (!Void_weight[in][History_ID[in][ih]]) // for type 1 or 3 signals, ID is the stream
                        P += IE_potential(delta_t, in, History_ID[in][ih]);
                }
            }
        }
    }
    if (P > Threshold[ilayer])
    { // Neuron will fire as spike contribution will bring it above threshold
        // compute fire time by looping more finely from t to t+tmax (tmax is peak time of EPSP)
        double this_t = t;
        do
        {
            P = P0;
            for (int ih = 1; ih < len; ih++)
            {
                double delta_t = this_t - History_time[in][ih];
                if (History_type[in][ih] == 1)
                { // EPSP
                    if (!Void_weight[in][History_ID[in][ih]])
                        P += Weight[in][History_ID[in][ih]] * EPS_potential(delta_t);
                }
                else if (History_type[in][ih] == 2)
                { // IPSP
                    int ilayer = Neuron_layer[in];
                    P += Inhibitory_potential(delta_t, ilayer);
                }
                else if (History_type[in][ih] == 3)
                { // IE
                    if (!Void_weight[in][History_ID[in][ih]])
                        P += IE_potential(delta_t, in, History_ID[in][ih]);
                }
            }
            this_t += 1 / (10000. * omega);
        } while (P < Threshold[ilayer] && this_t <= t + tmax);
        if (P >= Threshold[ilayer])
            return this_t;
    }
    return largenumber;
}

// Learning rate scheduler - this returns an oscillating, dampened function as a function of the epoch
// ---------------------------------------------------------------------------------------------------
float SNN::LR_Scheduler(double LR0, int epoch, int Nepochs)
{
    double par[3] = {-0.01, 0.2, 0.2};
    double x = 100. * epoch / Nepochs;
    return LR0 * exp(par[0] * x) * (par[1] + (1. - par[1]) * pow(cos(par[2] * x), 2));
}

// Calculate selectivity of set of neurons
// ---------------------------------------
float SNN::Compute_Selectivity(int level, int mode)
{
    double S = 0.;
    int inmin, inmax;
    if (level == 0)
    {
        inmin = 0;
        inmax = N_neuronsL[0];
    }
    else if (level == 1)
    {
        inmin = N_neuronsL[0];
        inmax = N_neurons;
    }
    if (mode == 0)
    { // Use additive rule
        for (int in = inmin; in < inmax; in++)
        {
            // select max efficiency class
            double maxeff = 0.;
            double sumeff = 0.;
            for (int ic = 0; ic < N_classes; ic++)
            {
                double e = Eff[ic + N_classes * in];
                if (e > maxeff)
                    maxeff = e;
                sumeff = sumeff + e;
            }
            if (sumeff > 0.)
                S += maxeff * N_classes / sumeff;
        }
        if (inmax - inmin > 0)
            S = S / (inmax - inmin);
    }
    else if (mode == 1)
    { // aim for collective effect (each neuron participates)
        S = 1.;
        for (int in = inmin; in < inmax; in++)
        {
            // select max efficiency class
            double maxeff = 0.;
            double sumeff = 0.;
            for (int ic = 0; ic < N_classes; ic++)
            {
                double e = Eff[ic + N_classes * in];
                if (e > maxeff)
                    maxeff = e;
                sumeff = sumeff + e;
            }
            if (sumeff > 0.)
                S *= maxeff / sumeff;
        }
    }
    else if (mode == 2)
    { // compute mutual information
        // I(N_neurons,N_classes) = Sum_i^N_n Sum_j^N_c Eff(i,j) log_2 [Eff(i,j)/Eff_i Eff_j)]
        // where  is the average efficiency of neuron i over classes, and Eff_j is the
        // average efficiency on class j over neurons
        S = 0.;
        double Effn[MaxNeurons];
        double Effc[MaxClasses];
        double sumeff = 0.;
        for (int in = inmin; in < inmax; in++)
        {
            sumeff = 0.;
            for (int ic = 0; ic < N_classes; ic++)
            {
                sumeff += Eff[ic + N_classes * in];
            }
            if (N_classes > 0)
                Effn[in] = sumeff / N_classes;
        }
        for (int ic = 0; ic < N_classes; ic++)
        {
            sumeff = 0.;
            for (int in = inmin; in < inmax; in++)
            {
                sumeff += Eff[ic + N_classes * in];
            }
            if (inmax > inmin)
                Effc[ic] = sumeff / (inmax - inmin);
        }
        for (int in = inmin; in < inmax; in++)
        {
            for (int ic = 0; ic < N_classes; ic++)
            {
                if (Effc[ic] * Effn[in] > 0.)
                    S += Eff[ic + N_classes * in] *
                         (log2(Eff[ic + N_classes * in] + epsilon) - log2(Effc[ic] * Effn[in]));
            }
        }
    }
    if (inmax > inmin)
        S /= N_classes * (inmax - inmin);
    return S;
}

// Compute Q-value
// ---------------
float SNN::Compute_Q(double eff, double acc, double sel)
{
    // efficiency saturates at eff_target
    if (eff > eff_target)
        eff = eff_target + pow(eff - eff_target, 1.5);
    // acceptance saturates at acc_target
    if (acc < acc_target)
        acc = acc_target - pow(acc_target - acc, 1.5);
    double Q0 = eff / sqrt(acc);
    double w = 5. * (exp(Q0 / 4.) - 1.) / (exp(1.) - 1.);
    return Q0 + w * sel;
}
