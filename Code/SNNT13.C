#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "Riostream.h"
#include "Snnt_constants.h"

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <string>

#include "class_implementation/SNN.h"

using namespace std;

// Constants and data used throughout the code
// -------------------------------------------
static int N_part; // Number of generated particles in an event
static float First_angle;
static float Eff[MaxNeurons * MaxClasses];       // Efficiency of each neuron to signals of different classes
static float Efftot[MaxClasses];                 // Global efficiency to a different class
float SumofSquaresofWeight[MaxNeurons] = {0};  // sum of squares synaptic weights for each neuron for RMS calc
float MeanofSquaresofWeight[MaxNeurons] = {0}; // mean of squares of synaptic weights for each neuron for RMS calc
float MaxWeight[MaxNeurons];
float MinWeight[MaxNeurons];
float RMSWeight[MaxNeurons];
float Efficiency_track[MaxNeurons][MaxClasses][10]; // to keep track of the efficiency evolution during training
float FK_track[MaxNeurons][MaxClasses][10];         // to keep track of the FakeRate evolution during training
static vector<float> PreSpike_Time;
static vector<int> PreSpike_Stream;
static vector<int> PreSpike_Signal; // 0 for background hit, 1 for signal hit, 2 for L1 neuron spike
static vector<int> neurons_index;   // contains neurons identifiers in random positions
static float Q_best;
static float SelL1_best;
static float Eff_best;
static float Acc_best;
static int indfile;
static char progress[53] = "[--10%--20%--30%--40%--50%--60%--70%--80%--90%-100%]"; // Progress bar
static long int ievent;

// New random number generator
static TRandom3 *myRNG = new TRandom3(23);

// Function that reads parameters from file, if you want to start from a previous optimization
// ----------------------------------------
/*


int Read_Parameters()
{
    string Path = "./MODE/SNNT/";
    ifstream tmpfile;
    indfile = -1;
    // Determine last available file number to read from, by attempting to open all files with same name and previous numbering
    do
    {
        if (indfile > -1)
            tmpfile.close();
        indfile++;
        stringstream tmpstring;
        tmpstring << "Params13_NL0=" << N_neuronsL[0] << "_NL1=" << N_neuronsL[1] << "_NCl=" << N_classes << "_CF01=" << CF01 << "_CFI0=" << CFI0 << "_CFI1=" << CFI1 << "_" << indfile;
        string tmpfilename = Path + tmpstring.str() + ".txt";
        tmpfile.open(tmpfilename);
    } while (tmpfile.is_open());

    if (indfile == -1)
    {
        cout << "     Warning, no file to read parameters from. " << endl;
        return -1;
    }
    ifstream parfile;
    stringstream sstr;
    char num[80];
    sprintf(num, "NL0=%d_NL1=%d_NCl=%d_CF01=%.2f_CFI0=%.2f_CFI1=%.2f_alfa=%.2f_%d", N_neuronsL[0], N_neuronsL[1], N_classes, CF01, CFI0, CFI1, alpha, indfile - 1); // we'll pick the last one in the list
    sstr << "Params13_";
    string nameparfile = Path + sstr.str() + num + ".txt";
    parfile.open(nameparfile);
    float e;
    int ie;
    parfile >> ie;
    if (ie != N_neurons)
        cout << "Warning, file " << nameparfile << " ie= " << ie << " N_neurons = " << N_neurons << " - input file not matching N_neurons" << endl;
    parfile >> ie;
    if (ie != N_streams)
        cout << "Warning, file " << nameparfile << " ie= " << ie << " N_streams = " << N_streams << " - input file not matching N_streams" << endl;
    parfile >> e;
    Threshold[0] = e;
    parfile >> e;
    Threshold[1] = e;
    parfile >> e;
    alpha = e;
    parfile >> e;
    L1inhibitfactor = e;
    parfile >> e;
    K = e;
    parfile >> e;
    K1 = e;
    parfile >> e;
    K2 = e;
    parfile >> e;
    IE_Pot_const = e;
    parfile >> e;
    IPSP_dt_dilation = e;
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
            parfile >> e;
            Delay[in][is] = e;
        }
    }
    bool b;
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
            parfile >> b;
            Void_weight[in][is] = b;
            if (in < N_neuronsL[0] && is >= N_InputStreams)
                Void_weight[in][is] = true; // just making sure
        }
    }
    parfile.close();
    return 0;
}

// Function that saves the layout data to file
// -------------------------------------------
void Write_Parameters()
{

    string Path = "./MODE/SNNT/";
    ifstream tmpfile;
    indfile = -1;
    // Determine first available file number to write, by attempting to open all files with same name and previous numbering
    do
    {
        if (indfile > -1)
            tmpfile.close();
        indfile++;
        stringstream tmpstring;
        tmpstring << "Params13_NL0=" << N_neuronsL[0] << "_NL1=" << N_neuronsL[1] << "_NCl=" << N_classes << "_CF01=" << CF01 << "_CFI0=" << CFI0 << "_CFI1=" << CFI1 << "_" << indfile;
        string tmpfilename = Path + tmpstring.str() + ".txt";
        tmpfile.open(tmpfilename);
    } while (tmpfile.is_open());

    ofstream parfile;
    stringstream sstr;
    char num[80];
    sprintf(num, "NL0=%d_NL1=%d_NCl=%d_CF01=%.2f_CFI0=%.2f_CFI1=%.2f_alfa=%.2f_%d", N_neuronsL[0], N_neuronsL[1], N_classes, CF01, CFI0, CFI1, alpha, indfile);
    sstr << "Params13_";
    string nameparfile = Path + sstr.str() + num + ".txt";
    parfile.open(nameparfile);
    parfile << N_neurons << " " << N_streams << endl;
    parfile << T0_best << " " << T1_best << " " << A_best << " "
            << L1if_best << " " << K_best << " " << K1_best << " " << K2_best << " " << IEPC_best << " " << IPSP_dt_dilation << endl;
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
            parfile << bestDelay[in][is] << endl;
        }
    }
    for (int in = 0; in < N_neurons; in++)
    {
        for (int is = 0; is < N_streams; is++)
        {
            parfile << bestVoid_weight[in][is] << endl;
        }
    }

    // Also write optimization output
    parfile << Q_best << " "
            << " " << Eff_best << " " << Acc_best << " " << SelL1_best << endl;

    // Finally, write complete set of hyperparameters and settings
    parfile << "                       L0 neurons: " << N_neuronsL[0] << endl;
    parfile << "                       L1 neurons: " << N_neuronsL[1] << endl;
    parfile << "            Connected L0-L1 frac.: " << CF01 << endl;
    parfile << "            Connected IN-L0 frac.: " << CFI0 << endl;
    parfile << "            Connected IN-L1 frac.: " << CFI1 << endl;
    parfile << "                    Track classes: " << N_classes << endl;
    parfile << "                     Total events: " << N_events << endl;
    parfile << "               Optimization loops: " << N_epochs << endl;
    parfile << "             Optimize SNN params.: ";
    if (update9)
    {
        parfile << "True" << endl;
    }
    else
    {
        parfile << "False" << endl;
    }
    parfile << "                  Optimize delays: ";
    if (updateDelays)
    {
        parfile << "True" << endl;
    }
    else
    {
        parfile << "False" << endl;
    }
    parfile << "             Optimize connections: ";
    if (updateConnections)
    {
        parfile << "True" << endl;
    }
    else
    {
        parfile << "False" << endl;
    }
    parfile << "                  Max mod. factor: " << MaxFactor << endl;
    parfile.close();
    return;
}
*/

// Reset synapse weights
// ---------------------

// void Reset_weights(SNN &snn_in, SNN &snn_old)
// {
//     for (int in = 0; in < snn_in.N_neurons; in++)
//     {
//         for (int is = 0; is < snn_in.N_streams; is++)
//         {
//             snn_in.Weight[in][is] = snn_old.Weight[in][is];
//         }
//     }
//     return;
// }



// clear hits vector
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

    // 10 bins in even positions are associated to a tracking layer
    // 11 bins in odd positions are associated to empty space among the former
    for (int i = 0; i < N_TrackingLayers; i++)
    {
        if (r_hit > Left_Layers[i] && r_hit < Right_Layers[i])
            return 2 * i + 1;
        else if (r_hit < Left_Layers[i])
            return 2 * i;
    }
    return N_bin_r - 1;
}

int GetBinZ(float z)
{
    float tmp = (z + z_range / 2.);
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
void Encode(float t_in)
{
    // sort by phi angle
    sort(hit_pos.begin(), hit_pos.end(), [](const Hit &h1, const Hit &h2)
         { return h1.phi < h2.phi; });

    for (auto &&row : hit_pos)
    {
        float time = t_in + row.phi / omega;
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
        float time = t_in + (row.phi + M_PI * 2.) / omega;

        int itl = GetStreamID(GetBinR(row.r), GetBinZ(row.z));

        PreSpike_Time.push_back(time);
        PreSpike_Stream.push_back(itl);
        PreSpike_Signal.push_back(row.id - 1); // 0,1,2 -> -1,0,1 respectively NoHit, Backgroung, Signal
    }
}

// Learning rate scheduler - this returns an oscillating, dampened function as a function of the epoch
// ---------------------------------------------------------------------------------------------------
float LR_Scheduler(float LR0, int epoch, int Nepochs)
{
    float par[3] = {-0.01, 0.2, 0.2};
    float x = 100. * epoch / Nepochs;
    return LR0 * exp(par[0] * x) * (par[1] + (1. - par[1]) * pow(cos(par[2] * x), 2));
}

float Compute_Selectivity(int level, int mode)
{return 0;}

// Calculate selectivity of set of neurons
// ---------------------------------------
/*

float Compute_Selectivity(int level, int mode)
{
    float S = 0.;
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
            float maxeff = 0.;
            float sumeff = 0.;
            for (int ic = 0; ic < N_classes; ic++)
            {
                float e = Eff[ic + N_classes * in];
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
            float maxeff = 0.;
            float sumeff = 0.;
            for (int ic = 0; ic < N_classes; ic++)
            {
                float e = Eff[ic + N_classes * in];
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
        float Effn[MaxNeurons];
        float Effc[MaxClasses];
        float sumeff = 0.;
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
*/
// Compute Q-value
// ---------------
float Compute_Q(float eff, float acc, float sel)
{
    // efficiency saturates at eff_target
    if (eff > eff_target)
        eff = eff_target + pow(eff - eff_target, 1.5);
    // acceptance saturates at acc_target
    if (acc < acc_target)
        acc = acc_target - pow(acc_target - acc, 1.5);
    float Q0 = eff / sqrt(acc);
    float w = 5. * (exp(Q0 / 4.) - 1.) / (exp(1.) - 1.);
    return Q0 + w * sel;
}

// To read our preprocessed file
void ReadFromProcessed(TTree *IT, TTree *OT, long int id_event_value)
{
    Reset_hits();
    pclass = 0;
    N_part = 0;
    First_angle = max_angle;

    float z;
    float r, phi;
    float id_event;
    float type;
    float cluster_pclass;

    IT->SetBranchAddress("cluster_z", &z);
    IT->SetBranchAddress("cluster_R", &r);
    IT->SetBranchAddress("cluster_phi", &phi);
    IT->SetBranchAddress("eventID", &id_event);
    IT->SetBranchAddress("cluster_type", &type);
    IT->SetBranchAddress("pclass", &cluster_pclass);

    // Loop over entries and find rows with the specified id_event value
    for (long int i = last_row_event_IT; i < IT->GetEntries(); ++i)
    {
        IT->GetEntry(i);

        if (static_cast<long int>(id_event) != id_event_value)
        {
            last_row_event_IT = i;
            break;
        }
        phi += M_PI;
        if (static_cast<int>(type) == 1)
        {
            type = SIG;
            pclass = (int)cluster_pclass;
            N_part = 1;
            phi += 2. * M_PI * ((int)(ievent / (NROOT))) * 1. / ((int)(N_events / NROOT) + 1);
            if (phi >= 2. * M_PI)
                phi -= 2. * M_PI;
            if (phi < First_angle)
                First_angle = phi;
        }
        else
        {
            type = BGR;
            if (phi >= 2. * M_PI)
                phi -= 2. * M_PI;
        }

        hit_pos.emplace_back(r, z, phi, static_cast<int>(type));
    }

    // OUT Tracker

    OT->SetBranchAddress("cluster_z", &z);
    OT->SetBranchAddress("cluster_R", &r);
    OT->SetBranchAddress("cluster_phi", &phi);
    OT->SetBranchAddress("eventID", &id_event);
    OT->SetBranchAddress("cluster_type", &type);

    for (long int i = last_row_event_OT; i < OT->GetEntries(); ++i)
    {
        OT->GetEntry(i);
        if (static_cast<long int>(id_event) != id_event_value)
        {
            last_row_event_OT = i;
            break;
        }
        phi += M_PI;
        if (static_cast<int>(type) == 1)
        {
            phi += 2. * M_PI * ((int)(ievent / (NROOT))) * 1. / ((int)(N_events / NROOT) + 1);
            type = SIG;
            if (phi >= 2. * M_PI)
                phi -= 2. * M_PI;
            if (phi < First_angle)
                First_angle = phi;
        }
        else
        {
            type = BGR;
            if (phi >= 2. * M_PI)
                phi -= 2. * M_PI;
        }

        hit_pos.emplace_back(r, z, phi, static_cast<int>(type));
    }

} 


void SNN_Tracking(SNN &snn_in)
{
    // Pass parameters:
    // ----------------
    // N_ev:      total number of simulated events
    // N_ep:      number of weight-learning cycles divido N_ev in N_ep gruppi e faccio girare il learning
    // NL0, NL1:  number of neurons performing track pattern recognition, organized in 2 layers
    // N_cl:      number of different signal classes (different particle momenta)
    // rootInput: name of the root file if you want to load the data. If not provided it will simulate the events
    // CF:        fraction of connected neurons between L0 and L1
    // ipspdf:    IPSP time dilation factor (to increase effect of inhibition)
    // Trainingcode: binary code to turn on updates of 9 pars, delays, voids
    // ReadPars:  whether to read parameters in from file
    // Thresh0:   threshold for firing at Layer 0
    // Thresh1:   threshold for firing at Layer 1
    // alpha:     parameter of signal model
    // l1if:      L1 inhibit factor
    // k, k1, k2: parameters of signal model
    // IEPC:      inhibition-excitation potential constant
    // IPSPdf:    IPSP dt dilation factor

    // The routine works as follows:
    // -----------------------------
    // Define tracker geometry ok
    // Initialize neuron potentials and synapse weights ok
    // Simulate background hits
    // Simulate hits from particle interactions
    // Encode hits in spike streams
    // Loop on optimization cycles (N_epochs)
    //   Loop on time steps in event-based fashion and modify neuron and synapse potentials
    //     Determine when an output spike will occur and jump there
    //     If spike:
    //       - record spike
    //       - update synapses
    //       - update neuron
    //       - inhibit other neurons
    //       - check latency with respect to signal patterns, record it
    // Dump statistics for run
    // ----------------------------------------------------------------------------------
    //
    // Every tracking layer provides an encoded stream of spikes into a dense layer of L0 neurons;
    // Every L0 neuron is densely connected to a L1 neuron, which also receives all tracking streams.
    // The topology is sketched below for 4 tracking layers, 3 L0 neurons, 2 L1 neurons
    //
    //  ------------|---- ....... -|---
    //                         O  ---|-
    //  ----|---|-------- ....... -----      O-|->
    //                     X   O  -----   X
    //  ------|----|----- ....... --|--      O--->
    //                         O  -|---
    //  --------------|-- ....... -----
    //
    //  Notes:
    //  - The optimization strategy is naive and needs to be improved
    //  - In addition, rather than re-learning weights every time parameters are modified, some better
    //    way of handling this should be implemented
    //  - IE spike model is to be improved
    // -----------------------------------------------------------------------------------------------

    // Pass parameters can't update static values, so we need to reassign the latter
    if (batch)
        gROOT->SetBatch(kTRUE);

    // initialization of neurons_index vector
    for (int i = 0; i < snn_in.N_neurons; i++)
        neurons_index.push_back(i);

    if (TrainingCode / 4 > 0)
    {
        update9 = true;
        TrainingCode -= 4;
    }
    if (TrainingCode / 2 > 0)
    {
        updateDelays = true;
        TrainingCode -= 4;
    }
    if (TrainingCode > 0)
    {
        updateConnections = true;
    }

    // Initial checks
    // --------------
    if (N_events > MaxEvents)
    {
        cout << "  Sorry, max # of events is 10,000,000. Terminating." << endl;
        return;
    }
    NevPerEpoch = N_events / N_epochs;
    if (NevPerEpoch < 1000)
    {
        cout << "  Too few events per epoch. Set to " << 10000 * N_epochs << endl;
        // N_ev = 10000*N_ep;
    }
    if (N_epochs < 1)
    {
        cout << "  Invalid N_epochs = " << N_epochs << ". Set to 1." << endl;
        N_epochs = 1;
    }
    if (snn_in.N_neuronsL[0] + snn_in.N_neuronsL[1] > MaxNeurons)
    {
        cout << "  Sorry, too many neurons. Terminating." << endl;
        return;
    }
    if (N_classes > MaxClasses)
    {
        cout << "  Sorry, too many classes (max is " << MaxClasses << "). Terminating." << endl;
        return;
    }

    // Welcome screen
    // --------------
    cout << endl
         << endl;
    cout << "                                 ------------------------------------" << endl;
    cout << endl;
    cout << "                                    S   N   N      T r a c k i n g" << endl;
    cout << endl;
    cout << "                                 ------------------------------------" << endl;
    cout << endl
         << endl
         << endl
         << endl;
    cout << "         ------------------------------------------------------------------------------------    " << endl;
    cout << "         Unsupervised search for tracks in 8-layer strip detector with spiking neural network    " << endl;
    cout << "                                                                             T.Dorigo, 3/2023    " << endl;
    cout << "         ------------------------------------------------------------------------------------    " << endl;
    cout << endl;
    cout << "         Run parameters: " << endl;
    cout << "         -----------------------------------" << endl;
    cout << "                       L0 neurons: " << snn_in.N_neuronsL[0]<< endl;
    cout << "                       L1 neurons: " << snn_in.N_neuronsL[1] << endl;
    cout << "            Connected L0-L1 frac.: " << snn_in.CF01 << endl;
    cout << "            Connected IN-L0 frac.: " << snn_in.CFI0 << endl;
    cout << "            Connected IN-L1 frac.: " << snn_in.CFI1 << endl;
    cout << "                    Track classes: " << N_classes << endl;
    cout << "                     Total events: " << N_events << endl;
    cout << "               Optimization loops: " << N_epochs << endl;
    cout << "             Optimize SNN params.: ";
    if (update9)
    {
        cout << "True" << endl;
    }
    else
    {
        cout << "False" << endl;
    }
    cout << "                  Optimize delays: ";
    if (updateDelays)
    {
        cout << "True" << endl;
    }
    else
    {
        cout << "False" << endl;
    }
    cout << "             Optimize connections: ";
    if (updateConnections)
    {
        cout << "True" << endl;
    }
    else
    {
        cout << "False" << endl;
    }
    cout << "                  Max mod. factor: " << MaxFactor << endl;

    // Suppress root warnings
    gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");
    gROOT->ProcessLine("gPrintViaErrorHandler = kTRUE;");

    // Histograms definition
    // ---------------------
    TH1F *SelectivityL0 = new TH1F("SelectivityL0", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *SelectivityL1 = new TH1F("SelectivityL1", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *Qvalue = new TH1F("Qvalue", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *Qmax = new TH1F("Qmax", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HEff = new TH1F("HEff", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HAcc = new TH1F("HAcc", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HT0 = new TH1F("HT0", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HT1 = new TH1F("HT1", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HA = new TH1F("HA", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HL1IF = new TH1F("HL1IF", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HK = new TH1F("HK", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HK1 = new TH1F("HK1", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HK2 = new TH1F("HK2", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HIEPC = new TH1F("HIEPC", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HIPSPdf = new TH1F("HIPSPdf", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH2F *EffMap = new TH2F("EffMap", "", snn_in.N_neurons, -0.5, snn_in.N_neurons - 0.5, N_classes, -0.5, N_classes - 0.5);
    SelectivityL1->SetLineColor(kBlack);
    Qmax->SetLineColor(2);
    HEff->SetMaximum(1.1);
    HEff->SetMinimum(0.);
    HAcc->SetLineColor(kRed);
    HAcc->SetMaximum(1.1);
    HAcc->SetMinimum(0.);
    HA->SetLineColor(kBlack);
    HT1->SetLineColor(kRed);
    HK1->SetLineColor(kRed);
    HK2->SetLineColor(kBlue);
    HT0->SetMinimum(0.);
    HK2->SetMinimum(0.);
    Qvalue->SetMinimum(0.);
    Qmax->SetMinimum(0.);
    HA->SetMinimum(0.);
    HIEPC->SetMinimum(0.);
    HIPSPdf->SetMinimum(0.);
    HL1IF->SetMinimum(0.);
    TH1F *HDelays = new TH1F("HDelays", "", 50, 0., snn_in.MaxDelay);
    TH2F *HVoidWs = new TH2F("HVoidWs", "", snn_in.N_neurons, -0.5, -0.5 + snn_in.N_neurons, snn_in.N_streams, -0.5, -0.5 + snn_in.N_streams);
    TH2F *Q_12 = new TH2F("Q_12", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *Q_34 = new TH2F("Q_34", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *Q_56 = new TH2F("Q_56", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *Q_78 = new TH2F("Q_78", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *Q_93 = new TH2F("Q_93", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *Q_MV = new TH2F("Q_MV", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *N_12 = new TH2F("N_12", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *N_34 = new TH2F("N_34", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *N_56 = new TH2F("N_56", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *N_78 = new TH2F("N_78", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *N_93 = new TH2F("N_93", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F *N_MV = new TH2F("N_MV", "", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);

    int N_bins = 100;
    TH2F *Latency[MaxNeurons * MaxClasses];
    char name[50];
    for (int i = 0; i < snn_in.N_neurons * N_classes; i++)
    {
        sprintf(name, "Latency%d", i);
        Latency[i] = new TH2F(name, name, N_bins, 0., (float)NevPerEpoch, max_angle + Empty_buffer, 0., (max_angle + Empty_buffer) / omega);
    }
    TH1F *HWeight[MaxNeurons * MaxStreams];
    TH1F *HRMSWeight[MaxNeurons];
    TH1F *HMaxWeight[MaxNeurons];
    TH1F *HMinWeight[MaxNeurons];
    TH1F *Efficiency[MaxNeurons * MaxClasses];
    TH1F *FakeRate[MaxNeurons];
    TH1F *Eff_totL0[MaxClasses];
    TH1F *Eff_totL1[MaxClasses];
    TH2F *StreamsS[10];
    TH2F *StreamsB[10];
    TH2F *StreamsN[10];
    TH1F *BestEff[MaxNeurons];
    TH1F *BestFR[MaxNeurons];
    TH1F *BestEtot[MaxNeurons];

    for (int i = 0; i < snn_in.N_neurons * snn_in.N_streams; i++)
    {
        sprintf(name, "HWeight%d", i);
        HWeight[i] = new TH1F(name, name, N_bins, 0., (float)NevPerEpoch);
    }
    for (int i = 0; i < snn_in.N_neurons; i++)
    {
        sprintf(name, "HRMSWeight%d", i);
        HRMSWeight[i] = new TH1F(name, name, N_bins, 0., (float)NevPerEpoch);
    }
    for (int i = 0; i < snn_in.N_neurons; i++)
    {
        sprintf(name, "HRMSWeight%d", i);
        HMaxWeight[i] = new TH1F(name, name, N_bins, 0., (float)NevPerEpoch);
    }
    for (int i = 0; i < snn_in.N_neurons; i++)
    {
        sprintf(name, "HRMSWeight%d", i);
        HMinWeight[i] = new TH1F(name, name, N_bins, 0., (float)NevPerEpoch);
    }
    for (int i = 0; i < snn_in.N_neurons * N_classes; i++)
    {
        sprintf(name, "Efficiency%d", i);
        Efficiency[i] = new TH1F(name, name, N_epochs, 0.5, 0.5 + N_epochs);
    }
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        sprintf(name, "FakeRate%d", in);
        FakeRate[in] = new TH1F(name, name, N_epochs, 0.5, 0.5 + N_epochs);
    }
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        sprintf(name, "BestEff%d", in);
        BestEff[in] = new TH1F(name, name, N_classes, -0.5, -0.5 + N_classes);
        sprintf(name, "BestFR%d", in);
        BestFR[in] = new TH1F(name, name, N_classes, -0.5, -0.5 + N_classes);
        sprintf(name, "BestEtot%d", in);
        BestEtot[in] = new TH1F(name, name, N_classes, -0.5, -0.5 + N_classes);
    }
    for (int ic = 0; ic < N_classes; ic++)
    {
        sprintf(name, "Eff_totL0%d", ic);
        Eff_totL0[ic] = new TH1F(name, name, N_epochs, 0.5, 0.5 + N_epochs);
        sprintf(name, "Eff_totL1%d", ic);
        Eff_totL1[ic] = new TH1F(name, name, N_epochs, 0.5, 0.5 + N_epochs);
    }

    for (int i = 0; i < 10; i++)
    {
        sprintf(name, "StreamsS%d", i);
        StreamsS[i] = new TH2F(name, name, (max_angle + Empty_buffer) * 500, 0., (max_angle + Empty_buffer) * 50. / omega, snn_in.N_InputStreams + snn_in.N_neurons, 0.5, snn_in.N_InputStreams + snn_in.N_neurons + 0.5);
        sprintf(name, "StreamsB%d", i);
        StreamsB[i] = new TH2F(name, name, (max_angle + Empty_buffer) * 500, 0., (max_angle + Empty_buffer) * 50. / omega, snn_in.N_InputStreams + snn_in.N_neurons, 0.5, snn_in.N_InputStreams + snn_in.N_neurons + 0.5);
        sprintf(name, "StreamsN%d", i);
        StreamsN[i] = new TH2F(name, name, (max_angle + Empty_buffer) * 500, 0., (max_angle + Empty_buffer) * 50. / omega, snn_in.N_neurons, 0.5, snn_in.N_neurons + 0.5);
    }

    // If requested, read in parameters
/*


    int okfile;
    if (ReadPars)
    {
        okfile = Read_Parameters();
        if (okfile == -1)
            return;
    }
*/
    // Final part of initial printout
    cout << "         -----------------------------------" << endl;
    cout << endl;
    cout << "         Starting values of parameters:" << endl;
    cout << "         -----------------------------------" << endl;
    cout << "                     L0 threshold: " << snn_in.Threshold[0] << endl;
    cout << "                     L1 threshold: " << snn_in.Threshold[1] << endl;
    cout << "                            alpha: " << snn_in.alpha << endl;
    cout << "                        L1inhibit: " << snn_in.L1inhibitfactor << endl;
    cout << "                                K: " << snn_in.K << endl;
    cout << "                               K1: " << snn_in.K1 << endl;
    cout << "                               K2: " << snn_in.K2 << endl;
    cout << "                 IE pot. constant: " << snn_in.IE_Pot_const << endl;
    cout << "                 IPSP dt dilation: " << snn_in.IPSP_dt_dilation << endl;
    cout << "         -----------------------------------" << endl;
    cout << endl;


    // Prime the event loop - we continuously sample detector readout and feed inputs to synapses
    // ------------------------------------------------------------------------------------------
    // TODO: MaxNeurons -> snn_in.N_neurons
    int N_fires[MaxNeurons];
    float LastP[MaxNeurons];
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        N_fires[in] = 0.;
        LastP[in] = 0.;
    }
    int fired_sum[MaxClasses][MaxNeurons];
    int random_fire[MaxNeurons];
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        random_fire[in] = 0;
        for (int ic = 0; ic < N_classes; ic++)
        {
            fired_sum[ic][in] = 0;
        }
    }
    bool not_fired_bgr = true;
    int atleastonefired = 0;
    int gen_sum[MaxClasses];
    int fired_anyL0[MaxClasses];
    int fired_anyL1[MaxClasses];
    for (int ic = 0; ic < N_classes; ic++)
    {
        gen_sum[ic] = 0;
        fired_anyL0[ic] = 0;
        fired_anyL1[ic] = 0;
    }
    bool doneL0[MaxClasses];
    bool doneL1[MaxClasses];
    bool Seen[MaxClasses][MaxNeurons];
    float selectivityL0 = 0.;
    float selectivityL1 = 0.;
    float averefftotL1 = 0.;
    float averacctotL1 = 0.;

    //TODO: implement clone method
    SNN snn_best(_NL0,  _NL1,
          _alpha,
          _CFI0, _CFI1, _CF01,
          _L1inhibitfactor,
          _K, _K1, _K2,
          _IE_Pot_const, _IPSP_dt_dilation,
          _MaxDelay,

          _tau_m,  _tau_s,  _tau_r,  _tau_plus,  _tau_minus,
          _a_plus,  _a_minus,

          _N_InputStreams,
          _Threshold0,  _Threshold1);
    SNN snn_old(_NL0,  _NL1,
          _alpha,
          _CFI0, _CFI1, _CF01,
          _L1inhibitfactor,
          _K, _K1, _K2,
          _IE_Pot_const, _IPSP_dt_dilation,
          _MaxDelay,

          _tau_m,  _tau_s,  _tau_r,  _tau_plus,  _tau_minus,
          _a_plus,  _a_minus,

          _N_InputStreams,
          _Threshold0,  _Threshold1);
        
    // Storing parameters subjected to random search
    snn_old.Threshold[0] = snn_in.Threshold[0];
    snn_old.Threshold[1] = snn_in.Threshold[1];
    snn_old.alpha = snn_in.alpha;
    snn_old.L1inhibitfactor = snn_in.L1inhibitfactor;
    snn_old.K = snn_in.K;
    snn_old.K1 = snn_in.K1;
    snn_old.K2 = snn_in.K2;
    snn_old.IE_Pot_const = snn_in.IE_Pot_const;
    snn_old.IPSP_dt_dilation = snn_in.IPSP_dt_dilation;
    float Weight_initial[snn_in.N_neurons][snn_in.N_streams];

    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        for (int is = 0; is < snn_in.N_streams; is++)
        {
            snn_old.Delay[in][is] = snn_in.Delay[in][is];
            snn_best.Delay[in][is] = snn_in.Delay[in][is];
        }
    }

    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        for (int is = 0; is < snn_in.N_streams; is++)
        {
            snn_old.Void_weight[in][is] = snn_in.Void_weight[in][is];
            snn_best.Void_weight[in][is] = snn_in.Void_weight[in][is];

            snn_old.Weight[in][is] = snn_in.Weight[in][is];
            snn_best.Weight[in][is] = snn_in.Weight[in][is];
            Weight_initial[in][is] = snn_in.Weight[in][is];
        }

    }
    float Q = 0.;
    float Q_old = 0.;
    Q_best = 0.;
    SelL1_best = 0.;
    Eff_best = 0.;
    Acc_best = 0.;
    
    snn_best.Threshold[0] = 0.;
    snn_best.Threshold[1] = 0.;
    snn_best.alpha = 0.;
    snn_best.L1inhibitfactor = 0.;
    snn_best.K = 0.;
    snn_best.K1 = 0.;
    snn_best.K2 = 0.;
    snn_best.IE_Pot_const = 0.;
    snn_best.IPSP_dt_dilation = 0.;

    float Optvar[9];
    float max_dx[9];
    float aver_dQ[9];
    Optvar[0] = 0.; // Threshold[0];
    Optvar[1] = 0.; // Threshold[1];
    Optvar[2] = 0.; // alpha;
    Optvar[3] = 0.; // L1inhibitfactor;
    Optvar[4] = 0.; // K;
    Optvar[5] = 0.; // K1;
    Optvar[6] = 0.; // K2;
    Optvar[7] = 0.; // IE_Pot_const;
    Optvar[8] = 0.; // IPSP_dt_dilation;
    for (int i = 0; i < 9; i++)
    {
        aver_dQ[i] = 0.;
        max_dx[i] = MaxFactor; // max factor of change in parameter values during optimization
    }
    float OptvarD[MaxNeurons * MaxStreams];
    float max_dxD[MaxNeurons * MaxStreams];
    float aver_dQD[MaxNeurons * MaxStreams];
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        for (int is = 0; is < snn_in.N_streams; is++)
        {
            int id = in * snn_in.N_streams + is;
            aver_dQD[id] = 0.;
            max_dxD[id] = 0.01; // max delay change factor during optimization
            OptvarD[id] = snn_in.Delay[in][is];
        }
    }
    int MaxdQHist = 50; // we do not want the full history, because we are moving in the par space; only last 10 points.
    vector<float> dQHist;
    vector<float> OptvarHist[9];
    vector<float> OptvarDHist[MaxNeurons * MaxStreams];
    int ibad = 0;
    float LR = MaxFactor;

    // Big loop on events
    // ------------------
    bool doprogress = true;
    int block = N_events / N_epochs / 50;
    if (block < 1)
        doprogress = false;
    if (doprogress)
        cout << "         " << progress[0];
    int currchar = 1;
    ievent = 0;
    int iev_thisepoch = 0;
    int iepoch = 0;
    int ind_qbest = 0;

    // Read the file with True Events and Generated BKG ------------------
    TFile *file = TFile::Open(rootInput, "READ");
    if (!file || file->IsZombie())
    {
        cerr << "Error: Cannot open file " << rootInput << endl;
        return;
    }

    TDirectoryFile *dirIT = dynamic_cast<TDirectoryFile *>(file->Get("clusterValidIT"));
    TDirectoryFile *dirOT = dynamic_cast<TDirectoryFile *>(file->Get("clusterValidOT"));

    if (!dirIT)
    {
        cerr << "Error: Cannot access directory clusterValidIT" << endl;
        file->Close();
        return;
    }

    if (!dirOT)
    {
        cerr << "Error: Cannot access directory clusterValidOT" << endl;
        file->Close();
        return;
    }

    TTree *IT = dynamic_cast<TTree *>(dirIT->Get("tree"));
    TTree *OT = dynamic_cast<TTree *>(dirOT->Get("tree"));

    if (!IT)
    {
        cerr << "Error: Cannot access tree in clusterValidIT" << endl;
        file->Close();
        return;
    }

    if (!OT)
    {
        cerr << "Error: Cannot access tree in clusterValidOT" << endl;
        file->Close();
        return;
    }

    IT->SetMaxVirtualSize(250000000);
    IT->LoadBaskets();

    OT->SetMaxVirtualSize(250000000);
    OT->LoadBaskets();

    // End of reading ----------------------------------------------
    
    
    // Create csv fout file
    ofstream fout;
    char csv_name[80];
    sprintf(csv_name, "MODE/CSV/NL0=%d_NL1=%d_NCl=%d_CF01=%.2f_CFI0=%.2f_CFI1=%.2f_alfa=%.2f_output.csv", snn_in.N_neuronsL[0], snn_in.N_neuronsL[1], N_classes, snn_in.CF01, snn_in.CFI0, snn_in.CFI1, snn_in.alpha);
    fout.open(csv_name);
    fout << "Event,ID,Stream,Time,Pclass" << endl;

    // Loop on events ----------------------------------------------
    bool insert = true;
    do
    {
        iev_thisepoch++;
        // if (ievent == 19531)
            // ievent = 19532;

        // if(ievent > 19000)
            // cout << "Event: " << ievent << endl;

        if(ievent % 1000 == 0)
            cout<<"Event: "<<ievent<<endl;
        
        if (doprogress)
        {
            if (ievent % block == 0)
            {
                // cout << "." << progress[currchar];
                // currchar++;
            }
        }

        if (ievent % NROOT == 0)
        {
            last_row_event_IT = 0;
            last_row_event_OT = 0;
        }

        ReadFromProcessed(IT, OT, ievent % NROOT);
        
        // See if we find with track with positive latency by at least one neuron
        for (int in = 0; in < snn_in.N_neurons; in++)
        {
            Seen[pclass][in] = false; // Becomes true if the neuron in has fired for the class pclass
        }
        doneL0[pclass] = false;
        doneL1[pclass] = false;
        not_fired_bgr = true;

        // Encode hits in spike streams
        // Here we encode the position of hits through the timing of a spike,
        // In the future we might think at how the analog charge readout in the hits could also be added
        float previous_firetime = 0;
        PreSpike_Time.clear();
        PreSpike_Stream.clear();
        PreSpike_Signal.clear();
        
        float t_in = ievent * (max_angle + Empty_buffer) / omega; // Initial time -> every event adds 25 ns
        Encode(t_in);

        // Keep track of latency calc for each neuron in this event
        bool not_filled[MaxNeurons];
        for (int in = 0; in < snn_in.N_neurons; in++)
        {
            not_filled[in] = true;
        }

        // Loop on spikes and modify neuron and synapse potentials
        // -------------------------------------------------------
        for (int ispike = 0; ispike < PreSpike_Time.size(); ispike++)
        {
            // By looping to size(), we can insert along the way and still make it to the end
            float t = PreSpike_Time[ispike];

            // Modify neuron potentials based on synapse weights
            // -------------------------------------------------
            float min_fire_time = snn_in.largenumber; // if no fire, neuron_firetime returns largenumber
            int in_first = -1;

            // Loop on neurons, but not in order to not favor any neuron
            // ---------------------------------------------------------

            // Shuffle order
            auto rng = default_random_engine{};
            shuffle(neurons_index.begin(), neurons_index.end(), rng);

            for (auto in : neurons_index)
            {
                // Compute future fire times of neurons and their order
                float fire_time = snn_in.Neuron_firetime(in, t);

                if (fire_time < min_fire_time)
                {
                    in_first = in;
                    min_fire_time = fire_time;
                }
            }
            if (in_first == -1)
                insert = true;

            // Ok, neuron in_first is going to fire next.
            // Peek at next event in list, to see if it comes before in_first fires
            // --------------------------------------------------------------------
            else
            {
                float latency = 0.;
                N_fires[in_first]++;
                snn_in.Fire_time[in_first].push_back(min_fire_time);

                // Learn weights with spike-time-dependent plasticity: long-term synaptic potentiation
                snn_in.LTP(in_first, min_fire_time, nearest_spike_approx, snn_old);

                // Reset history of this neuron
                snn_in.History_time[in_first].clear();
                snn_in.History_type[in_first].clear();
                snn_in.History_ID[in_first].clear();
                snn_in.History_time[in_first].push_back(min_fire_time);
                snn_in.History_type[in_first].push_back(0);
                snn_in.History_ID[in_first].push_back(0); // ID is not used for type 0 history events

                // IPSP for all others at relevant layer
                for (int in2 = 0; in2 < snn_in.N_neurons; in2++)
                {
                    if (in2 != in_first)
                    {
                        if (snn_in.Neuron_layer[in2] == snn_in.Neuron_layer[in_first])
                        { // inhibitions within layer or across
                            snn_in.History_time[in2].push_back(min_fire_time);
                            snn_in.History_type[in2].push_back(2);
                            snn_in.History_ID[in2].push_back(in_first);
                        }
                    }
                }
                // Create EPS signal in L0 neuron-originated streams
                if (snn_in.Neuron_layer[in_first] == 0)
                {   // this is a Layer-0 neuron
                    for(int in = snn_in.N_neuronsL[0]; in < snn_in.N_neurons; in++){
                        snn_in.History_time[in].push_back(min_fire_time);
                        snn_in.History_type[in].push_back(1);
                        snn_in.History_ID[in].push_back(snn_in.N_InputStreams + in_first);
                        snn_in.LTD(in, in_first, min_fire_time, nearest_spike_approx, snn_old); 

                    }
                }
                
                // Fill spikes train histogram
                if (ievent >= N_events - 500.)
                {
                    int is = (ievent - N_events + 500) / 50;
                    float time = min_fire_time - (max_angle + Empty_buffer) / omega * (ievent / 50) * 50;
                    StreamsN[is]->Fill(time, in_first + 1);
                    if (N_part > 0)
                    {
                        fout << ievent << ", " << 2 << ", " << in_first << "," << time << "," << pclass << endl;
                    }
                }

                // Fill latency histogram
                if (N_part > 0)
                {
                    // How long did it take the first neuron to fire with respect to the arrival time of the first hit?
                    latency = min_fire_time - t_in - First_angle / omega;
                    if (latency >= 0. && not_filled[in_first])
                    {
                        if (iepoch == N_epochs - 1)
                            Latency[in_first * N_classes + pclass]->Fill(0.5 + iev_thisepoch, latency);
                        Seen[pclass][in_first] = true;
                        not_filled[in_first] = false;
                    }
                }
                else
                {
                    if (not_filled[in_first] && iev_thisepoch > NevPerEpoch * 0.9)
                    {
                        random_fire[in_first]++;
                        not_filled[in_first] = false;
                    }
                    if (in_first >= snn_in.N_neuronsL[0] && iev_thisepoch > NevPerEpoch * 0.9)
                    { // for Q-value calculations
                        if (not_fired_bgr)
                        {
                            atleastonefired++;
                            not_fired_bgr = false;
                        }
                    }
                }

                ispike-=1;
                insert=false;

                //take a step back and search for another activation
                previous_firetime = min_fire_time;

            }// end if in_first fires
            //insert the new spike for the next iteration
            if(insert){
                // Save information on hit-based streams for last 500 events to histograms
                if (ievent >= N_events - 500.)
                {
                    // dividing N_events in 10 groups
                    int is = (ievent - N_events + 500) / 50;
                    // time = tin + thit - tin(First event of the group)
                    float time = PreSpike_Time[ispike] - (max_angle + Empty_buffer) / omega * (ievent / 50) * 50;

                    // Histograms
                    if (PreSpike_Signal[ispike] == 1)
                    {
                        StreamsS[is]->Fill(time, PreSpike_Stream[ispike] + 1);
                        fout << ievent << ", " << PreSpike_Signal[ispike] << ", " << PreSpike_Stream[ispike] + 1 << "," << time << "," << pclass << endl;
                    }
                    else if (PreSpike_Signal[ispike] == 0)
                    {
                        StreamsB[is]->Fill(time, PreSpike_Stream[ispike] + 1);
                    }
                }
                int is = PreSpike_Stream[ispike];
                for (auto in : neurons_index)
                {
                    //  We implement a scheme where input streams produce an IE signal into L0, an EPS into L1, and L0 neurons EPS into L1
                    //  Add to neuron history, masking out L1 spikes for L0 neurons
                    if (!snn_in.Void_weight[in][is])
                    { // otherwise stream "is" does not lead to neuron "in"
                            snn_in.History_time[in].push_back(t);
                            // All input spikes lead to EPSP
                            snn_in.History_type[in].push_back(1);
                            snn_in.History_ID[in].push_back(is);
                            snn_in.LTD(in, is, t, nearest_spike_approx, snn_old); 
                    }
                }
            }
        }     // end ispike loop, ready to start over

        // Fill info for efficiency calculations
        if (N_part > 0)
        {
            // efficiency calculations only in the last 10% of the epoch
            if (iev_thisepoch > NevPerEpoch * 0.9)
            {
                gen_sum[pclass]++;
                for (int in = 0; in < snn_in.N_neurons; in++)
                {
                    if (Seen[pclass][in])
                    {
                        fired_sum[pclass][in]++;
                        if (in < snn_in.N_neuronsL[0])
                        {
                            if (!doneL0[pclass])
                            {
                                doneL0[pclass] = true;
                                fired_anyL0[pclass]++;
                            }
                        }
                        else
                        {
                            if (!doneL1[pclass])
                            {
                                doneL1[pclass] = true;
                                fired_anyL1[pclass]++;
                            }
                        }
                    }
                }
            }
        }

        // Write histograms of weights
        if (iepoch == N_epochs - 1)
        {
            int bin = (int)(1000. * (float)iev_thisepoch / NevPerEpoch);
            for (int in = 0; in < snn_in.N_neurons; in++)
            {
                SumofSquaresofWeight[in] = 0;
                MaxWeight[in] = -0.1; // this for the purpose of finding the maxima
                MinWeight[in] = 1.1;  // for the purpose of finding the minima

                for (int is = 0; is < snn_in.N_streams; is++)
                {
                    // int bin = (int)(1000. * (float)iev_thisepoch / NevPerEpoch);
                    if (!snn_in.Void_weight[in][is])
                    {
                        HWeight[in * snn_in.N_streams + is]->SetBinContent(bin, snn_in.Weight[in][is]);
                        SumofSquaresofWeight[in] += snn_in.Weight[in][is] * snn_in.Weight[in][is]; // for RMS calculation
                        if (snn_in.Weight[in][is] > MaxWeight[in])
                            MaxWeight[in] = snn_in.Weight[in][is]; // finding maxima
                        else if (snn_in.Weight[in][is] < MinWeight[in])
                            MinWeight[in] = snn_in.Weight[in][is]; // finding minima
                    }
                    else
                        HWeight[in * snn_in.N_streams + is]->SetBinContent(bin, -1);
                }

                // RMS Plot
                MeanofSquaresofWeight[in] = SumofSquaresofWeight[in] / (snn_in.N_neurons);
                RMSWeight[in] = sqrt(MeanofSquaresofWeight[in] - 1. / (snn_in.N_streams * snn_in.N_streams));
                HRMSWeight[in]->SetBinContent(bin, RMSWeight[in]);

                // MaxWeight plot
                HMaxWeight[in]->SetBinContent(bin, MaxWeight[in]);

                // MinWeight Plot
                HMinWeight[in]->SetBinContent(bin, MinWeight[in]);
            }
        }

        // prespike_time.push_back(time) -> time associated to an hit or to a spike coming from L0
        // prespike_Stream -> stream id associated to the hit bin or to the L0 neuron
        // prespike_Signal -> spike type: 0 if BKG, 1 if Track, 2 if NeuronFire

        // Fill efficiency histograms every NevPerEpoch events, compute Q value and Selectivity, modify parameters
        // ---------------------------------------------------------------------------------------------------
        
        if (iev_thisepoch == NevPerEpoch)
        { // we did NevPerEpoch events

            // Reset counter that inhibits efficiency and Q calculations until we reach steady state with weights
            iev_thisepoch = 0;
            iepoch++;
            // End of progress bar
            if (doprogress)
                cout << progress[51] << endl;

            for (int in = 0; in < snn_in.N_neurons; in++)
            {
                for (int ic = 0; ic < N_classes; ic++)
                {
                    int combind = ic + N_classes * in;
                    Eff[combind] = fired_sum[ic][in];
                    if (gen_sum[ic] > 0)
                        Eff[combind] /= gen_sum[ic];
                    Efficiency[combind]->SetBinContent(iepoch, Eff[combind]);
                }
                float fakerate = random_fire[in] * 2. / NevPerEpoch; // there are NevPerEpoch/2 events with no tracks, where we compute random_fire per neuron
                FakeRate[in]->SetBinContent(iepoch, fakerate);
            }
            float Efftot[MaxClasses];
            for (int ic = 0; ic < N_classes; ic++)
            {
                float etl0 = fired_anyL0[ic];
                if (gen_sum[ic] > 0)
                    etl0 /= gen_sum[ic];
                Eff_totL0[ic]->SetBinContent(iepoch, etl0);
                float etl1 = fired_anyL1[ic]; // L1 efficiency is what counts.
                if (gen_sum[ic] > 0)
                    etl1 /= gen_sum[ic];
                Eff_totL1[ic]->SetBinContent(iepoch, etl1);
                Efftot[ic] = etl1;
            }

            selectivityL0 = Compute_Selectivity(0, 2);
            SelectivityL0->Fill(iepoch, selectivityL0);
            selectivityL1 = Compute_Selectivity(1, 2);
            SelectivityL1->Fill(iepoch, selectivityL1);

            // Q value is average efficiency divided by sqrt (aver eff plus aver acceptance)
            // -----------------------------------------------------------------------------
            averacctotL1 = atleastonefired * (2. / NevPerEpoch * 10.); // total acceptance, computed with 0.1*NevPerEpoch/2 events with no tracks
            averefftotL1 = 0.;
            for (int ic = 0; ic < N_classes; ic++)
            {
                averefftotL1 += Efftot[ic];
            }
            averefftotL1 /= N_classes;

            Q = Compute_Q(averefftotL1, averacctotL1, selectivityL1);

            // Fix maximum excursion of parameters with a schedule
            LR = LR_Scheduler(MaxFactor, iepoch, N_epochs);
            for (int i = 0; i < 9; i++)
            {
                max_dx[i] = LR;
            }
            for (int id = 0; id < snn_in.N_neurons * snn_in.N_streams; id++)
            {
                max_dxD[id] = 0.1 * LR;
            }

            // Re-initialize neurons
            snn_in.Init_neurons();
            // Reset hits
            Reset_hits();
            // Reset weights to initial conditions before new investigation
            // snn_in.Reset_weights(Weight_initial);
            // Init delays
            if (!updateDelays && !ReadPars && !learnDelays)
                snn_in.Init_delays(); // This unlike void connections, because we can opt to learn these at each cycle too

            cout << "         Ev. # " << ievent + 1 << " - LR = " << LR << "; Selectivity L0 = " << selectivityL0 << " L1 = " << selectivityL1
                 << "; Eff = " << averefftotL1 << " Acc = " << averacctotL1 << "; Firings: ";

            for (int in = 0; in < snn_in.N_neurons; in++)
            {
                cout << N_fires[in] << " ";
            }
            cout << endl;

            // Keep a history of recent delta Q values, so that we know how to sample next
            dQHist.push_back(Q);
            for (int i = 0; i < 9; i++)
            {
                OptvarHist[i].push_back(Optvar[i]);
            }
            for (int id = 0; id < snn_in.N_neurons * snn_in.N_streams; id++)
            {
                OptvarDHist[id].push_back(OptvarD[id]);
            }

            // Keep only the last MaxdQHist values of dQ history
            if (dQHist.size() > MaxdQHist)
            {
                dQHist.erase(dQHist.begin());
                for (int i = 0; i < 9; i++)
                {
                    OptvarHist[i].erase(OptvarHist[i].begin());
                }
                for (int id = 0; id < snn_in.N_neurons * snn_in.N_streams; id++)
                {
                    OptvarDHist[id].erase(OptvarDHist[id].begin());
                }
            }
            // Find running weighted average of dQ-values in par space
            float dQ_max = 0.;
            for (int i = 0; i < 9; i++)
            {
                float sum_dQdx = 0.;
                float sum_dx = 0.;
                for (int j = 1; j < dQHist.size(); j++)
                {
                    float dx = OptvarHist[i][j] - OptvarHist[i][j - 1];
                    sum_dQdx += (dQHist[j] - dQHist[j - 1]) * dx;
                    sum_dx += dx;
                }
                aver_dQ[i] = 0.;
                if (sum_dx != 0.)
                    aver_dQ[i] = sum_dQdx / sum_dx;
                if (fabs(aver_dQ[i]) > dQ_max)
                    dQ_max = fabs(aver_dQ[i]);
            }
            // Same, for delays
            float dQ_maxD = 0.;
            for (int id = 0; id < snn_in.N_neurons * snn_in.N_streams; id++)
            {
                float sum_dQdxD = 0.;
                float sum_dxD = 0.;
                for (int j = 1; j < dQHist.size(); j++)
                {
                    float dxD = OptvarDHist[id][j] - OptvarDHist[id][j - 1];
                    sum_dQdxD += (dQHist[j] - dQHist[j - 1]) * dxD;
                    sum_dxD += dxD;
                }
                aver_dQD[id] = 0.;
                if (sum_dxD != 0.)
                    aver_dQD[id] = sum_dQdxD / sum_dxD;
                if (fabs(aver_dQD[id]) > dQ_maxD)
                    dQ_maxD = fabs(aver_dQD[id]);
            }

            // Fill debugging graphs of Q as a function of parameters
            int ibin, jbin;
            float n, cont, newcont;
            if (iepoch > 1 || N_epochs == 1)
            { // only do it from end of second epoch onwards, as we are filling delta values
                ibin = 1 + (int)(20. * (snn_in.Threshold[0] / snn_old.Threshold[0] - 1. + MaxFactor) / (2. * MaxFactor));
                jbin = 1 + (int)(20. * (snn_in.Threshold[1] / snn_old.Threshold[1] - 1. + MaxFactor) / (2. * MaxFactor));
                if (ibin > 0 && ibin < 21 && jbin > 0 && jbin < 21)
                {
                    cont = Q_12->GetBinContent(ibin, jbin);
                    n = N_12->GetBinContent(ibin, jbin);
                    newcont = (Q - Q_old) / (n + 1) + cont * n / (n + 1);
                    Q_12->SetBinContent(ibin, jbin, newcont);
                    N_12->SetBinContent(ibin, jbin, n + 1);
                }
                ibin = 1 + (int)(20. * (snn_in.alpha / snn_old.alpha - 1. + MaxFactor) / (2 * MaxFactor));
                jbin = 1 + (int)(20. * (snn_in.L1inhibitfactor / snn_old.L1inhibitfactor - 1. + MaxFactor) / (2. * MaxFactor));
                if (ibin > 0 && ibin < 21 && jbin > 0 && jbin < 21)
                {
                    cont = Q_34->GetBinContent(ibin, jbin);
                    n = N_34->GetBinContent(ibin, jbin);
                    newcont = (Q - Q_old) / (n + 1) + cont * n / (n + 1);
                    Q_34->SetBinContent(ibin, jbin, newcont);
                    N_34->SetBinContent(ibin, jbin, n + 1);
                }
                ibin = 1 + (int)(20. * (snn_in.K / snn_old.K - 1. + MaxFactor) / (2. * MaxFactor));
                jbin = 1 + (int)(20. * (snn_in.K1 / snn_old.K1 - 1. + MaxFactor) / (2. * MaxFactor));
                if (ibin > 0 && ibin < 21 && jbin > 0 && jbin < 21)
                {
                    cont = Q_56->GetBinContent(ibin, jbin);
                    n = N_56->GetBinContent(ibin, jbin);
                    newcont = (Q - Q_old) / (n + 1) + cont * n / (n + 1);
                    Q_56->SetBinContent(ibin, jbin, newcont);
                    N_56->SetBinContent(ibin, jbin, n + 1);
                }
                ibin = 1 + (int)(20. * (snn_in.K2 / snn_old.K2 - 1. + MaxFactor) / (2. * MaxFactor));
                jbin = 1 + (int)(20. * (snn_in.IE_Pot_const / snn_old.IE_Pot_const - 1. + MaxFactor) / (2. * MaxFactor));
                if (ibin > 0 && ibin < 21 && jbin > 0 && jbin < 21)
                {
                    cont = Q_78->GetBinContent(ibin, jbin);
                    n = N_78->GetBinContent(ibin, jbin);
                    newcont = (Q - Q_old) / (n + 1) + cont * n / (n + 1);
                    Q_78->SetBinContent(ibin, jbin, newcont);
                    N_78->SetBinContent(ibin, jbin, n + 1);
                }
                jbin = 1 + (int)(20. * (snn_in.alpha / snn_old.alpha - 1. + MaxFactor) / (2. * MaxFactor));
                ibin = 1 + (int)(20. * (snn_in.IPSP_dt_dilation / snn_old.IPSP_dt_dilation - 1. + MaxFactor) / (2. * MaxFactor));
                if (ibin > 0 && ibin < 21 && jbin > 0 && jbin < 21)
                {
                    cont = Q_93->GetBinContent(ibin, jbin);
                    n = N_93->GetBinContent(ibin, jbin);
                    newcont = (Q - Q_old) / (n + 1) + cont * n / (n + 1);
                    Q_93->SetBinContent(ibin, jbin, newcont);
                    N_93->SetBinContent(ibin, jbin, n + 1);
                }
                // Now graph of mean vs sqm of Delay distribution
                float meanDelay = 0.;
                float sqmDelay = 0.;
                float meanOldDelay = 0.;
                float sqmOldDelay = 0.;
                for (int in = 0; in < snn_in.N_neurons; in++)
                {
                    for (int is = 0; is < snn_in.N_streams; is++)
                    {
                        meanDelay += snn_in.Delay[in][is];
                        sqmDelay += pow(snn_in.Delay[in][is], 2);
                        meanOldDelay += snn_old.Delay[in][is];
                        sqmOldDelay += pow(snn_old.Delay[in][is], 2);
                    }
                }
                meanDelay /= snn_in.N_neurons * snn_in.N_streams;
                sqmDelay = sqmDelay / (snn_in.N_neurons * snn_in.N_streams) - meanDelay * meanDelay;
                sqmDelay = sqrt(sqmDelay);
                meanOldDelay /= snn_old.N_neurons * snn_old.N_streams;
                sqmOldDelay = sqmOldDelay / (snn_old.N_neurons * snn_old.N_streams) - meanOldDelay * meanOldDelay;
                sqmOldDelay = sqrt(sqmOldDelay);
                ibin = 1 + (int)(20. * (meanDelay / (meanOldDelay + snn_in.epsilon) - 1. + MaxFactor) / (2. * MaxFactor));
                jbin = 1 + (int)(20. * (sqmDelay / (sqmOldDelay + snn_in.epsilon) - 1. + MaxFactor) / (2. * MaxFactor));
                if (ibin > 0 && ibin < 21 && jbin > 0 && jbin < 21)
                {
                    cont = Q_MV->GetBinContent(ibin, jbin);
                    n = N_MV->GetBinContent(ibin, jbin);
                    newcont = (Q - Q_old) / (n + 1) + cont * n / (n + 1);
                    Q_MV->SetBinContent(ibin, jbin, newcont);
                    N_MV->SetBinContent(ibin, jbin, n + 1);
                }
                TCanvas *QQ = new TCanvas("QQ", "", 800, 600);
                QQ->Divide(3, 2);
                QQ->cd(1);
                Q_12->Draw("COL4");
                QQ->cd(2);
                Q_34->Draw("COL4");
                QQ->cd(3);
                Q_56->Draw("COL4");
                QQ->cd(4);
                Q_78->Draw("COL4");
                QQ->cd(5);
                Q_93->Draw("COL4");
                QQ->cd(6);
                Q_MV->Draw("COL4");
                QQ->Update();
            }

            // Fill histograms with delays
            HDelays->Reset();
            HVoidWs->Reset();
            for (int in = 0; in < snn_in.N_neurons; in++)
            {
                for (int is = 0; is < snn_in.N_streams; is++)
                {
                    HDelays->Fill(snn_in.Delay[in][is]);
                    if (snn_in.Void_weight[in][is])
                    {
                        HVoidWs->SetBinContent(in + 1, is + 1, 0.);
                    }
                    else
                    {
                        HVoidWs->SetBinContent(in + 1, is + 1, 1.);
                    }
                }
            }

            // Is this Q factor not larger than before?
            cout << "         Q = " << Q << " Old = " << Q_old << " Best = " << Q_best << " ib = " << ibad;

            // Update histograms with current parameter values and optimization metrics
            Qvalue->SetBinContent(iepoch, Q);
            if (Q > Q_best)
            {
                ind_qbest = iepoch;
                Q_best = Q;
                
                SelL1_best = selectivityL1;
                Eff_best = averefftotL1;
                Acc_best = averacctotL1;
                
                snn_best.Threshold[0] = snn_in.Threshold[0];
                snn_best.Threshold[1] = snn_in.Threshold[1];
                snn_best.alpha = snn_in.alpha;
                snn_best.L1inhibitfactor = snn_in.L1inhibitfactor;
                snn_best.K = snn_in.K;
                snn_best.K1 = snn_in.K1;
                snn_best.K2 = snn_in.K2;
                snn_best.IE_Pot_const = snn_in.IE_Pot_const;
                snn_best.IPSP_dt_dilation = snn_in.IPSP_dt_dilation;
                for (int in = 0; in < snn_in.N_neurons; in++)
                {
                    for (int is = 0; is < snn_in.N_streams; is++)
                    {
                        snn_best.Delay[in][is] = snn_in.Delay[in][is];
                        snn_best.Void_weight[in][is] = snn_in.Void_weight[in][is];
                    }
                }
            }
            Qmax->SetBinContent(iepoch, Q_best);
            HEff->SetBinContent(iepoch, averefftotL1);
            HAcc->SetBinContent(iepoch, averacctotL1);
            HT0->SetBinContent(iepoch, snn_in.Threshold[0]);
            HT1->SetBinContent(iepoch, snn_in.Threshold[1]);
            HA->SetBinContent(iepoch, snn_in.alpha);
            HL1IF->SetBinContent(iepoch, snn_in.L1inhibitfactor);
            HK->SetBinContent(iepoch, snn_in.K);
            HK1->SetBinContent(iepoch, snn_in.K1);
            HK2->SetBinContent(iepoch, snn_in.K2);
            HIEPC->SetBinContent(iepoch, snn_in.IE_Pot_const);
            HIPSPdf->SetBinContent(iepoch, snn_in.IPSP_dt_dilation);
            TCanvas *CU = new TCanvas("CU", "", 1600, 700);
            CU->Divide(5, 2);
            CU->cd(1);
            Qvalue->Draw();
            Qmax->Draw("SAME");
            CU->cd(2);
            HEff->Draw();
            HAcc->Draw("SAME");
            float h;
            float hmax = -snn_in.largenumber;
            for (int ibin = 1; ibin <= N_epochs; ibin++)
            {
                h = HT0->GetBinContent(ibin);
                if (hmax < h)
                    hmax = h;
                h = HT1->GetBinContent(ibin);
                if (hmax < h)
                    hmax = h;
            }
            HT0->SetMaximum(hmax + 0.1 * fabs(hmax));
            HT1->SetMaximum(hmax + 0.1 * fabs(hmax));
            CU->cd(3);
            HT0->Draw();
            HT1->Draw("SAME");
            HA->Draw("SAME");
            CU->cd(4);
            HL1IF->Draw();
            hmax = -snn_in.largenumber;
            for (int ibin = 1; ibin <= N_epochs; ibin++)
            {
                h = HK->GetBinContent(ibin);
                if (hmax < h)
                    hmax = h;
                h = HK1->GetBinContent(ibin);
                if (hmax < h)
                    hmax = h;
                h = HK2->GetBinContent(ibin);
                if (hmax < h)
                    hmax = h;
            }
            HK->SetMaximum(hmax + 0.1 * fabs(hmax));
            HK1->SetMaximum(hmax + 0.1 * fabs(hmax));
            HK2->SetMaximum(hmax + 0.1 * fabs(hmax));
            CU->cd(5);
            HK2->Draw();
            HK->Draw("SAME");
            HK1->Draw("SAME");
            CU->cd(6);
            HIEPC->Draw();
            CU->cd(7);
            HIPSPdf->Draw();
            CU->cd(8);
            HDelays->Draw();
            CU->cd(9);
            HVoidWs->Draw("COL4");
            CU->cd(10);
            hmax = -snn_in.largenumber;
            for (int ibin = 1; ibin <= N_epochs; ibin++)
            {
                h = SelectivityL0->GetBinContent(ibin);
                if (hmax < h)
                    hmax = h;
                h = SelectivityL1->GetBinContent(ibin);
                if (hmax < h)
                    hmax = h;
            }
            SelectivityL0->SetMaximum(hmax + 0.1 * fabs(hmax));
            SelectivityL1->SetMaximum(hmax + 0.1 * fabs(hmax));
            SelectivityL1->Draw();
            SelectivityL0->Draw("SAME");
            CU->Update();

            if (iepoch == 1)
            { // The first time we modify at random the parameters

                // Store previous values
                if (update9)
                {
                    snn_old.Threshold[0] = snn_in.Threshold[0];
                    snn_old.Threshold[1] = snn_in.Threshold[1];
                    snn_old.alpha = snn_in.alpha;
                    snn_old.L1inhibitfactor = snn_in.L1inhibitfactor;
                    snn_old.K = snn_in.K;
                    snn_old.K1 = snn_in.K1;
                    snn_old.K2 = snn_in.K2;
                    snn_old.IE_Pot_const = snn_in.IE_Pot_const;
                    snn_old.IPSP_dt_dilation = snn_in.IPSP_dt_dilation;
                    for (int i = 0; i < 9; i++)
                    {
                        Optvar[i] = myRNG->Uniform(-max_dx[i], max_dx[i]);
                    }
                    snn_in.Threshold[0] *= 1. + Optvar[0];
                    snn_in.Threshold[1] *= 1. + Optvar[1];
                    snn_in.alpha *= 1. + Optvar[2];
                    snn_in.L1inhibitfactor *= 1. + Optvar[3];
                    snn_in.K *= 1. + Optvar[4];
                    snn_in.K1 *= 1. + Optvar[5];
                    snn_in.K2 *= 1. + Optvar[6];
                    snn_in.IE_Pot_const *= 1. + Optvar[7];
                    snn_in.IPSP_dt_dilation *= 1. + Optvar[8];
                }
                if (updateDelays)
                {
                    for (int in = 0; in < snn_in.N_neurons; in++)
                    {
                        for (int is = 0; is < snn_in.N_streams; is++)
                        {
                            int id = in * snn_in.N_streams + is;
                            snn_old.Delay[in][is] = snn_in.Delay[in][is];
                            OptvarD[id] = myRNG->Uniform(-max_dxD[id], max_dxD[id]);
                            snn_in.Delay[in][is] += OptvarD[id];
                            if (snn_in.Delay[in][is] > snn_in.MaxDelay)
                                snn_in.Delay[in][is] = snn_in.MaxDelay;
                            if (snn_in.Delay[in][is] < 0.)
                                snn_in.Delay[in][is] = 0.;
                        }
                    }
                }
                if (updateConnections)
                {
                    for (int in = 0; in < snn_in.N_neurons; in++)
                    {
                        for (int is = 0; is < snn_in.N_streams; is++)
                        {
                            if (in >= snn_in.N_neuronsL[0] || is < snn_in.N_InputStreams)
                            {
                                snn_old.Void_weight[in][is] = snn_in.Void_weight[in][is];
                                if (!snn_in.Void_weight[in][is])
                                {
                                    if (myRNG->Uniform() < ProbWSwitchDown)
                                    {
                                        snn_in.Void_weight[in][is] = !snn_in.Void_weight[in][is];
                                    }
                                }
                                else
                                {
                                    if (myRNG->Uniform() < ProbWSwitchUp)
                                    {
                                        snn_in.Void_weight[in][is] = !snn_in.Void_weight[in][is];
                                    }
                                }
                            }
                        }
                    }
                }
                Q_old = Q;
                ibad = 0;
            }
            else
            { // We are in second or larger epoch, can look at history of improvement for directions

                if (Q <= Q_old)
                { // We did a step in the wrong direction, need to take it back and modify at random, sampling wisely...
                  // ibad counts how many trials we make from last improved Q. After 3 unsuccessful trials we allow for a step away
                    if (ibad < 2)
                    {
                        // Store previous values
                        if (update9)
                        {
                            snn_in.Threshold[0] = snn_old.Threshold[0];
                            snn_in.Threshold[1] = snn_old.Threshold[1];
                            snn_in.alpha = snn_old.alpha;
                            snn_in.L1inhibitfactor = snn_old.L1inhibitfactor;
                            snn_in.K = snn_old.K;
                            snn_in.K1 = snn_old.K1;
                            snn_in.K2 = snn_old.K2;
                            snn_in.IE_Pot_const = snn_old.IE_Pot_const;
                            snn_in.IPSP_dt_dilation = snn_old.IPSP_dt_dilation;
                        }
                        if (updateDelays)
                        {
                            for (int in = 0; in < snn_in.N_neurons; in++)
                            {
                                for (int is = 0; is < snn_in.N_streams; is++)
                                {
                                    snn_in.Delay[in][is] = snn_old.Delay[in][is];
                                }
                            }
                        }
                        if (updateConnections)
                        {
                            for (int in = 0; in < snn_in.N_neurons; in++)
                            {
                                for (int is = 0; is < snn_in.N_streams; is++)
                                {
                                    snn_in.Void_weight[in][is] = snn_old.Void_weight[in][is];
                                }
                            }
                        }
                        // The calculations below do the following:
                        // - find a number between exp(-2.) and exp(+2.) to fix the slope of the dx cumulant distribution
                        // - find the multiplier of each parameter as a function of how much dx increases Q
                        //   The random number is converted by the function pow(r,lambda) to a number between -max_dx and max_dx
                        //   which distributes uniformly for no slope of dq vs dx, and peaky at the extrema for larger correlation
                        for (int i = 0; i < 9; i++)
                        {
                            float lambda = 1.;
                            if (dQ_max != 0.)
                                lambda = exp(2. * aver_dQ[i] / dQ_max);
                            float r = myRNG->Uniform();
                            Optvar[i] = -max_dx[i] + 2. * max_dx[i] * pow(r, lambda);
                        }
                        if (update9)
                        {
                            snn_in.Threshold[0] *= 1. + Optvar[0];
                            snn_in.Threshold[1] *= 1. + Optvar[1];
                            snn_in.alpha *= 1. + Optvar[2];
                            snn_in.L1inhibitfactor *= 1. + Optvar[3];
                            snn_in.K *= 1. + Optvar[4];
                            snn_in.K1 *= 1. + Optvar[5];
                            snn_in.K2 *= 1. + Optvar[6];
                            snn_in.IE_Pot_const *= 1. + Optvar[7];
                            snn_in.IPSP_dt_dilation *= 1. + Optvar[8];
                        }
                        if (updateDelays)
                        {
                            float lambda;
                            for (int in = 0; in < snn_in.N_neurons; in++)
                            {
                                for (int is = 0; is < snn_in.N_streams; is++)
                                {
                                    int id = in * snn_in.N_streams + is;
                                    lambda = 1.;
                                    if (dQ_maxD > 0.)
                                        lambda = exp(2. * aver_dQD[id] / dQ_maxD);
                                    float r = myRNG->Uniform();
                                    OptvarD[id] = -max_dxD[id] + 2. * max_dxD[id] * pow(r, lambda);
                                    snn_in.Delay[in][is] += OptvarD[id];
                                    if (snn_in.Delay[in][is] > snn_in.MaxDelay)
                                        snn_in.Delay[in][is] = snn_in.MaxDelay;
                                    if (snn_in.Delay[in][is] < 0.)
                                        snn_in.Delay[in][is] = 0.;
                                }
                            }
                        }
                        if (updateConnections)
                        {
                            for (int in = 0; in < snn_in.N_neurons; in++)
                            {
                                for (int is = 0; is < snn_in.N_streams; is++)
                                {
                                    if (in >= snn_in.N_neuronsL[0] || is < snn_in.N_InputStreams)
                                    {
                                        if (!snn_in.Void_weight[in][is])
                                        {
                                            if (myRNG->Uniform() < ProbWSwitchDown)
                                            {
                                                snn_in.Void_weight[in][is] = !snn_in.Void_weight[in][is];
                                            }
                                        }
                                        else
                                        {
                                            if (myRNG->Uniform() < ProbWSwitchUp)
                                            {
                                                snn_in.Void_weight[in][is] = !snn_in.Void_weight[in][is];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        ibad++;
                    }
                    else
                    { // ibad=2, need to reset

                        if (update9)
                        {
                            snn_in.Threshold[0] = snn_best.Threshold[0];
                            snn_in.Threshold[1] = snn_best.Threshold[1];
                            snn_in.alpha = snn_best.alpha;
                            snn_in.L1inhibitfactor = snn_best.L1inhibitfactor;
                            snn_in.K = snn_best.K;
                            snn_in.K1 = snn_best.K1;
                            snn_in.K2 = snn_best.K2;
                            snn_in.IE_Pot_const = snn_best.IE_Pot_const;
                            snn_in.IPSP_dt_dilation = snn_best.IPSP_dt_dilation;
                            snn_old.Threshold[0] = snn_in.Threshold[0];
                            snn_old.Threshold[1] = snn_in.Threshold[1];
                            snn_old.alpha = snn_in.alpha;
                            snn_old.L1inhibitfactor = snn_in.L1inhibitfactor;
                            snn_old.K = snn_in.K;
                            snn_old.K1 = snn_in.K1;
                            snn_old.K2 = snn_in.K2;
                            snn_old.IE_Pot_const = snn_in.IE_Pot_const;
                            snn_old.IPSP_dt_dilation = snn_in.IPSP_dt_dilation;
                            for (int i = 0; i < 9; i++)
                            {
                                Optvar[i] = myRNG->Uniform(-max_dx[i], max_dx[i]);
                            }
                            snn_in.Threshold[0] *= 1. + Optvar[0];
                            snn_in.Threshold[1] *= 1. + Optvar[1];
                            snn_in.alpha *= 1. + Optvar[2];
                            snn_in.L1inhibitfactor *= 1. + Optvar[3];
                            snn_in.K *= 1. + Optvar[4];
                            snn_in.K1 *= 1. + Optvar[5];
                            snn_in.K2 *= 1. + Optvar[6];
                            snn_in.IE_Pot_const *= 1. + Optvar[7];
                            snn_in.IPSP_dt_dilation *= 1. + Optvar[8];
                        }
                        if (updateDelays)
                        {
                            for (int in = 0; in < snn_in.N_neurons; in++)
                            {
                                for (int is = 0; is < snn_in.N_streams; is++)
                                {
                                    int id = in * snn_in.N_streams + is;
                                    snn_in.Delay[in][is] = snn_best.Delay[in][is];
                                    snn_old.Delay[in][is] = snn_in.Delay[in][is];
                                    OptvarD[id] = myRNG->Uniform(-max_dxD[id], max_dxD[id]);
                                    snn_in.Delay[in][is] += OptvarD[id];
                                    if (snn_in.Delay[in][is] > snn_in.MaxDelay)
                                        snn_in.Delay[in][is] = snn_in.MaxDelay;
                                    if (snn_in.Delay[in][is] < 0.)
                                        snn_in.Delay[in][is] = 0.;
                                }
                            }
                        }
                        if (updateConnections)
                        {
                            for (int in = 0; in < snn_in.N_neurons; in++)
                            {
                                for (int is = 0; is < snn_in.N_streams; is++)
                                {
                                    if (in >= snn_in.N_neuronsL[0] || is < snn_in.N_InputStreams)
                                    {
                                        snn_in.Void_weight[in][is] = snn_best.Void_weight[in][is];
                                        snn_old.Void_weight[in][is] = snn_in.Void_weight[in][is];
                                        if (!snn_in.Void_weight[in][is])
                                        {
                                            if (myRNG->Uniform() < ProbWSwitchDown)
                                            {
                                                snn_in.Void_weight[in][is] = !snn_in.Void_weight[in][is];
                                            }
                                        }
                                        else
                                        {
                                            if (myRNG->Uniform() < ProbWSwitchUp)
                                            {
                                                snn_in.Void_weight[in][is] = !snn_in.Void_weight[in][is];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        Q_old = Q_best;
                        ibad = 0;
                    }
                }
                else if (Q > Q_old)
                { // The direction was ok

                    // Find parameter multipliers, to continue in the same direction that improved Q (with some added momentum and stochasticity)
                    if (update9)
                    {
                        float RT0 = myRNG->Gaus(1.1, 0.1) * snn_in.Threshold[0] /  snn_old.Threshold[0];
                        float RT1 = myRNG->Gaus(1.1, 0.1) * snn_in.Threshold[1] /  snn_old.Threshold[1];
                        float Ra = myRNG->Gaus(1.1, 0.1) * snn_in.alpha /  snn_old.alpha;
                        float RI = myRNG->Gaus(1.1, 0.1) * snn_in.L1inhibitfactor /  snn_old.L1inhibitfactor;
                        float RK = myRNG->Gaus(1.1, 0.1) * snn_in.K /  snn_old.K;
                        float RK1 = myRNG->Gaus(1.1, 0.1) * snn_in.K1 /  snn_old.K1;
                        float RK2 = myRNG->Gaus(1.1, 0.1) * snn_in.K2 /  snn_old.K2;
                        float RIE = myRNG->Gaus(1.1, 0.1) * snn_in.IE_Pot_const /  snn_old.IE_Pot_const;
                        float RID = myRNG->Gaus(1.1, 0.1) * snn_in.IPSP_dt_dilation /  snn_old.IPSP_dt_dilation;
                        // Store previous values
                        snn_old.Threshold[0]= snn_in.Threshold[0];
                        snn_old.Threshold[1] = snn_in.Threshold[1];
                        snn_old.alpha = snn_in.alpha;
                        snn_old.L1inhibitfactor = snn_in.L1inhibitfactor;
                        snn_old.K = snn_in.K;
                        snn_old.K1 = snn_in.K1;
                        snn_old.K2 = snn_in.K2;
                        snn_old.IE_Pot_const = snn_in.IE_Pot_const;
                        snn_old.IPSP_dt_dilation = snn_in.IPSP_dt_dilation;
                        // And update parameters
                        Optvar[0] = RT0 - 1.;
                        Optvar[1] = RT1 - 1.;
                        Optvar[2] = Ra - 1.;
                        Optvar[3] = RI - 1.;
                        Optvar[4] = RK - 1.;
                        Optvar[5] = RK1 - 1.;
                        Optvar[6] = RK2 - 1.;
                        Optvar[7] = RIE - 1.;
                        Optvar[8] = RID - 1.;
                        snn_in.Threshold[0] *= RT0;
                        snn_in.Threshold[1] *= RT1;
                        snn_in.alpha *= Ra;
                        snn_in.L1inhibitfactor *= RI;
                        snn_in.K *= RK;
                        snn_in.K1 *= RK1;
                        snn_in.K2 *= RK2;
                        snn_in.IE_Pot_const *= RIE;
                        snn_in.IPSP_dt_dilation *= RID;
                    }
                    // Same story, for delays
                    if (updateDelays)
                    {
                        for (int in = 0; in < snn_in.N_neurons; in++)
                        {
                            for (int is = 0; is < snn_in.N_streams; is++)
                            {
                                float R = myRNG->Gaus(1.1, 0.1) * (snn_in.Delay[in][is] - snn_old.Delay[in][is]);
                                snn_old.Delay[in][is] = snn_in.Delay[in][is];
                                if (R > 0.)
                                {
                                    snn_in.Delay[in][is] += R;
                                    OptvarD[in * snn_in.N_streams + is] = R;
                                }
                                if (snn_in.Delay[in][is] > snn_in.MaxDelay)
                                    snn_in.Delay[in][is] = snn_in.MaxDelay;
                                if (snn_in.Delay[in][is] < 0.)
                                    snn_in.Delay[in][is] = 0.;
                            }
                        }
                    }
                    if (updateConnections)
                    {
                        // Need to do nothing, no "direction" to go further in
                    }
                    Q_old = Q;
                }

            } // if not iepoch==1

            if (ievent < N_events - 1)
            { // Otherwise we graciously exit loop
                if (update9)
                {
                    cout << " - Try TL0 = " << snn_in.Threshold[0]
                         << " TL1 = " << snn_in.Threshold[1] << " a = " << snn_in.alpha << " L1inh = " << snn_in.L1inhibitfactor
                         << " K = " << snn_in.K << " K1 = " << snn_in.K1 << " K2 = " << snn_in.K2 << " IEPC = " << snn_in.IE_Pot_const << " IPSPdf = " << snn_in.IPSP_dt_dilation << endl
                         << endl;
                }
                else
                {
                    cout << endl
                         << endl;
                }

                // Reset a few counters
                for (int in = 0; in < snn_in.N_neurons; in++)
                {
                    N_fires[in] = 0.;
                }
                for (int in = 0; in < snn_in.N_neurons; in++)
                {
                    random_fire[in] = 0;
                    for (int ic = 0; ic < N_classes; ic++)
                    {
                        fired_sum[ic][in] = 0;
                    }
                }
                for (int ic = 0; ic < N_classes; ic++)
                {
                    gen_sum[ic] = 0;
                    fired_anyL0[ic] = 0;
                    fired_anyL1[ic] = 0;
                }
                atleastonefired = 0;

                // Reset progress bar
                if (doprogress)
                {
                    cout << "         " << progress[0];
                    currchar = 1;
                }
            }

        } // if ievent+1%NevPerEpoch = 0
        ievent++; // only go to next event if we did a backward pass too
    } while (ievent < N_events);

    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        snn_in.sumweight[in]=0;
        for (int is = 0; is < snn_in.N_streams; is++)
        {
            if(!snn_in.Void_weight[in][is]) snn_in.sumweight[in]+=snn_in.Weight[in][is];
        }

        cout << in << " " << snn_in.sumweight[in] << endl;
    }

    // closing the input file
    delete IT;
    delete OT;
    delete dirIT;
    delete dirOT;

    fout.close();
    file->Close();
    delete file;

    // Draw histograms

    cout << "Drawing histos" << endl;
    TCanvas *S = new TCanvas("S", "", 3000, 600);
    S->Divide(5, 2);
    for (int i = 0; i < 10; i++)
    {
        // if i is even, drow B and S
        if (i < 5)
        {
            S->cd(i + 1);
            StreamsB[i]->SetLineColor(kRed);
            StreamsB[i]->Draw("BOX");
            StreamsS[i]->SetLineColor(kBlue);
            StreamsS[i]->Draw("BOXSAME");
        }
        // if i is odd, draw B and N
        else
        {
            S->cd(i + 1);
            StreamsN[i - 5]->SetLineColor(kGreen);
            StreamsN[i - 5]->Draw("BOXSAME");
        }
    }

    TCanvas *C = new TCanvas("C", "", 1000, 1000);
    C->Divide(N_classes, snn_in.N_neurons);
    for (int i = 0; i < snn_in.N_neurons * N_classes; i++)
    {
        C->cd(i + 1);
        Latency[i]->Draw("COL4");
    }

    TCanvas *E0 = new TCanvas("E0", "", 800, 800);
    E0->Divide(N_classes, snn_in.N_neuronsL[0]);
    for (int i = 0; i < snn_in.N_neuronsL[0] * N_classes; i++)
    {
        E0->cd(i + 1);
        Efficiency[i]->SetMaximum(1.1);
        Efficiency[i]->SetMinimum(0.);
        Efficiency[i]->Draw("");
        int in = i / N_classes;
        FakeRate[in]->SetMarkerColor(2);
        FakeRate[in]->SetLineColor(2);
        FakeRate[in]->Draw("SAME");
        int ic = i % N_classes;
        Eff_totL0[ic]->SetMarkerColor(3);
        Eff_totL0[ic]->SetLineColor(3);
        Eff_totL0[ic]->Draw("SAME");
        Efficiency[i]->Draw("SAME");
    }

    TCanvas *E1 = new TCanvas("E1", "", 800, 800);
    E1->Divide(N_classes, snn_in.N_neuronsL[1]);
    for (int i = snn_in.N_neuronsL[0] * N_classes; i < snn_in.N_neurons * N_classes; i++)
    {
        E1->cd(i + 1 - snn_in.N_neuronsL[0] * N_classes);
        Efficiency[i]->SetMaximum(1.1);
        Efficiency[i]->SetMinimum(0.);
        Efficiency[i]->Draw("");
        int in = i / N_classes;
        FakeRate[in]->SetMarkerColor(2);
        FakeRate[in]->SetLineColor(2);
        FakeRate[in]->Draw("SAME");
        int ic = i % N_classes;
        Eff_totL1[ic]->SetMarkerColor(3);
        Eff_totL1[ic]->SetLineColor(3);
        Eff_totL1[ic]->Draw("SAME");
        Efficiency[i]->Draw("SAME");
    }

    // Plot the efficiencies and acceptances for the best q-value run
    // --------------------------------------------------------------
    for (int i = 0; i < snn_in.N_neurons * N_classes; i++)
    {
        int in = i / N_classes;
        int ic = i % N_classes;
        BestEff[in]->SetBinContent(ic + 1, Efficiency[i]->GetBinContent(ind_qbest));
        BestFR[in]->SetBinContent(ic + 1, FakeRate[in]->GetBinContent(ind_qbest));
        if (in < snn_in.N_neuronsL[0])
        {
            BestEtot[in]->SetBinContent(ic + 1, Eff_totL0[ic]->GetBinContent(ind_qbest));
        }
        else
        {
            BestEtot[in]->SetBinContent(ic + 1, Eff_totL1[ic]->GetBinContent(ind_qbest));
        }
    }
    TCanvas *BE = new TCanvas("BE", "", 400, 800);
    BE->Divide(1, snn_in.N_neurons);
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        BE->cd(in + 1);
        BestEff[in]->SetMaximum(1.1);
        BestEff[in]->SetMinimum(0.);
        BestEff[in]->Draw("");
        BestFR[in]->SetMarkerColor(2);
        BestFR[in]->SetLineColor(2);
        BestFR[in]->Draw("SAME");
        BestEtot[in]->SetMarkerColor(3);
        BestEtot[in]->SetLineColor(3);
        BestEtot[in]->Draw("SAME");
        BestEff[in]->Draw("SAME");
    }
    TCanvas *MW = new TCanvas("MW", "", 400, 800);
    MW->Divide(1, snn_in.N_neurons);
    MW->SetLogy();
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        MW->cd(in + 1);
        // HRMSWeight[in]->SetMaximum(1.1);
        // HRMSWeight[in]->SetMinimum(0);
        HRMSWeight[in]->SetMarkerColor(kGreen);
        HRMSWeight[in]->SetLineColor(kGreen);
        HRMSWeight[in]->Draw("");
        HMaxWeight[in]->SetMarkerColor(kRed);
        HMaxWeight[in]->SetLineColor(kRed);
        HMaxWeight[in]->Draw("SAME");
        HMinWeight[in]->SetMarkerColor(kBlue);
        HMinWeight[in]->SetLineColor(kBlue);
        HMinWeight[in]->Draw("SAME");
        HRMSWeight[in]->Draw("SAME");
    }

    TCanvas *SE = new TCanvas("SE", "", 800, 400);
    SE->Divide(3, 1);
    SE->cd(1);
    SelectivityL0->Draw();
    SE->cd(2);
    SelectivityL1->Draw();
    SE->cd(3);
    Qvalue->Draw();
    Qmax->Draw("SAME");

    TCanvas *W = new TCanvas("W", "", 500, 800);
    int nrow = snn_in.N_neurons / 2;
    if (snn_in.N_neurons % 2 != 0)
        nrow += 1;
    W->Divide(2, nrow);
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        W->cd(in + 1);
        for (int is = 0; is < snn_in.N_streams; is++)
        {
            // HWeight[in * N_streams + is]->SetMaximum(1.1);
            HWeight[in * snn_in.N_streams + is]->SetMinimum(0.);
            int color = is + 1;
            if (color == 8)
                color = 9;
            HWeight[in * snn_in.N_streams + is]->SetLineColor(color);
            if (is == 0)
            {
                HWeight[in * snn_in.N_streams + is]->Draw();
            }
            else
            {
                HWeight[in * snn_in.N_streams + is]->Draw("SAME");
            }
        }
    }

    // Draw final Efficiency and acceptance maps
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        for (int ic = 0; ic < N_classes; ic++)
        {
            EffMap->SetBinContent(in + 1, ic + 1, Efficiency[ic + in * N_classes]->GetBinContent(ind_qbest));
        }
    }
    TCanvas *Y = new TCanvas("Y", "", 600, 900);
    Y->cd();
    EffMap->Draw("COL4");

    // Final Statistics
    // ----------------
    cout << endl
         << endl;
    cout << "         Run parameters" << endl;
    cout << "         -----------------------------------" << endl;
    cout << "                       L0 neurons: " << snn_in.N_neuronsL[0] << endl;
    cout << "                       L1 neurons: " << snn_in.N_neuronsL[1] << endl;
    cout << "            Connected L0-L1 frac.: " << snn_in.CF01 << endl;
    cout << "            Connected IN-L0 frac.: " << snn_in.CFI0 << endl;
    cout << "            Connected IN-L1 frac.: " << snn_in.CFI1 << endl;
    cout << "                    Track classes: " << N_classes << endl;
    cout << endl;
    cout << "         Optimization results" << endl;
    cout << "         -----------------------------------" << endl;
    cout << "               Average efficiency: " << Eff_best << endl;
    cout << "                Average fake rate: " << Acc_best << endl;
    cout << "                  Maximum Q value: " << Q_best << endl;
    cout << "                   L1 selectivity: " << SelL1_best << endl;
    cout << endl;
    cout << "         Optimized parameter values" << endl;
    cout << "         -----------------------------------" << endl;
    cout << "                     L0 threshold: " << snn_best.Threshold[0] << endl;
    cout << "                     L1 threshold: " << snn_best.Threshold[1] << endl;
    cout << "                            alpha: " << snn_best.alpha << endl;
    cout << "                        L1inhibit: " << snn_best.L1inhibitfactor << endl;
    cout << "                                K: " << snn_best.K << endl;
    cout << "                               K1: " << snn_best.K1 << endl;
    cout << "                               K2: " << snn_best.K2 << endl;
    cout << "                 IE pot. constant: " << snn_best.IE_Pot_const << endl;
    cout << "                 IPSP dt dilation: " << snn_best.IPSP_dt_dilation << endl;
    cout << "         -----------------------------------" << endl;
    cout << endl;

    // Dump to file optimized parameters and results
    // ---------------------------------------------
    /*
    TODO
     if (N_epochs > 1)
    {
        Write_Parameters(); // This also defines indfile, used below
    }
    */
   

    // Dump histograms to root file
    string Path = "./MODE/SNNT/";
    std::stringstream sstr;
    char num[80];
    sprintf(num, "NL0=%d_NL1=%d_NCl=%d_CF01=%.2f_CFI0=%.2f_CFI1=%.2f_alfa=%.2f_%d", snn_in.N_neuronsL[0], snn_in.N_neuronsL[1], N_classes, snn_in.CF01, snn_in.CFI0, snn_in.CFI1, snn_in.alpha, indfile);
    sstr << "Histos13_";
    string namerootfile = Path + sstr.str() + num + ".root";
    TFile *rootfile = new TFile(namerootfile.c_str(), "RECREATE");
    rootfile->cd();

    // Write canvases first
    cout << "Writing canvas" << endl;
    S->Write();
    C->Write();
    E0->Write();
    E1->Write();
    BE->Write();
    SE->Write();
    W->Write();
    Y->Write();
    MW->Write();

    cout << "Saving pdf" << endl;
    S->SaveAs("./pdf/S.pdf", "pdf");
    C->SaveAs("./pdf/C.pdf", "pdf");
    E0->SaveAs("./pdf/E0.pdf", "pdf");
    E1->SaveAs("./pdf/E1.pdf", "pdf");
    BE->SaveAs("./pdf/BE.pdf", "pdf");
    SE->SaveAs("./pdf/SE.pdf", "pdf");
    W->SaveAs("./pdf/W.pdf", "pdf");
    Y->SaveAs("./pdf/Y.pdf", "pdf");
    MW->SaveAs("./pdf/MWlog.pdf", "pdf");

    // Then histograms
    SelectivityL0->Write();
    SelectivityL1->Write();
    Qvalue->Write();
    Qmax->Write();
    HEff->Write();
    HAcc->Write();
    HT0->Write();
    HT1->Write();
    HA->Write();
    HL1IF->Write();
    HK->Write();
    HK1->Write();
    HK2->Write();
    HIEPC->Write();
    HVoidWs->Write();
    HDelays->Write();
    Q_12->Write();
    Q_34->Write();
    Q_56->Write();
    Q_78->Write();
    Q_MV->Write();
    N_12->Write();
    N_34->Write();
    N_56->Write();
    N_78->Write();
    N_MV->Write();
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        for (int ic = 0; ic < N_classes; ic++)
        {
            int id = in * N_classes + ic;
            Latency[id]->Write();
            Efficiency[id]->Write();
        }
        for (int is = 0; is < snn_in.N_streams; is++)
        {
            int id = in * snn_in.N_streams + is;
            HWeight[id]->Write();
        }
        FakeRate[in]->Write();
        BestEff[in]->Write();
        BestFR[in]->Write();
        BestEtot[in]->Write();
        // HRMSWeight[in]->Write();
        // HMaxWeight[in]->Write();
        // HMinWeight[in]->Write();
    }
    for (int ic = 0; ic < N_classes; ic++)
    {
        Eff_totL0[ic]->Write();
        Eff_totL1[ic]->Write();
    }
    for (int i = 0; i < 10; i++)
    {
        StreamsS[i]->Write();
        StreamsB[i]->Write();
        StreamsN[i]->Write();
    }
    EffMap->Write();

    rootfile->Write();

    // End of program
    rootfile->Close();

    TFile *rootfile2 = new TFile("output.root", "RECREATE");
    MW->Write();
    for (int in = 0; in < snn_in.N_neurons; in++)
    {

        HRMSWeight[in]->Write();
        HMaxWeight[in]->Write();
        HMinWeight[in]->Write();
    }
    rootfile2->Close();
    gROOT->Time();


    return;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// Main routine
// ------------
int main(int argc, char *argv[])
{
    // Loop through the command-line arguments
    for (int i = 1; i < argc; ++i)
    {
        const char* arg = argv[i];
        cout << arg << endl;
        if (strcmp(arg,"--NL0")==0)
            _NL0 = stoi(argv[i + 1]);
        else if (arg == "--NL1")
            _NL1 = stoi(argv[i + 1]);


        else if (arg == "--alpha")
            _alpha = stof(argv[i + 1]);


        else if (arg == "--CF01")
            _CF01 = stof(argv[i + 1]);
        else if (arg == "--CFI0")
            _CFI0 = stof(argv[i + 1]);
        else if (arg == "--CFI1")
            _CFI1 = stof(argv[i + 1]);


        else if (arg == "--L1inhibitfactor")
            _L1inhibitfactor = stof(argv[i + 1]);


        else if (arg == "--K")
            _K = stof(argv[i + 1]);
        else if (arg == "--K1")
            _K1 = stof(argv[i + 1]);
        else if (arg == "--K2")
            _K2 = stof(argv[i + 1]);


        else if (arg == "--IE_Pot_const")
            _IE_Pot_const = stof(argv[i + 1]);
        else if (arg == "--IPSP_dt_dilation")
            _IPSP_dt_dilation = stof(argv[i + 1]);
        

        else if (arg == "--MaxDelay")
            _MaxDelay = stof(argv[i + 1]);


        else if (arg == "--tau_m")
            _tau_m = stof(argv[i + 1]);
        else if (arg == "--tau_s")
            _tau_s = stof(argv[i + 1]);
        else if (arg == "--tau_r")
            _tau_r = stof(argv[i + 1]);
        else if (arg == "--tau_plus")
            _tau_plus = stof(argv[i + 1]);
        else if (arg == "--tau_minus")
            _tau_minus = stof(argv[i + 1]);


        else if (arg == "--a_plus")
            _a_plus = stof(argv[i + 1]);
        else if (arg == "--a_minus")
            _a_minus = stof(argv[i + 1]);


        else if (arg == "--N_InputStreams")
            _N_InputStreams = stoi(argv[i + 1]);



        else if (arg == "--Threshold")
            _Threshold0 = stof(argv[i + 1]);
        else if (arg == "--ThresholdL1")
            _Threshold1 = stof(argv[i + 1]);



        // Main's parameters 
        else if (arg == "--N_ev")
            N_events = stoi(argv[i + 1]);
        else if (arg == "--N_ep")
            N_epochs = stoi(argv[i + 1]);
        else if (arg == "--batch")
            batch = stoi(argv[i + 1]);
        else if (arg == "--rootInput"){
            cout << argv[i+1] <<endl;
            rootInput = argv[i + 1];              
        }
        else if(arg == "--N_classes")
            N_classes = stoi(argv[i + 1]);
        else if(arg == "--TrainingCode")
            TrainingCode = stoi(argv[i + 1]);
        else if(arg == "--ReadPars")
            ReadPars = stoi(argv[i + 1]);
        else if(arg == "--NROOT")
            NROOT = stoi(argv[i + 1]);
        else
        {
            cerr << "Invalid flag" << endl;
            exit(1);
        }
        
        
    }
    SNN S(_NL0,  _NL1,
    _alpha,
          _CFI0, _CFI1, _CF01,
          _L1inhibitfactor,
          _K, _K1, _K2,
          _IE_Pot_const, _IPSP_dt_dilation,
          _MaxDelay,

          _tau_m,  _tau_s,  _tau_r,  _tau_plus,  _tau_minus,
          _a_plus,  _a_minus,

          _N_InputStreams,
          _Threshold0,  _Threshold1);
    
    SNN_Tracking(S);
    return 0;
}