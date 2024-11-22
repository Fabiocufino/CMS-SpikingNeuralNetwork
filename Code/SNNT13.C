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
#include <tuple>
#include <stdexcept>

#include "Class/SNN.h"

using json = nlohmann::json;
using namespace std;

static bool insert = true;

// Constants and data used throughout the code
// -------------------------------------------
static int N_part; // Number of generated particles in an event
static double First_angle;
static float *Eff; // Efficiency of each neuron to signals of different classes
static float *Eff_window; // Efficiency of each neuron to signals of different classes
static vector<double> PreSpike_Time;
static vector<int> PreSpike_Stream;
static vector<int> PreSpike_Signal; // 0 for background hit, 1 for signal hit, 2 for L1 neuron spike
static vector<int> PreSpike_Class;
static vector<int> neurons_index;   // contains neurons identifiers in random positions
static float Q_best_L0;
static float SelL0_best;
static float Eff_best_L0;
static float Acc_best_L0;
static float Q_best_L1;
static float SelL1_best;
static float SelTOT_best;
static float Eff_best_L1;
static float Acc_best_L1;
static int indfile;
static char progress[53] = "[--10%--20%--30%--40%--50%--60%--70%--80%--90%-100%]"; // Progress bar
static int ievent;
static float max_fake = 0;

// New random number generator
static TRandom3 *myRNG = new TRandom3(static_cast<unsigned int>(time(0)));

// clear hits vector
void Reset_hits()
{
    hit_pos.clear();
    return;
}

// Start binning r-z plane ---------------
int GetBinR(double r_hit)
{
    if (r_hit < 0 || r_hit > max_R)
        return N_bin_r - 1;

    // 10 bins are associated to a tracking layer
    for (int i = 0; i < N_TrackingLayers; i++)
    {
        if (r_hit > Left_Layers[i] && r_hit < Right_Layers[i])
            return i;
        else if (r_hit < Left_Layers[i+1])
            return i;
    }
    return N_bin_r - 1;
}

int GetBinZ(double z)
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
        PreSpike_Signal.push_back(row.id); //1,2 -> respectively Backgroung, Signal
        PreSpike_Class.push_back(row.pclass);
    }

    // rescan from [0, delta]
    for (auto &&row : hit_pos)
    {
        if (row.phi > delta)
            break;
        double time = t_in + (row.phi + M_PI * 2.) / omega;

        int itl = GetStreamID(GetBinR(row.r), GetBinZ(row.z));

        PreSpike_Time.push_back(time);
        PreSpike_Stream.push_back(itl);
        PreSpike_Signal.push_back(row.id); // 1,2 -> respectively Backgroung, Signal
        if(row.pclass >= 0) PreSpike_Class.push_back(row.pclass + N_classes); //ghost particle -> I indicate it with a shifted pclass
        else PreSpike_Class.push_back(row.pclass);
    }
}

// Calculate selectivity of set of neurons
// ---------------------------------------
float Compute_Selectivity(int level, int mode, SNN &snn)
{
    float S = 0.;
    int inmin, inmax;
    if (level == 0)
    {
        inmin = 0;
        inmax = snn.N_neuronsL[0];
    }
    else if (level == 1)
    {
        inmin = snn.N_neuronsL[0];
        inmax = snn.N_neurons;
    }
    else if (level == -1){
        inmin = 0;
        inmax = snn.N_neurons;
    }

    if (mode == 2)
    {   // compute mutual information
        // I(N_neurons,N_classes) = Sum_i^N_n Sum_j^N_c Eff(i,j) log_2 [Eff(i,j)/Eff_i Eff_j)]
        // where  is the average efficiency of neuron i over classes, and Eff_j is the
        // average efficiency on class j over neurons
        S = 0.;
        float Effn[snn.N_neurons];
        float Effc[N_ev_classes];
        float sumeff = 0.;
        for (int in = inmin; in < inmax; in++)
        {
            sumeff = 0.;
            for (int ic = 0; ic < N_ev_classes; ic++)
            {
                sumeff += Eff[ic + N_ev_classes * in];
            }
            if (N_ev_classes > 0)
                Effn[in] = sumeff / N_ev_classes;
        }
        for (int ic = 0; ic < N_ev_classes; ic++)
        {
            sumeff = 0.;
            for (int in = inmin; in < inmax; in++)
            {
                sumeff += Eff[ic + N_ev_classes * in];
            }
            if (inmax > inmin)
                Effc[ic] = sumeff / (inmax - inmin);
        }
        for (int in = inmin; in < inmax; in++)
        {
            for (int ic = 0; ic < N_ev_classes; ic++)
            {
                if (Effc[ic] * Effn[in] > 0.)
                    S += Eff[ic + N_ev_classes * in] *
                         (log2(Eff[ic + N_ev_classes * in] + epsilon) - log2(Effc[ic] * Effn[in]));
            }
        }
    }
    if (inmax > inmin)
        S /= N_ev_classes * (inmax - inmin);
    return S;
}

float computeMutualInformation(
    int level,
    int *gen_sum,
    int **fired_sum,
    SNN &snn,
    bool equiprobable_classes = true)
{
    int inmin, inmax;
    if (level == 0)
    {
        inmin = 0;
        inmax = snn.N_neuronsL[0];
    }
    else if (level == 1)
    {
        inmin = snn.N_neuronsL[0];
        inmax = snn.N_neurons;
    }
    else if (level == -1){
        inmin = 0;
        inmax = snn.N_neurons;
    }
    const float epsilon = 1e-10f;
    float N_total = 0.0f;
    
    // Compute total number of events
    for (int ic = 0; ic < N_ev_classes; ic++) {
        N_total += gen_sum[ic];
    }

    // Compute joint probabilities P(neuron fires, class ic)
    vector<vector<float>> P_joint(snn.N_neurons, vector<float>(N_ev_classes, 0.0f));
    for (int in = inmin; in < inmax; in++) {
        for (int ic = 0; ic < N_ev_classes; ic++) {
            P_joint[in][ic] = fired_sum[ic][in] / N_total;
        }
    }

    // Compute marginal probabilities for neurons P(neuron fires)
    vector<float> P_neuron(snn.N_neurons, 0.0f);
    for (int in = inmin; in < inmax; in++) {
        for (int ic = 0; ic < N_ev_classes; ic++) {
            P_neuron[in] += P_joint[in][ic];
        }
    }

    // Compute marginal probabilities for classes P(class ic)
    vector<float> P_class(N_ev_classes, 0.0f);
    if (equiprobable_classes) {
        for (int ic = 0; ic < N_ev_classes; ic++) {
            P_class[ic] = 1.0f / N_ev_classes;
        }
    } else {
        for (int ic = 0; ic < N_ev_classes; ic++) {
            P_class[ic] = gen_sum[ic] / N_total;
        }
    }

    // Compute mutual information
    float I = 0.0f;
    for (int in = inmin; in < inmax; in++) {
        for (int ic = 0; ic < N_ev_classes; ic++) {
            float P_joint_value = P_joint[in][ic];
            float P_neuron_value = P_neuron[in];
            float P_class_value = P_class[ic];

            if (P_joint_value > 0 && P_neuron_value > 0 && P_class_value > 0) {
                I += P_joint_value * (log2(P_joint_value + epsilon) - log2(P_neuron_value + epsilon) - log2(P_class_value + epsilon));
            }
        }
    }

    return I;
}

vector<float> computeMutualInformationPerNeuron(
    int N_ev_classes,
    int N_neurons,
    const vector<float>& gen_sum,
    const vector<vector<float>>& fired_sum,
    bool equiprobable_classes = true)
{
    const float epsilon = 1e-10f;
    float N_total = 0.0f;

    // Compute total number of events
    for (int ic = 0; ic < N_ev_classes; ic++) {
        N_total += gen_sum[ic];
    }

    // Compute joint probabilities P(neuron fires, class ic)
    vector<vector<float>> P_joint(N_neurons, vector<float>(N_ev_classes, 0.0f));
    for (int in = 0; in < N_neurons; in++) {
        for (int ic = 0; ic < N_ev_classes; ic++) {
            P_joint[in][ic] = fired_sum[ic][in] / N_total;
        }
    }

    // Compute marginal probabilities for neurons P(neuron fires)
    vector<float> P_neuron(N_neurons, 0.0f);
    for (int in = 0; in < N_neurons; in++) {
        for (int ic = 0; ic < N_ev_classes; ic++) {
            P_neuron[in] += P_joint[in][ic];
        }
    }

    // Compute marginal probabilities for classes P(class ic)
    vector<float> P_class(N_ev_classes, 0.0f);
    if (equiprobable_classes) {
        for (int ic = 0; ic < N_ev_classes; ic++) {
            P_class[ic] = 1.0f / N_ev_classes;
        }
    } else {
        for (int ic = 0; ic < N_ev_classes; ic++) {
            P_class[ic] = gen_sum[ic] / N_total;
        }
    }

    // Compute mutual information per neuron
    vector<float> MI_per_neuron(N_neurons, 0.0f);
    for (int in = 0; in < N_neurons; in++) {
        float I = 0.0f;
        for (int ic = 0; ic < N_ev_classes; ic++) {
            float P_joint_value = P_joint[in][ic];
            float P_neuron_value = P_neuron[in];
            float P_class_value = P_class[ic];

            if (P_joint_value > 0 && P_neuron_value > 0 && P_class_value > 0) {
                I += P_joint_value * (log2(P_joint_value + epsilon) - log2(P_neuron_value + epsilon) - log2(P_class_value + epsilon));
            }
        }
        MI_per_neuron[in] = I;
    }

    return MI_per_neuron;
}

// Compute Q-value
// ---------------
float Compute_Q(float eff, float acc, float sel)
{
    float Q0 = eff / sqrt(acc);
    float w = 5. * (exp(Q0 / 4.) - 1.) / (exp(1.) - 1.);
    return Q0 + w * sel;
}

tuple<TFile*, TTree*, TTree*, TTree*> readRootFile(const string &rootInput)
{
    // Open the ROOT file in read mode
    TFile *file = TFile::Open(rootInput.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        throw runtime_error("Error: Cannot open file " + rootInput);
    }

    // Access the required directories
    TDirectoryFile *dirIT = dynamic_cast<TDirectoryFile *>(file->Get("clusterValidIT"));
    TDirectoryFile *dirOT = dynamic_cast<TDirectoryFile *>(file->Get("clusterValidOT"));
    TDirectoryFile *dirEV = dynamic_cast<TDirectoryFile *>(file->Get("classification"));

    if (!dirIT)
    {
        file->Close();
        throw runtime_error("Error: Cannot access directory clusterValidIT in file " + rootInput);
    }

    if (!dirOT)
    {
        file->Close();
        throw runtime_error("Error: Cannot access directory clusterValidOT in file " + rootInput);
    }

    if (!dirEV)
    {
        file->Close();
        throw runtime_error("Error: Cannot access directory classification in file " + rootInput);
    }

    // Extract the trees from the directories
    TTree *IT = dynamic_cast<TTree *>(dirIT->Get("tree"));
    TTree *OT = dynamic_cast<TTree *>(dirOT->Get("tree"));
    TTree *ET = dynamic_cast<TTree *>(dirEV->Get("event_tree"));

    if (!IT)
    {
        file->Close();
        throw runtime_error("Error: Cannot access tree in clusterValidIT in file " + rootInput);
    }

    if (!OT)
    {
        file->Close();
        throw runtime_error("Error: Cannot access tree in clusterValidOT in file " + rootInput);
    }

    if (!ET)
    {
        file->Close();
        throw runtime_error("Error: Cannot access event_tree in classification in file " + rootInput);
    }

    // Do not close the file here; it will be managed by the caller
    return {file, IT, OT, ET};
}

// To read our preprocessed file
void ReadFromProcessed(TTree *IT, TTree *OT, TTree *ET, int id_event_value)
{
    Reset_hits();
    First_angle = max_angle;

    //retrieve the information about the eventClass
    int eventClass;
    ET->SetBranchAddress("eventClass", &eventClass);
    ET->GetEntry(id_event_value-1);
    pclass = eventClass;
    
    //TODO: generalize to more than 2 particles with an integer division
    if(eventClass < 0) {N_part=0; pclass=0;}  
    else if(eventClass >= 0 && eventClass < N_classes) N_part=1;
    else N_part = 2;
    
    // Reading the Inner Tracker

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

    if (ievent % NROOT == 0)
    {
        last_row_event_IT = 0;
        last_row_event_OT = 0;
    }

    // Loop over entries and find rows with the specified id_event value
    for (int i = last_row_event_IT; i < IT->GetEntries(); ++i)
    {
        IT->GetEntry(i);

        if (static_cast<int>(id_event) != id_event_value)
        {
            last_row_event_IT = i;
            break;
        }
        phi += M_PI;
        if (static_cast<int>(type) == 1)
        {
            type = SIG;
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

        hit_pos.emplace_back(r, z, phi, static_cast<int>(type), cluster_pclass);
    }

    // OUT Tracker

    OT->SetBranchAddress("cluster_z", &z);
    OT->SetBranchAddress("cluster_R", &r);
    OT->SetBranchAddress("cluster_phi", &phi);
    OT->SetBranchAddress("eventID", &id_event);
    OT->SetBranchAddress("cluster_type", &type);
    OT->SetBranchAddress("pclass", &cluster_pclass);

    for (int i = last_row_event_OT; i < OT->GetEntries(); ++i)
    {
        OT->GetEntry(i);
        if (static_cast<int>(id_event) != id_event_value)
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

        hit_pos.emplace_back(r, z, phi, static_cast<int>(type), cluster_pclass);
    }
}

// Plot potential stuff
void ReadWeights(TFile *file, SNN &P)
{
    vector<TH1F *> Hvec;
    int iw = 0;
    const char *name = "HWeight";
    while (true)
    {
        TH1F *hist = nullptr;
        char buffer[50];
        sprintf(buffer, "%s%d", name, iw);

        hist = dynamic_cast<TH1F *>(file->Get(buffer));
        if (hist == nullptr)
            break;
        Hvec.push_back(hist);
        iw++;
    }
    cout << "Loaded " << Hvec.size() << " histograms" << endl;
    cout << "Extracting last weights configuration " << endl;
    for (int in = 0; in < P.N_neurons; in++)
    {
        for (int is = 0; is < P.N_streams; is++)
        {
            // Get the number of bins in the x-axis
            int lastBin = Hvec[in * P.N_streams + is]->GetNbinsX();
            // Get the content of the last bin
            float lastBinValue = Hvec[in * P.N_streams + is]->GetBinContent(lastBin);
            // weight = -1 -> inexisting connection
            P.Void_weight[in][is] = (lastBinValue == -1);
            P.Weight[in][is] = lastBinValue;
        }
    }
    cout << "Weights loaded successfully" << endl;
}


void appendToJson(const string& filename, float Efficiency, float Fake_rate, float Selectivity, float Q_value) {
    ifstream file_in(filename);
    if (!file_in.is_open()) {
        cerr << "Unable to open file: " << filename << endl;
        return;
    }

    json j;
    file_in >> j; // Load the JSON content into j
    file_in.close();

    // Append new values
    j["Efficiency"] = Efficiency;
    j["Fake_rate"] = Fake_rate;
    j["Q"] = Q_value;
    j["Selectivity"] = Selectivity;

    // Open the file for writing (overwrite with updated JSON)
    ofstream file_out(filename);
    if (!file_out.is_open()) {
        cerr << "Unable to open file for writing: " << filename << endl;
        return;
    }

    file_out << j.dump(4); 
    file_out.close();
}

// plot neuron potentials as a function of time
//TODO adapt it to multiple particles
void PlotPotentials(string rootInput, SNN &P, int _N_events)
{
    insert = true;
    vector<int> neurons_index;
    // initialization of neurons_index vector
    for (int i = 0; i < P.N_neurons; i++)
        neurons_index.push_back(i);
    N_events = _N_events;

    // vectors to plot
    vector<double> Time[N_events];
    vector<float> Potential[N_events][P.N_neurons];
    int fire_count[P.N_neurons];
    for (int ic = 0; ic < P.N_neurons; ic++)
        fire_count[ic] = 0;

    cout << "Initializaing the plot SNN" << endl;

    auto [file, IT, OT, ET] = readRootFile(rootInput);
    
    int ievent = 0;
    double previous_firetime = 0;
    P.Init_neurons(-1);
    // Loop on events ----------------------------------------------
    do
    {
        PreSpike_Time.clear();
        PreSpike_Stream.clear();
        PreSpike_Signal.clear();
        PreSpike_Class.clear();

        if (ievent % NROOT == 0)
        {
            last_row_event_IT = 0;
            last_row_event_OT = 0;
        }
        ReadFromProcessed(IT, OT, ET, ievent % NROOT + 1);

        double t_in = ievent * (max_angle + Empty_buffer) / omega; // Initial time -> every event adds 25 ns
        Encode(t_in);

        // Loop on spikes and modify neuron and synapse potentials
        // -------------------------------------------------------
        for (int ispike = 0; ispike < PreSpike_Time.size(); ispike++)
        {
            // By looping to size(), we can insert along the way and still make it to the end
            double t = PreSpike_Time[ispike];

            // Modify neuron potentials based on synapse weights
            // -------------------------------------------------
            double min_fire_time = P.largenumber; // if no fire, neuron_firetime returns largenumber
            int in_first = -1;

            // Loop on neurons, but not in order to not favor any neuron
            // ---------------------------------------------------------

            // Shuffle order
            auto rng = default_random_engine{};
            shuffle(neurons_index.begin(), neurons_index.end(), rng);

            for (auto in : neurons_index)
            {

                // Compute future fire times of neurons and their order
                double fire_time = P.Neuron_firetime(in, t);

                if (fire_time < min_fire_time)
                {
                    in_first = in;
                    min_fire_time = fire_time;
                }
            }

            // no neuron is firing
            if (in_first == -1)
            {
                //TODO: check
                // finish the plot for the previous spike
                double t_prime;
                double delta_t;
                if (!insert)
                {
                    t_prime = previous_firetime;
                    delta_t = (t - t_prime) / 11;

                    for (int inc = 0; inc < 11; inc++)
                    {
                        Time[ievent].push_back(t_prime + inc * delta_t);
                        for (auto in : neurons_index)
                            Potential[ievent][in].push_back(P.Neuron_Potential(in, t_prime + inc * delta_t, false));
                    }
                }
                // plot the result of the incoming spike
                else
                {
                    if (ispike == 0)
                        t_prime = t_in;
                    else
                        t_prime = PreSpike_Time[ispike - 1];

                    delta_t = (t - t_prime) / 11;
                    for (int inc = 0; inc < 11; inc++)
                    {
                        Time[ievent].push_back(t_prime + inc * delta_t);
                        for (auto in : neurons_index)
                            Potential[ievent][in].push_back(P.Neuron_Potential(in, t_prime + inc * delta_t, false));
                    }
                }

                insert = true;
            }

            // Ok, neuron in_first is going to fire next.
            // Peek at next event in list, to see if it comes before in_first fires
            // --------------------------------------------------------------------
            else
            {
                double t_prime;
                double delta_t;
                // handle firing of neuron in_first
                // That means that we are handling the first neruon activation bewtween two EPSP spike

                if (insert)
                {
                    if (ispike == 0)
                        t_prime = min(t_in, min_fire_time);
                    else
                        t_prime = PreSpike_Time[ispike - 1];
                }
                else
                    t_prime = previous_firetime;
                delta_t = (min_fire_time - t_prime) / 11;

                for (int inc = 0; inc < 11; inc++)
                {
                    Time[ievent].push_back(t_prime + inc * delta_t);
                    for (auto in : neurons_index)
                        Potential[ievent][in].push_back(P.Neuron_Potential(in, t_prime + inc * delta_t, false));
                }

                //Handle the activation of the neuron in_first
                P.Activate_Neuron(in_first, min_fire_time);
                
                // IPSP for all others at relevant layer
                for (int in2 = 0; in2 < P.N_neurons; in2++)
                {
                    if (in2 != in_first)
                    {
                        if (P.Neuron_layer[in2] == P.Neuron_layer[in_first])
                        { // inhibitions within layer or across
                            P.insert_spike(in2, min_fire_time, P.IPSP, P.N_InputStreams + in_first, P.NOCLASS, ievent);
                        }
                    }
                }

                // Create EPS signal in L0 neuron-originated streams
                if (P.Neuron_layer[in_first] == 0)
                { // this is a Layer-0 neuron
                    for (int in = P.N_neuronsL[0]; in < P.N_neurons; in++)
                    {
                        int is = P.N_InputStreams + in_first;
                        if(!P.Void_weight[in][is]){
                            P.insert_spike(in, min_fire_time+ P.Delay[in][is], P.EPSP, is, P.NOCLASS, ievent);
                        }
                        
                    }
                }
                ispike -= 1;
                insert = false;

                // take a step back and search for another activation
                previous_firetime = min_fire_time;
            } // end if in_first fires

            // insert the new spike for the next iteration
            if (insert)
            {
                for (auto in : neurons_index)
                {
                    //  We implement a scheme where input streams produce an IE signal into L0, an EPS into L1, and L0 neurons EPS into L1
                    //  Add to neuron history, masking out L1 spikes for L0 neurons
                    int is = PreSpike_Stream[ispike];
                    if (!P.Void_weight[in][is])
                    {   // otherwise stream "is" does not lead to neuron "in"
                        // All input spikes lead to EPSP
                        P.insert_spike(in, t+ P.Delay[in][is], P.EPSP, is, PreSpike_Class[ispike], ievent);
                    }
                }
            }
        }
        ievent++; // only go to next event if we did a backward pass too
    } while (ievent < N_events);

    // dump the potentials inside a csv file
    ofstream outfile;
    outfile.open("MODE/potentials.csv");

    // Header
    outfile << "Event,Time";
    for (int in = 0; in < P.N_neurons; in++)
    {
        outfile << ",V(t)_" << in;
    }
    outfile << endl;

    // content
    for (ievent = 0; ievent < N_events; ievent++)
    {
        for (int it = 0; it < Time[ievent].size(); it++)
        {
            outfile << ievent << "," << Time[ievent][it];

            for (int in = 0; in < P.N_neurons; in++)
            {
                outfile << "," << Potential[ievent][in][it];
            }
            outfile << endl;
        }
    }

    // closing the input file
    delete IT;
    delete OT;

    file->Close();
    outfile.close();
    delete file;
}

vector<int> generate_unique_integers(int count, int range_start, int range_end) {
    if (count > (range_end - range_start + 1)) {
        throw invalid_argument("Count exceeds the range size.");
    }
    
    vector<int> numbers(range_end - range_start + 1);
    iota(numbers.begin(), numbers.end(), range_start);  // Fill with sequential numbers
    random_device rd;
    mt19937 g(rd());
    shuffle(numbers.begin(), numbers.end(), g);
    
    return vector<int>(numbers.begin(), numbers.begin() + count);  // Select first `count` elements
}

vector<int> generate_excluded_neurons(int exclude_L0, int exclude_L1, int N_neuronsL0, int N_neurons) {
    vector<int> excluded_L0 = generate_unique_integers(exclude_L0, 0, N_neuronsL0 - 1);
    vector<int> excluded_L1 = generate_unique_integers(exclude_L1, N_neuronsL0, N_neurons - 1);

    // Combine results if desired
    vector<int> result = excluded_L0;
    result.insert(result.end(), excluded_L1.begin(), excluded_L1.end());
    return result;
}

void SNN_Tracking(SNN &snn_in, int file_id_GS = -1)
{
    // The routine works as follows:
    // -----------------------------
    // Define tracker geometry ok
    // Initialize neuron potentials and synapse weights ok
    // Encode hits in spike streams
    // Loop on optimization cycles (N_epochs)
    //   Loop on time steps in event-based fashion and modify neuron and synapse potentials
    //     Determine when an output spike will occur and jump there
    //     If spike:
    //       - record spike
    //       - update synapses
    //       - update neuron
    //       - inhibit other neurons
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
    // -----------------------------------------------------------------------------------------------

    // Pass parameters can't update static values, so we need to reassign the latter
    if (batch)
        gROOT->SetBatch(kTRUE);

    // initialization of neurons_index vector
    for (int i = 0; i < snn_in.N_neurons; i++)
        neurons_index.push_back(i);

    // Initial checks
    // --------------
    if (N_events > MaxEvents)
    {
        cout << "  Sorry, max # of events is 10,000,000. Terminating." << endl;
        return;
    }
    NevPerEpoch = N_events / N_epochs;

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
    if (N_ev_classes > MaxClasses)
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
    cout << "                                    S   N   N      T r a c k i n g   " << endl;
    cout << endl;
    cout << "                                 ------------------------------------" << endl;
    cout << endl
         << endl
         << endl
         << endl;
    cout << "         ------------------------------------------------------------------------------------    " << endl;
    cout << "             Unsupervised search for tracks in Phase2 CMS Tracker with spiking neural network    " << endl;
    cout << "                                                                                      12/2024    " << endl;
    cout << "                               Muhammad Awais, Emanuele Coradin, Fabio Cufino, Tommaso Dorigo    " << endl;
    cout << "                           Enrico Lupi, Eleonora Porcu, Jinu Raj, Fredrik Sandin and Mia Tosi    " << endl;
    cout << "         ------------------------------------------------------------------------------------    " << endl;
    cout << endl;
    cout << "         Run parameters: " << endl;
    cout << "         -----------------------------------" << endl;
    cout << "                       L0 neurons: " << snn_in.N_neuronsL[0] << endl;
    cout << "                       L1 neurons: " << snn_in.N_neuronsL[1] << endl;
    cout << "            Connected L0-L1 frac.: " << snn_in.CF01 << endl;
    cout << "            Connected IN-L0 frac.: " << snn_in.CFI0 << endl;
    cout << "            Connected IN-L1 frac.: " << snn_in.CFI1 << endl;
    cout << "                    Track classes: " << N_classes << endl;
    cout << "                     Total events: " << N_events << endl;
    cout << "               Optimization loops: " << N_epochs << endl;

    // Suppress root warnings
    gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");
    gROOT->ProcessLine("gPrintViaErrorHandler = kTRUE;");

    // Histograms definition
    // ---------------------
    TH1F *SelectivityL0 = new TH1F("SelectivityL0", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *SelectivityL1 = new TH1F("SelectivityL1", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *Qmax = new TH1F("Qmax", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HEff = new TH1F("HEff", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH1F *HAcc = new TH1F("HAcc", "", N_epochs, 0.5, 0.5 + N_epochs);
    TH2F *EffMap = new TH2F("EffMap", "", snn_in.N_neurons, -0.5, snn_in.N_neurons - 0.5, N_ev_classes, -0.5, N_ev_classes - 0.5);
    TH2F *EffMap_window = new TH2F("EffMap_window", "", snn_in.N_neurons, -0.5, snn_in.N_neurons - 0.5, N_ev_classes, -0.5, N_ev_classes - 0.5);
    SelectivityL1->SetLineColor(kBlack);
    Qmax->SetLineColor(2);
    HEff->SetMaximum(1.1);
    HEff->SetMinimum(0.);
    HAcc->SetLineColor(kRed);
    HAcc->SetMaximum(1.1);
    HAcc->SetMinimum(0.);
    Qmax->SetMinimum(0.);
    TH1D *HistDelays = new TH1D("HistDelays", "", 50, 0., snn_in.MaxDelay);
    TH2F *HVoidWs = new TH2F("HVoidWs", "", snn_in.N_neurons, -0.5, -0.5 + snn_in.N_neurons, snn_in.N_streams, -0.5, -0.5 + snn_in.N_streams);

    int N_bins = 100;
    char name[50];
    TH1F *HWeight[snn_in.N_neurons * snn_in.N_streams];
    TH1D *HDelay[snn_in.N_neurons * snn_in.N_streams];
    TH1F *HRMSWeight[snn_in.N_neurons];
    TH1F *HMaxWeight[snn_in.N_neurons];
    TH1F *HMinWeight[snn_in.N_neurons];
    TH1F *Efficiency[snn_in.N_neurons * N_ev_classes];
    TH1F *Efficiency_window[snn_in.N_neurons * N_ev_classes];
    TH1F *FakeRate[snn_in.N_neurons];
    TH1F *FakeRate_window[snn_in.N_neurons];
    TH1F *Eff_totL0[N_ev_classes];
    TH1F *Eff_totL1[N_ev_classes];
    TH2D *StreamsS[10];
    TH2D *StreamsB[10];
    TH2D *StreamsN[10];
    TH1F *BestEff[snn_in.N_neurons];
    TH1F *BestFR[snn_in.N_neurons];
    TH1F *BestEtot[snn_in.N_neurons];

    for (int i = 0; i < snn_in.N_neurons * snn_in.N_streams; i++)
    {
        sprintf(name, "HWeight%d", i);
        HWeight[i] = new TH1F(name, name, N_bins, 0., (float)NevPerEpoch);
        sprintf(name, "HDelay%d", i);
        HDelay[i] = new TH1D(name, name, N_bins, 0., (float)NevPerEpoch);
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
    for (int i = 0; i < snn_in.N_neurons * N_ev_classes; i++)
    {
        sprintf(name, "Efficiency%d", i);
        Efficiency[i] = new TH1F(name, name, N_epochs, 0.5, 0.5 + N_epochs);
        sprintf(name, "Efficiency_window%d", i);
        Efficiency_window[i] = new TH1F(name, name, N_epochs, 0.5, 0.5 + N_epochs);
    }
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        sprintf(name, "FakeRate%d", in);
        FakeRate[in] = new TH1F(name, name, N_epochs, 0.5, 0.5 + N_epochs);
        sprintf(name, "FakeRate_window%d", in);
        FakeRate_window[in] = new TH1F(name, name, N_epochs, 0.5, 0.5 + N_epochs);
    }
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        sprintf(name, "BestEff%d", in);
        BestEff[in] = new TH1F(name, name, N_ev_classes, -0.5, -0.5 + N_ev_classes);
        sprintf(name, "BestFR%d", in);
        BestFR[in] = new TH1F(name, name, N_ev_classes, -0.5, -0.5 + N_ev_classes);
        sprintf(name, "BestEtot%d", in);
        BestEtot[in] = new TH1F(name, name, N_ev_classes, -0.5, -0.5 + N_ev_classes);
    }
    for (int ic = 0; ic < N_ev_classes; ic++)
    {
        sprintf(name, "Eff_totL0%d", ic);
        Eff_totL0[ic] = new TH1F(name, name, N_epochs, 0.5, 0.5 + N_epochs);
        sprintf(name, "Eff_totL1%d", ic);
        Eff_totL1[ic] = new TH1F(name, name, N_epochs, 0.5, 0.5 + N_epochs);
        
    }

    for (int i = 0; i < 10; i++)
    {
        sprintf(name, "StreamsS%d", i);
        StreamsS[i] = new TH2D(name, name, (max_angle + Empty_buffer) * N_display, 0., (max_angle + Empty_buffer) * N_display/10. / omega, snn_in.N_InputStreams, 0.5, snn_in.N_InputStreams + 0.5);
        sprintf(name, "StreamsB%d", i);
        StreamsB[i] = new TH2D(name, name, (max_angle + Empty_buffer) * N_display, 0., (max_angle + Empty_buffer) * N_display/10. / omega, snn_in.N_InputStreams, 0.5, snn_in.N_InputStreams + 0.5);
        sprintf(name, "StreamsN%d", i);
        StreamsN[i] = new TH2D(name, name, (max_angle + Empty_buffer) * N_display, 0., (max_angle + Empty_buffer) * N_display/10. / omega, snn_in.N_neurons, 0.5, snn_in.N_neurons + 0.5);
    }

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
    cout << "                 IPSP dt dilation: " << snn_in.IPSP_dt_dilation << endl;
    cout << "         -----------------------------------" << endl;
    cout << endl;

    // Prime the event loop - we continuously sample detector readout and feed inputs to synapses
    // ------------------------------------------------------------------------------------------
    int N_fires[snn_in.N_neurons];
    double LastP[snn_in.N_neurons];
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        N_fires[in] = 0.;
        LastP[in] = 0.;
    }
    int **fired_sum;
    fired_sum = new int*[N_ev_classes];
    for (int i = 0; i < N_ev_classes; ++i) {
        fired_sum[i] = new int[snn_in.N_neurons];
    }

    int random_fire[snn_in.N_neurons];
    int fired_sum_window[N_ev_classes][snn_in.N_neurons];
    int random_fire_window[snn_in.N_neurons];

    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        random_fire[in] = 0;
        random_fire_window[in] = 0;
        for (int ic = 0; ic < N_ev_classes; ic++)
        {
            fired_sum[ic][in] = 0;
            fired_sum_window[ic][in] = 0;
        }
    }
    
    bool not_fired_bgr = true;
    bool not_fired_bgr_L0 = true;
    int atleastonefired = 0;
    int atleastonefired_L0 = 0;

    int gen_sum[N_ev_classes];
    int fired_anyL0[N_ev_classes];
    int fired_anyL1[N_ev_classes];
    int N_train = NevPerEpoch * Train_fraction;
    int N_test  = NevPerEpoch - N_train;
    
    for (int ic = 0; ic < N_ev_classes; ic++)
    {
        gen_sum[ic] = 0;
        fired_anyL0[ic] = 0;
        fired_anyL1[ic] = 0;
    }
    bool doneL0[N_ev_classes];
    bool doneL1[N_ev_classes];
    bool Seen[N_ev_classes][snn_in.N_neurons];
    float selectivityL0  = 0.;
    float selectivityL1  = 0.;
    float selectivityTOT = 0.;
    float averefftotL0 = 0.;
    float averacctotL0 = 0.;
    float averefftotL1 = 0.;
    float averacctotL1 = 0.;
    vector<vector<vector<pair<int, int>>>> History_ev_class(snn_in.N_neurons);

    Eff = new float[snn_in.N_neurons * N_ev_classes];
    Eff_window = new float[snn_in.N_neurons * N_ev_classes];
    float SumofSquaresofWeight[snn_in.N_neurons] = {0};  // sum of squares synaptic weights for each neuron for RMS calc
    float MeanofSquaresofWeight[snn_in.N_neurons] = {0}; // mean of squares of synaptic weights for each neuron for RMS calc
    float MaxWeight[snn_in.N_neurons];
    float MinWeight[snn_in.N_neurons];
    float RMSWeight[snn_in.N_neurons];

    int count_classes[N_ev_classes+1] = {}; 

    SNN snn_best;
    SNN snn_old;

    // Storing parameters subjected to random search
    snn_old.copy_from(snn_in);
    snn_best.copy_from(snn_in);
    
    
    float Q = 0.;
    float Q_L0 = 0.;

    float Q_old = 0.;
    
    Q_best_L0 = 0.;
    SelL0_best = 0.;
    Eff_best_L0 = 0.;
    Acc_best_L0 = 0.;

    Q_best_L1 = 0.;
    SelL1_best = 0.;
    SelTOT_best = 0.;
    Eff_best_L1 = 0.;
    Acc_best_L1 = 0.;

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
    // Read the ROOT file and get the trees
    auto [file, IT, OT, ET] = readRootFile(rootInput);
    

    // End of reading ----------------------------------------------
    
    // Create csv fout file
    //TODO: RETHINK THE CSV
    /*
    ofstream fout;
    char csv_name[80];
    if (file_id_GS == -1)
        sprintf(csv_name, "MODE/CSV/NL0=%d_NL1=%d_NCl=%d_CF01=%.2f_CFI0=%.2f_CFI1=%.2f_alfa=%.2f_output.csv", snn_in.N_neuronsL[0], snn_in.N_neuronsL[1], N_ev_classes, snn_in.CF01, snn_in.CFI0, snn_in.CFI1, snn_in.alpha);
    else
    {
        sprintf(csv_name, "%i.csv", file_id_GS);
    }
    fout.open(csv_name);
    fout << "Event,ID,Stream,Time,Pclass" << endl;
    */
    // Loop on events ----------------------------------------------
    bool insert = true;
    do
    {
        iev_thisepoch++;

        if (doprogress)
        {
            if (iev_thisepoch % block == 0)
            {
                cout << progress[currchar] << flush;
                currchar++;
            }
        }

        if (ievent % NROOT == 0)
        {
            last_row_event_IT = 0;
            last_row_event_OT = 0;
        }

        ReadFromProcessed(IT, OT, ET, ievent % NROOT + 1);

        // See if we find with track with positive latency by at least one neuron
        for (int in = 0; in < snn_in.N_neurons; in++)
        {
            Seen[pclass][in] = false; // Becomes true if the neuron in has fired for the class pclass
        }
        doneL0[pclass] = false;
        doneL1[pclass] = false;
        not_fired_bgr = true;
        not_fired_bgr_L0 = true;

        // Encode hits in spike streams
        // Here we encode the position of hits through the timing of a spike,
        // In the future we might think at how the analog charge readout in the hits could also be added
        double previous_firetime = 0;

        PreSpike_Time.clear();
        PreSpike_Stream.clear();
        PreSpike_Signal.clear();
        PreSpike_Class.clear();

        double t_in = ievent * (max_angle + Empty_buffer) / omega; // Initial time -> every event adds 25 ns
        Encode(t_in);

        // Keep track of latency calc for each neuron in this event
        bool not_filled[snn_in.N_neurons];
        for (int in = 0; in < snn_in.N_neurons; in++)
        {
            not_filled[in] = true;
        }

        // Loop on spikes and modify neuron and synapse potentials
        // -------------------------------------------------------
        for (int ispike = 0; ispike < PreSpike_Time.size(); ispike++)
        {
            // By looping to size(), we can insert along the way and still make it to the end
            double t = PreSpike_Time[ispike];

            // Modify neuron potentials based on synapse weights
            // -------------------------------------------------
            double min_fire_time = snn_in.largenumber; // if no fire, neuron_firetime returns largenumber
            int in_first = -1;

            // Loop on neurons, but not in order to not favor any neuron
            // ---------------------------------------------------------

            // Shuffle order
            auto rng = default_random_engine{};
            shuffle(neurons_index.begin(), neurons_index.end(), rng);

            for (auto in : neurons_index)
            {
                // Compute future fire times of neurons and their order
                double fire_time = snn_in.Neuron_firetime(in, t);

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
                double latency = 0.;
                
                // Learn weights with spike-time-dependent plasticity: long-term synaptic potentiation
                if(iev_thisepoch < N_train){
                    snn_in.LTD_weights(in_first, min_fire_time, nearest_spike_approx_weights, snn_old);
                    snn_in.LTP_weights(in_first, min_fire_time, nearest_spike_approx_weights, snn_old);
                    // DEcomment for delay 
                    snn_in.LTD_delays(in_first, min_fire_time, nearest_spike_approx_delays, snn_old);
                    snn_in.LTP_delays(in_first, min_fire_time, nearest_spike_approx_delays, snn_old);
                } 
                //we are in the test phase
                else{
                    //loop back on the history of spikes to check if the neuron is reacting to noise or signal
                    auto PreActivation_History = snn_in.Inspect_History(in_first, min_fire_time, window);

                    // Use push_back() to add the retrieved history as a new activation entry
                    History_ev_class[in_first].push_back(PreActivation_History);

                }

                N_fires[in_first]++;
                snn_in.Fire_time[in_first].push_back(min_fire_time);
                
                // Reset history of this neuron
                snn_in.Activate_Neuron(in_first, min_fire_time);
                
                // IPSP for all others at relevant layer
                for (int in2 = 0; in2 < snn_in.N_neurons; in2++)
                {
                    if (in2 != in_first)
                    {
                        if (snn_in.Neuron_layer[in2] == snn_in.Neuron_layer[in_first])
                        { // inhibitions within layer or across
                            snn_in.insert_spike(in2, min_fire_time, snn_in.IPSP, snn_in.N_InputStreams + in_first, snn_in.NOCLASS, iev_thisepoch);
                        }
                    }
                }

                // Create EPS signal in L0 neuron-originated streams
                if (snn_in.Neuron_layer[in_first] == 0)
                { // this is a Layer-0 neuron
                    for (int in = snn_in.N_neuronsL[0]; in < snn_in.N_neurons; in++)
                    {
                        int is = snn_in.N_InputStreams + in_first;
                        if(!snn_in.Void_weight[in][is]){
                            snn_in.insert_spike(in, min_fire_time + snn_in.Delay[in][is], snn_in.EPSP, is, snn_in.NOCLASS, iev_thisepoch);

                        }
                    }
                }

                // Fill spikes train histogram
                if (ievent >= N_events - N_display)
                {
                    int is = (ievent - N_events + N_display) / (N_display/10);
                    double time = min_fire_time - (max_angle + Empty_buffer) / omega * (ievent / (N_display/10)) * (N_display/10);
                    StreamsN[is]->Fill(time, in_first + 1);
                    if (N_part > 0)
                    {
                        //fout << ievent << ", " << 2 << ", " << in_first << "," << time << "," << pclass << endl;
                    }
                }

                // Fill latency histogram
                if (N_part > 0)
                {
                    // How long did it take the first neuron to fire with respect to the arrival time of the first hit?
                    latency = min_fire_time - t_in - First_angle / omega;
                    if (latency >= 0. && not_filled[in_first])
                    {
                        Seen[pclass][in_first] = true;
                        not_filled[in_first] = false;
                    }
                }
                else
                {
                    if (not_filled[in_first] && iev_thisepoch >= N_train)
                    {
                        random_fire[in_first]++;
                        not_filled[in_first] = false;
                    }
                    if (in_first >= snn_in.N_neuronsL[0] && iev_thisepoch >= N_train)
                    { // for Q-value calculations
                        if (not_fired_bgr)
                        {
                            atleastonefired++;
                            not_fired_bgr = false;
                        }
                    }
                    else if (in_first <= snn_in.N_neuronsL[0] && iev_thisepoch >= N_train)
                    { // for Q-value calculations
                        if (not_fired_bgr_L0)
                        {
                            atleastonefired_L0++;
                            not_fired_bgr_L0 = false;
                        }
                    }
                    
                }

                ispike -= 1;
                insert = false;

                // take a step back and search for another activation
                previous_firetime = min_fire_time;

            } // end if in_first fires
            // insert the new spike for the next iteration
            if (insert)
            {
                // Save information on hit-based streams for last N_display events to histograms
                if (ievent >= N_events - N_display)
                {
                    // dividing N_events in 10 groups
                    int is = (ievent - N_events + N_display) / (N_display/10);
                    // time = tin + thit - tin(First event of the group)
                    double time = PreSpike_Time[ispike] - (max_angle + Empty_buffer) / omega * (ievent / (N_display/10)) * (N_display/10);

                    // Histograms
                    if (PreSpike_Signal[ispike] == SIG)
                    {
                        StreamsS[is]->Fill(time, PreSpike_Stream[ispike] + 1);
                        //fout << ievent << ", " << PreSpike_Signal[ispike] << ", " << PreSpike_Stream[ispike] + 1 << "," << time << "," << pclass << endl;
                    }
                    else if (PreSpike_Signal[ispike] == BGR)
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
                        // All input spikes lead to EPSP
                        snn_in.insert_spike(in, t+ snn_in.Delay[in][is], snn_in.EPSP, is, PreSpike_Class[ispike], iev_thisepoch);
                    }
                }
            }
        } // end ispike loop, ready to start over

        // Fill info for efficiency calculations
        if (N_part > 0)
        {
            // efficiency calculations only in the last 10% of the epoch
            if (iev_thisepoch >= N_train)
            {
                //second stat way
                
                //TODO: MODIFY STAT
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
            int bin = (int)(N_bins * (float)iev_thisepoch / NevPerEpoch);
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
                        HDelay[in * snn_in.N_streams + is]->SetBinContent(bin, snn_in.Delay[in][is]);
                        SumofSquaresofWeight[in] += snn_in.Weight[in][is] * snn_in.Weight[in][is]; // for RMS calculation
                        if (snn_in.Weight[in][is] > MaxWeight[in])
                            MaxWeight[in] = snn_in.Weight[in][is]; // finding maxima
                        else if (snn_in.Weight[in][is] < MinWeight[in])
                            MinWeight[in] = snn_in.Weight[in][is]; // finding minima
                    }
                    else{
                        HWeight[in * snn_in.N_streams + is]->SetBinContent(bin, -1);
                        HDelay[in * snn_in.N_streams + is]->SetBinContent(bin, -1);
                    }
                        
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
        // PreSpike_Class  -> particle class of the hit

        // Fill efficiency histograms every NevPerEpoch events, compute Q value and Selectivity, modify parameters
        // ---------------------------------------------------------------------------------------------------

        if (iev_thisepoch == NevPerEpoch)
        { // we did NevPerEpoch events
            //classic efficiency calculations
            // Reset counter that inhibits efficiency and Q calculations until we reach steady state with weights
            iev_thisepoch = 0;
            iepoch++;
            // End of progress bar
            if (doprogress)
                cout << progress[51] << endl;

            for (int in = 0; in < snn_in.N_neurons; in++)
            {
                for (int ic = 0; ic < N_ev_classes; ic++)
                {   
                    int combind = ic + N_ev_classes * in;
                    Eff[combind] = fired_sum[ic][in];
                    if (gen_sum[ic] > 0)
                        Eff[combind] /= gen_sum[ic];
                    Efficiency[combind]->SetBinContent(iepoch, Eff[combind]);
                }
                
                float fakerate = random_fire[in] * 2. / NevPerEpoch / (1.-Train_fraction); // there are NevPerEpoch/2 events with no tracks, where we compute random_fire per neuron
                FakeRate[in]->SetBinContent(iepoch, fakerate);
                
            }
            max_fake = *max_element(random_fire, random_fire + (sizeof(random_fire) / sizeof(random_fire[0]))) * 2. / NevPerEpoch / (1.-Train_fraction);
            float Efftot[N_ev_classes];
            float Efftot_L0[N_ev_classes];
            
            for (int ic = 0; ic < N_ev_classes; ic++)
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
                Efftot_L0[ic] = etl0;
            }

            selectivityL0 = computeMutualInformation(0, gen_sum, fired_sum, snn_in, false);
            SelectivityL0-> Fill(iepoch, selectivityL0);
            selectivityL1 = computeMutualInformation(1, gen_sum, fired_sum, snn_in, false);
            SelectivityL1-> Fill(iepoch, selectivityL1);
            selectivityTOT = computeMutualInformation(-1, gen_sum, fired_sum, snn_in, false);

            // Q value is average efficiency divided by sqrt (aver eff plus aver acceptance)
            // -----------------------------------------------------------------------------
            averacctotL1 = atleastonefired *    (2. / NevPerEpoch / (1.-Train_fraction)); // total acceptance, computed with N_Test*NevPerEpoch/2 events with no tracks
            averacctotL0 = atleastonefired_L0 * (2. / NevPerEpoch / (1.-Train_fraction));
            
            for (int ic = 0; ic < N_ev_classes; ic++)
            {
                averefftotL1 += Efftot[ic];
                averefftotL0 += Efftot_L0[ic];
            }
            averefftotL1 /= N_ev_classes;
            averefftotL0 /= N_ev_classes;

            Q = Compute_Q(averefftotL1, averacctotL1, selectivityL1);
            Q_L0 = Compute_Q(averefftotL0, averacctotL0, selectivityL0);

            //--------------- New method to calculate efficiency, fake rate, Q value ----------------ù
            
            for(int in = 0; in < snn_in.N_neurons; in++){
                vector<vector<bool>> Check_class(N_test, vector<bool>(N_ev_classes, false));
                
                for (int iactivation = 0; iactivation < History_ev_class[in].size(); ++iactivation) {
                    bool fake_fire = true;

                    for (auto&& ev_class : History_ev_class[in][iactivation]) {
                        int event_index = ev_class.first - N_train;
                        int id_class = ev_class.second;

                        // Check if the event index is within the valid range
                        if (event_index < 0 || event_index >= N_test) {
                            cout << "Warning: event_index " << event_index << " is out of range for neuron " << in << endl;
                            continue; // Skip to the next event if it's out of range
                        }

                        if(id_class >= snn_in.SIGCLASS){
                            if (id_class >= N_classes) id_class-=N_classes;
                            fake_fire = false;
                            //If it's the first time that he's firing for this class
                            if(!Check_class[event_index][id_class]){
                                //Update the flag
                                Check_class[event_index][id_class] = true;
                                //Increment for the efficiency
                                fired_sum_window[id_class][in]++;
                            }
                        }                    
                    }
                    //If it didn't receive any particle hit in the window -> increment Fake Fire
                    if(fake_fire) random_fire_window[in]++;
                }

                cout << "Efficiency calculation with window method: " << in << endl;    
                //produce the metrics
                int total_fire = 0;
                for (int ic = 0; ic < N_ev_classes; ic++)
                {
                    cout << ic << ": " << fired_sum_window[ic][in] << endl;
                    total_fire+=fired_sum_window[ic][in];
                    int combind = ic + N_ev_classes * in;
                    Eff_window[combind] = fired_sum_window[ic][in];
                    if (gen_sum[ic] > 0)
                        Eff_window[combind] /= gen_sum[ic];
                    Efficiency_window[combind]->SetBinContent(iepoch, Eff_window[combind]);
                }
                float FP_rate = 0;
                total_fire+=random_fire_window[in];
                if(total_fire>0) FP_rate = 1.*random_fire_window[in]/total_fire;
                FakeRate_window[in]->SetBinContent(iepoch, FP_rate);

                cout << "Neuron " << in << ":" << endl;
                cout << "   - FP rate : " << FP_rate * 100 << endl;

            }

            // ------------------------------------------------------------------------
 
            // Re-initialize neurons
            snn_in.Init_neurons(iev_thisepoch);
            // Reset hits
            Reset_hits();

            /*
            ACHTUNG: SHOULDN'T THIS BE PUT AFTER ALL?

            // Reset weights to initial conditions before new investigation
            snn_in.Reset_weights();
            // Init delays
            if (!false && !ReadPars && !learnDelays)
                snn_in.Init_delays_uniform(); // This unlike void connections, because we can opt to learn these at each cycle too
            */

            cout << "         Ev. # " << ievent + 1 << "; Selectivity L0 = " << selectivityL0 << " L1 = " << selectivityL1 << " totale: " << selectivityTOT
                 << "; Eff L0 = " << averefftotL0 << " Acc L0 = " << averacctotL0 << "; Eff L1 = " << averefftotL1 << " Acc L1 = " << averacctotL1 <<  "Acc max: " << max_fake << "; Firings: ";

            cout << "Total firings per neuron: ";
            for (int in = 0; in < snn_in.N_neurons; in++)
            {
                cout << N_fires[in] << " ";
            }
            cout << endl;

            // Fill debugging graphs of Q as a function of parameters
            int ibin, jbin;
            float n, cont, newcont;
            if (iepoch > 1 || N_epochs == 1)
            { // only do it from end of second epoch onwards, as we are filling delta values
                // Now graph of mean vs sqm of Delay distribution
                double meanDelay = 0.;
                double sqmDelay = 0.;
                double meanOldDelay = 0.;
                double sqmOldDelay = 0.;
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
            }

            // Fill histograms with delays
            HistDelays->Reset();
            HVoidWs->Reset();
            for (int in = 0; in < snn_in.N_neurons; in++)
            {
                for (int is = 0; is < snn_in.N_streams; is++)
                {
                    HistDelays->Fill(snn_in.Delay[in][is]);
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
            cout << "         Q = " << Q << " Old = " << Q_old << " Best = " << Q_best_L1 << endl;

            // Update histograms with current parameter values and optimization metrics
            if (Q > Q_best_L1)
            {
                ind_qbest = iepoch;
                Q_best_L1 = Q;
                Q_best_L0 = Q_L0;

                SelL0_best = selectivityL0;
                Eff_best_L0 = averefftotL0;
                Acc_best_L0 = averacctotL0;

                SelL1_best = selectivityL1;
                Eff_best_L1 = averefftotL1;
                Acc_best_L1 = averacctotL1;

                SelTOT_best = selectivityTOT;

                snn_best.Threshold[0] = snn_in.Threshold[0];
                snn_best.Threshold[1] = snn_in.Threshold[1];
                snn_best.alpha = snn_in.alpha;
                snn_best.L1inhibitfactor = snn_in.L1inhibitfactor;
                snn_best.K = snn_in.K;
                snn_best.K1 = snn_in.K1;
                snn_best.K2 = snn_in.K2;
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
            Qmax->SetBinContent(iepoch, Q_best_L1);
            HEff->SetBinContent(iepoch, averefftotL1);
            HAcc->SetBinContent(iepoch, averacctotL1);
            TCanvas *CU = new TCanvas("CU", "", 1600, 700);
            CU->Divide(5, 2);
            CU->cd(1);
            Qmax->Draw("SAME");
            CU->cd(2);
            HEff->Draw();
            HAcc->Draw("SAME");
            CU->cd(3);
            CU->cd(8);
            HistDelays->Draw();
            CU->cd(9);
            HVoidWs->Draw("COL4");
            CU->cd(10);
            float h = 0;
            float hmax = -snn_in.largenumber;
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

            if (ievent < N_events - 1)
            { // Otherwise we graciously exit loop
                // Reset a few counters
                for (int in = 0; in < snn_in.N_neurons; in++)
                {
                    N_fires[in] = 0.;
                }
                for (int in = 0; in < snn_in.N_neurons; in++)
                {
                    random_fire[in] = 0;
                    for (int ic = 0; ic < N_ev_classes; ic++)
                    {
                        fired_sum[ic][in] = 0;
                    }
                }
                for (int ic = 0; ic < N_ev_classes; ic++)
                {
                    gen_sum[ic] = 0;
                    fired_anyL0[ic] = 0;
                    fired_anyL1[ic] = 0;
                }
                atleastonefired = 0;
                atleastonefired_L0 = 0;

                not_fired_bgr = true;
                not_fired_bgr_L0 = true;

                // Reset progress bar
                if (doprogress)
                {
                    cout << "         " << progress[0];
                    currchar = 1;
                }
            }

        }         // if ievent+1%NevPerEpoch = 0
        ievent++; // only go to next event if we did a backward pass too
    } while (ievent < N_events);
    // closing the input file
    delete IT;
    delete OT;

    //fout.close();
    file->Close();
    delete file;

    // Draw histograms
    cout << "Saving histograms" << endl;
    TCanvas *S = new TCanvas("S", "", 3000, 600);
    S->Divide(5, 4);

    for (int irow = 0; irow < 4; irow++)
    {
        for (int icol = 0; icol < 5; icol++)
        {
            S->cd(icol + irow * 5 + 1);
            if (irow % 2 == 0)
            {
                StreamsB[icol + (irow) * 5 / 2]->SetLineColor(kRed);
                StreamsB[icol + (irow) * 5 / 2]->Draw("BOX");
                StreamsS[icol + (irow) * 5 / 2]->SetLineColor(kBlue);
                StreamsS[icol + (irow) * 5 / 2]->Draw("BOXSAME");
            }
            else
            {
                StreamsN[icol + (irow - 1) * 5 / 2]->SetLineColor(kGreen);
                StreamsN[icol + (irow - 1) * 5 / 2]->Draw("BOX");
            }
        }
    }

    TCanvas *E0 = new TCanvas("E0", "", 800, 800);
    E0->Divide(N_ev_classes, snn_in.N_neuronsL[0]);
    for (int i = 0; i < snn_in.N_neuronsL[0] * N_ev_classes; i++)
    {
        E0->cd(i + 1);
        Efficiency[i]->SetMaximum(1.1);
        Efficiency[i]->SetMinimum(0.);
        Efficiency[i]->Draw("");
        
        int in = i / N_ev_classes;
        FakeRate[in]->SetMarkerColor(2);
        FakeRate[in]->SetLineColor(2);
        FakeRate[in]->Draw("SAME");
        int ic = i % N_ev_classes;
        Eff_totL0[ic]->SetMarkerColor(3);
        Eff_totL0[ic]->SetLineColor(3);
        Eff_totL0[ic]->Draw("SAME");
        Efficiency[i]->Draw("SAME");
    }

    TCanvas *E1 = new TCanvas("E1", "", 800, 800);
    E1->Divide(N_ev_classes, snn_in.N_neuronsL[1]);
    for (int i = snn_in.N_neuronsL[0] * N_ev_classes; i < snn_in.N_neurons * N_ev_classes; i++)
    {
        E1->cd(i + 1 - snn_in.N_neuronsL[0] * N_ev_classes);
        Efficiency[i]->SetMaximum(1.1);
        Efficiency[i]->SetMinimum(0.);
        Efficiency[i]->Draw("");
        int in = i / N_ev_classes;
        FakeRate[in]->SetMarkerColor(2);
        FakeRate[in]->SetLineColor(2);
        FakeRate[in]->Draw("SAME");
        int ic = i % N_ev_classes;
        Eff_totL1[ic]->SetMarkerColor(3);
        Eff_totL1[ic]->SetLineColor(3);
        Eff_totL1[ic]->Draw("SAME");
        Efficiency[i]->Draw("SAME");
    }

    TCanvas *E0_window = new TCanvas("E0_window", "", 800, 800);
    E0_window->Divide(N_ev_classes, snn_in.N_neuronsL[0]);
    for (int i = 0; i < snn_in.N_neuronsL[0] * N_ev_classes; i++)
    {
        E0_window->cd(i + 1);
        Efficiency_window[i]->SetMaximum(1.1);
        Efficiency_window[i]->SetMinimum(0.);
        Efficiency_window[i]->Draw("");
        
        int in = i / N_ev_classes;
        FakeRate_window[in]->SetMarkerColor(2);
        FakeRate_window[in]->SetLineColor(2);
        FakeRate_window[in]->Draw("SAME");
    }

    TCanvas *E1_window = new TCanvas("E1_window", "", 800, 800);
    E1_window->Divide(N_ev_classes, snn_in.N_neuronsL[1]);
    for (int i = snn_in.N_neuronsL[0] * N_ev_classes; i < snn_in.N_neurons * N_ev_classes; i++)
    {
        E1_window->cd(i + 1 - snn_in.N_neuronsL[0] * N_ev_classes);
        Efficiency_window[i]->SetMaximum(1.1);
        Efficiency_window[i]->SetMinimum(0.);
        Efficiency_window[i]->Draw("");
        int in = i / N_ev_classes;
        FakeRate_window[in]->SetMarkerColor(2);
        FakeRate_window[in]->SetLineColor(2);
        FakeRate_window[in]->Draw("SAME");
    }

    // Plot the efficiencies and acceptances for the best q-value run
    // --------------------------------------------------------------
    for (int i = 0; i < snn_in.N_neurons * N_ev_classes; i++)
    {
        int in = i / N_ev_classes;
        int ic = i % N_ev_classes;
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
            HWeight[in * snn_in.N_streams + is]->SetMaximum((HMaxWeight[in]->GetMaximum()) * 1.1);
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

    TCanvas *D = new TCanvas("D", "", 500, 800);

    D->Divide(2, nrow);
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        D->cd(in + 1);
        for (int is = 0; is < snn_in.N_streams; is++)
        {
            HDelay[in * snn_in.N_streams + is]->SetMaximum(snn_in.MaxDelay * 1.1);
            HDelay[in * snn_in.N_streams + is]->SetMinimum(0.);
            int color = is + 1;
            if (color == 8)
                color = 9;
            HDelay[in * snn_in.N_streams + is]->SetLineColor(color);
            if (is == 0)
            {
                HDelay[in * snn_in.N_streams + is]->Draw();
            }
            else
            {
                HDelay[in * snn_in.N_streams + is]->Draw("SAME");
            }
        }
    }


    // Draw final Efficiency and acceptance maps
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        for (int ic = 0; ic < N_ev_classes; ic++)
        {
            EffMap->SetBinContent(in + 1, ic + 1, Efficiency[ic + in * N_ev_classes]->GetBinContent(ind_qbest));
        }
    }
    TCanvas *Y = new TCanvas("Y", "", 600, 900);
    Y->cd();
    EffMap->Draw("COL4");

    // Draw final Efficiency and acceptance maps
    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        for (int ic = 0; ic < N_ev_classes; ic++)
        {
            EffMap_window->SetBinContent(in + 1, ic + 1, Efficiency_window[ic + in * N_ev_classes]->GetBinContent(ind_qbest));
        }
    }
    TCanvas *Y_window = new TCanvas("Y_window", "", 600, 900);
    Y_window->cd();
    EffMap_window->Draw("COL4");

    //Draw delays

    // Create a 2D histogram to store the delays
    TH2D* delayHistogram = new TH2D("Delay Histogram", "Delay Histogram", snn_in.N_neurons, 0, snn_in.N_neurons, snn_in.N_streams, 0, snn_in.N_streams);

    // Fill the histogram with delay values
    for (int in = 0; in < snn_in.N_neurons; in++) {
        for (int is = 0; is < snn_in.N_streams; is++) {
            delayHistogram->SetBinContent(in + 1, is + 1, snn_in.Delay[in][is]);
        }
    }

    // Create a canvas to draw the histogram
    TCanvas* delay_canvas = new TCanvas("Delay", "Delay", 800, 600);
    delay_canvas->cd();
    gStyle->SetOptStat(0);
    delayHistogram->SetOption("goff");
    delayHistogram->Draw("colz"); // Draw the histogram with a color palette
    
    //Draw delta delays
    // Create a 2D histogram to store the delays
    TH2D* delta_delayHistogram = new TH2D("Delta Delay Histogram", "Delta Delay Histogram", snn_in.N_neurons, 0, snn_in.N_neurons, snn_in.N_streams, 0, snn_in.N_streams);

    // Fill the histogram with delay values
    for (int in = 0; in < snn_in.N_neurons; in++) {
        for (int is = 0; is < snn_in.N_streams; is++) {
            delta_delayHistogram->SetBinContent(in + 1, is + 1, snn_in.Delay[in][is] - snn_in.Delay_initial[in][is]);
        }
    }

    // Create a canvas to draw the histogram
    TCanvas* delta_delay_canvas = new TCanvas("Delta Delay", "Delta Delay", 800, 600);
    delta_delay_canvas->cd();
    gStyle->SetOptStat(0);
    delta_delayHistogram->SetOption("goff");
    delta_delayHistogram->Draw("colz"); // Draw the histogram with a color palette


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
    cout << "               Average efficiency L0: " << Eff_best_L0 << endl;
    cout << "                Average fake rate L0: " << Acc_best_L0 << endl;
    cout << "                  Maximum Q value L0: " << Q_best_L0 << endl;
    cout << "                   L0 selectivity: " << SelL0_best << endl;
    cout << endl;
    cout << "               Average efficiency L1: " << Eff_best_L1 << endl;
    cout << "                Average fake rate L1: " << Acc_best_L1 << endl;
    cout << "                  Maximum Q value L1: " << Q_best_L1 << endl;
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
    cout << "                 IPSP dt dilation: " << snn_best.IPSP_dt_dilation << endl;
    cout << "         -----------------------------------" << endl;
    cout << endl;

    // Dump to file optimized parameters and results
    // ---------------------------------------------

    // Dump SNN parameters to JSON file
    string Path = SNN_PATH + "/Code/MODE/JSON/";
    stringstream sstr;
    char num[80];

    if (file_id_GS == -1)
        sprintf(num, "NL0=%d_NL1=%d_NCl=%d_CF01=%.2f_CFI0=%.2f_CFI1=%.2f_alfa=%.2f_%d", snn_in.N_neuronsL[0], snn_in.N_neuronsL[1], N_ev_classes, snn_in.CF01, snn_in.CFI0, snn_in.CFI1, snn_in.alpha, indfile);
    else
        sprintf(num, "%i", file_id_GS);

    sstr << "Parameters_";
    string namejsonfile = Path + sstr.str() + num + ".json";
    
    snn_in.dumpToJson(namejsonfile);
    appendToJson(namejsonfile, Eff_best_L1, max_fake, SelTOT_best, Q_best_L1);
    
    // Dump histograms to root file
    Path = SNN_PATH + "/Code/MODE/SNNT/";
    sstr.str(string());
    
    sstr << "Histos13_";
    string namerootfile = Path + sstr.str() + num + ".root";
    TFile *rootfile = new TFile(namerootfile.c_str(), "RECREATE");
    rootfile->cd();
    
    // Write canvases first
    cout << "Writing canvas" << endl;
    S->Write();
    E0->Write();
    E1->Write();
    BE->Write();
    SE->Write();
    W->Write();
    D->Write();
    Y->Write();
    MW->Write();
    E0_window->Write();
    E1_window->Write();
    Y_window->Write();
    delay_canvas->Write();
    delta_delay_canvas->Write();

    cout << "Saving pdf" << endl;
    S->SaveAs((SNN_PATH + "/Code/pdf/S.pdf").c_str(), "pdf");
    E0->SaveAs((SNN_PATH + "/Code/pdf/E0.pdf").c_str(), "pdf");
    E1->SaveAs((SNN_PATH + "/Code/pdf/E1.pdf").c_str(), "pdf");
    BE->SaveAs((SNN_PATH + "/Code/pdf/BE.pdf").c_str(), "pdf");
    SE->SaveAs((SNN_PATH + "/Code/pdf/SE.pdf").c_str(), "pdf");
    W->SaveAs((SNN_PATH + "/Code/pdf/W.pdf").c_str(), "pdf");
    D->SaveAs((SNN_PATH + "/Code/pdf/D.pdf").c_str(), "pdf");
    Y->SaveAs((SNN_PATH + "/Code/pdf/Y.pdf").c_str(), "pdf");
    MW->SaveAs((SNN_PATH + "/Code/pdf/MWlog.pdf").c_str(), "pdf");
    delay_canvas->SaveAs((SNN_PATH + "/Code/pdf/delay_canvas.pdf").c_str(), "pdf");
    delta_delay_canvas->SaveAs((SNN_PATH + "/Code/pdf/delta_delay_canvas.pdf").c_str(), "pdf");

    // Then histograms
    SelectivityL0->Write();
    SelectivityL1->Write();
    Qmax->Write();
    HEff->Write();
    HAcc->Write();
    delayHistogram->Write();
    delta_delayHistogram->Write();
    HVoidWs->Write();
    HistDelays->Write();

    for (int in = 0; in < snn_in.N_neurons; in++)
    {
        for (int ic = 0; ic < N_ev_classes; ic++)
        {
            int id = in * N_ev_classes + ic;
            Efficiency[id]->Write();
            Efficiency_window[id]->Write();
        }
        for (int is = 0; is < snn_in.N_streams; is++)
        {
            int id = in * snn_in.N_streams + is;
            HWeight[id]->Write();
            HDelay[id]->Write();
        }
        FakeRate[in]->Write();
        FakeRate_window[in]->Write();
        BestEff[in]->Write();
        BestFR[in]->Write();
        BestEtot[in]->Write();
        HRMSWeight[in]->Write();
        HMaxWeight[in]->Write();
        HMinWeight[in]->Write();
    }
    for (int ic = 0; ic < N_ev_classes; ic++)
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
    EffMap_window->Write();

    MW->Write();
    rootfile->Write();

    // End of program
    rootfile->Close();
    gROOT->Time();

    // Deallocate memory to prevent memory leaks
    for (int i = 0; i < N_ev_classes; ++i) {
        delete[] fired_sum[i]; // Delete each row
    }
    delete[] fired_sum;
    delete SelectivityL0;
    delete SelectivityL1;
    delete Qmax;
    delete HEff;
    delete HAcc;
    delete EffMap;
    delete EffMap_window;
    delete HistDelays;
    delete HVoidWs;

    for (int i = 0; i < snn_in.N_neurons * snn_in.N_streams; ++i) {
        delete HWeight[i];
        delete HDelay[i];
    }
    for (int i = 0; i < snn_in.N_neurons * N_ev_classes; ++i) {
        delete Efficiency[i];
        delete Efficiency_window[i];
    }
    for (int i = 0; i < snn_in.N_neurons; ++i) {
        delete HRMSWeight[i];
        delete HMaxWeight[i];
        delete HMinWeight[i];
        delete FakeRate[i];
        delete FakeRate_window[i];
        delete BestEff[i];
        delete BestFR[i];
        delete BestEtot[i];
    }
    for (int i = 0; i < N_ev_classes; ++i) {
        delete Eff_totL0[i];
        delete Eff_totL1[i];
    }
    for (int i = 0; i < 10; ++i) {
        delete StreamsS[i];
        delete StreamsB[i];
        delete StreamsN[i];
    }



    return;
}

void PrintHelp()
{
    cout << "This is the list of the available parameters that you can pass to the program:" << endl
         << endl;

    cout << "SNN parameters" << endl;
    cout << "   --NL0:           number of neurons in layer 0" << endl;
    cout << "   --NL1:           number of neurons in layer 1" << endl;
    cout << "   --CF01:          average fraction of active connections between layer 0 and 1. Default 1. Float" << endl;
    cout << "   --CFI0:          average fraction of active connections between the tracker layer 0. Default 1. Float" << endl;
    cout << "   --CFI1:          average fraction of active connections between the tracker layer 1. Default 1. Float" << endl;
    cout << "   --split_layer0:  possibility of splitting layer 0 in two indipendent parts with different lateral inhibitions. Default false. Bool" << endl;
    cout << "   --exclude_L0:    number of neurons of layer 0 to randomly can be turned off. Default 0. Integer" << endl;
    cout << "   --exclude_L1:    number of neurons of layer 0 to randomly can be turned off. Default 0. Integer" << endl;

    cout << endl;

    cout << "Neuron related parameters" << endl;
    cout << "   --alpha:            Inhibition strength. Float." << endl;
    cout << "   --L1inhibitfactor:  To further amplify the inhibition strength in layer 1. Float." << endl;
    cout << "   --K:                Regulate EPSP amplitude." << endl;
    cout << "   --K1:               Activation potential shape constant" << endl;
    cout << "   --K2:               Activation potential shape constant" << endl;
    cout << "   --IPSP_dt_dilation: Regulate the inhibition duration." << endl;
    cout << "   --MaxDelay:         Maximum delay" << endl;
    cout << "   --tau_m:            Membrane time constant" << endl;
    cout << "   --tau_s:            Synapse time constant" << endl;
    cout << "   --tau_r:            Refractory period: time for which the neuron can't fire again after its activation." << endl;
    cout << "   --TH0:              Threshold for  theneurons belonging to layer 0" << endl;
    cout << "   --TH1:              Threshold for  theneurons belonging to layer 1" << endl;
    cout << "   --sparsity:         Regulate the spreadness of the initialization of the delays" << endl;

    cout << endl;

    cout << "Learning related parameters" << endl;
    cout << "   --a_plus:       Regulate the strength of the potentiation of the synaptic weights through STDP." << endl;
    cout << "   --a_minus:      Regulate the strength of the depression of the synaptic weights through STDP." << endl;
    cout << "   --tau_plus:     Time constant of the STDP potentiation of the synaptic weights" << endl;
    cout << "   --tau_minus:    Time constant of the STDP depression of the synaptic weights" << endl;
    cout << "   --d_plus:       Regulate the strength of the potentiation of the synaptic delays through STDP." << endl;
    cout << "   --d_minus:      Regulate the strength of the depression of the synaptic delays through STDP." << endl;
    cout << "   --taud_plus:    First time constant of the STDP potentiation of the synaptic delays" << endl;
    cout << "   --taud_minus:   First time constant of the STDP depression of the synaptic delays" << endl;
    cout << "   --taud_plus_2:  Second time constant of the STDP potentiation of the synaptic delays. Condition: taud_plus_2 < taud_plus" << endl;
    cout << "   --taud_minus_2: Second time constant of the STDP depression of the synaptic delays. Condition: taud_minus_2 < taud_minus" << endl;

    cout << endl;

    cout << "Main program parameters" << endl;
    cout << "   --N_ev:           Number of events" << endl;
    cout << "   --N_ep:           DEPRECATED. Number of epochs" << endl;
    cout << "   --batch:          To run in batch mode" << endl;
    cout << "   --rootInput:      File root containing data. Must start with /Code" << endl;
    cout << "   --N_classes:      Number of different classes of particles" << endl;
    cout << "   --N_ev_classes:   Number of different kinds of events to handle multi-tracks scenarios" << endl;
    cout << "   --ReadPars:       If setted it indicates the path to the json file to load a predefined network" << endl;
    cout << "   --NROOT:          Number of events contained in the root file" << endl;
    cout << "   --N_display:      Number of events to be displayed in the plots" << endl;
    cout << "   --Train_fraction: Fraction of events to be used during the training phase" << endl;
    cout << "   --window:         Set the time-window for the calculation of the efficiency" << endl;           
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// Main routine
// ------------
int main(int argc, char *argv[])
{
    SNN_PATH = string(getenv("SNN_PATH"));
    int file_id_GS = -1;
    SNN *S;
    
    // Pre scan to look for ReadPars flag
    for (int i = 1; i < argc; i++)
    {
        const char *arg = argv[i];
        if (strcmp(arg, "--ReadPars") == 0){
            ReadPars = argv[i + 1];
            cout << "-- WARNING --: you are fetching the network from a previous configuration!" << endl 
                 << "               The effect of the following flags will be ignored: " << endl 
                 << "                   - NL0 "  << endl
                 << "                   - NL1 "  << endl
                 << "                   - CF01 " << endl
                 << "                   - CFI0 " << endl
                 << "                   - CFI1 " << endl
                 << "                   - sparsity " << endl
                 << "                   - split_layer0 " << endl;
            break;
        }
            
    }

    // update the default values
    if (ReadPars != "none") {
        S = new SNN();
        cout << "Retrieving parameters from file" << endl;
        S->loadFromJson(ReadPars);

        // Set loaded values to variables for potential command-line override
        _NL0                = S->N_neuronsL[0];
        _NL1                = S->N_neuronsL[1];
        _alpha              = S->alpha;
        _CFI0               = S->CFI0;
        _CFI1               = S->CFI1;
        _CF01               = S->CF01;
        _L1inhibitfactor    = S->L1inhibitfactor;
        _K                  = S->K;
        _K1                 = S->K1;
        _K2                 = S->K2;
        _IPSP_dt_dilation   = S->IPSP_dt_dilation;
        _MaxDelay           = S->MaxDelay;
        _tau_m              = S->tau_m;
        _tau_s              = S->tau_s;
        _tau_r              = S->tau_r;
        _tau_plus           = S->tau_plus;
        _tau_minus          = S->tau_minus;
        _a_plus             = S->a_plus;
        _a_minus            = S->a_minus;
        _taud_plus          = S->taud_plus;
        _taud_minus         = S->taud_minus;
        _taud_plus_2        = S->taud_plus_2;
        _taud_minus_2       = S->taud_minus_2;
        _d_plus             = S->d_plus;
        _d_minus            = S->d_minus;
        _N_InputStreams     = S->N_InputStreams;
        _Threshold0         = S->Threshold[0];
        _Threshold1         = S->Threshold[1];
        _sparsity           = S->sparsity;
        _split_layer0       = S->split_layer0;
    }

    // Loop through the command-line arguments
    for (int i = 1; i < argc; i++)
    {
        const char *arg = argv[i];
        if (ReadPars == "none") {
            if (strcmp(arg, "--NL0") == 0)
                _NL0 = stoi(argv[i + 1]);
            else if (strcmp(arg, "--NL1") == 0)
                _NL1 = stoi(argv[i + 1]);
            else if (strcmp(arg, "--CF01") == 0)
                _CF01 = stof(argv[i + 1]);
            else if (strcmp(arg, "--CFI0") == 0)
                _CFI0 = stof(argv[i + 1]);
            else if (strcmp(arg, "--CFI1") == 0)
                _CFI1 = stof(argv[i + 1]);
             else if (strcmp(arg, "--sparsity") == 0)
                _sparsity = stof(argv[i + 1]);
            else if (strcmp(arg, "--split_layer0") == 0)
                _split_layer0 = stoi(argv[i + 1]);
        }
        if (strcmp(arg, "--alpha") == 0)
            _alpha = stof(argv[i + 1]);

        else if (strcmp(arg, "--L1inhibitfactor") == 0)
            _L1inhibitfactor = stof(argv[i + 1]);

        else if (strcmp(arg, "--K") == 0)
            _K = stof(argv[i + 1]);
        else if (strcmp(arg, "--K1") == 0)
            _K1 = stof(argv[i + 1]);
        else if (strcmp(arg, "--K2") == 0)
            _K2 = stof(argv[i + 1]);

        else if (strcmp(arg, "--IPSP_dt_dilation") == 0)
            _IPSP_dt_dilation = stof(argv[i + 1]);

        else if (strcmp(arg, "--MaxDelay") == 0)
            _MaxDelay = stof(argv[i + 1]);

        else if (strcmp(arg, "--tau_m") == 0)
            _tau_m = stof(argv[i + 1]);
        else if (strcmp(arg, "--tau_s") == 0)
            _tau_s = stof(argv[i + 1]);
        else if (strcmp(arg, "--tau_r") == 0)
            _tau_r = stof(argv[i + 1]);
        else if (strcmp(arg, "--tau_plus") == 0)
            _tau_plus = stof(argv[i + 1]);
        else if (strcmp(arg, "--tau_minus") == 0)
            _tau_minus = stof(argv[i + 1]);

        else if (strcmp(arg, "--a_plus") == 0)
            _a_plus = stof(argv[i + 1]);
        else if (strcmp(arg, "--a_minus") == 0)
            _a_minus = stof(argv[i + 1]);

        else if (strcmp(arg, "--taud_plus") == 0)
            _taud_plus = stof(argv[i + 1]);
        else if (strcmp(arg, "--taud_minus") == 0)
            _taud_minus = stof(argv[i + 1]);
        else if (strcmp(arg, "--taud_plus_2") == 0)
            _taud_plus_2 = stof(argv[i + 1]);
        else if (strcmp(arg, "--taud_minus_2") == 0)
            _taud_minus_2 = stof(argv[i + 1]);

        else if (strcmp(arg, "--d_plus") == 0)
            _d_plus = stof(argv[i + 1]);
        else if (strcmp(arg, "--d_minus") == 0)
            _d_minus = stof(argv[i + 1]);

        else if (strcmp(arg, "--TH0") == 0)
            _Threshold0 = stof(argv[i + 1]);
        else if (strcmp(arg, "--TH1") == 0)
            _Threshold1 = stof(argv[i + 1]);

        // Main's parameters
        else if (strcmp(arg, "--N_ev") == 0)
            N_events = stoi(argv[i + 1]);
        else if (strcmp(arg, "--N_ep") == 0)
            N_epochs = stoi(argv[i + 1]);
        else if (strcmp(arg, "--N_display") == 0)
            N_display = stoi(argv[i + 1]);
        else if (strcmp(arg, "--Train_fraction") == 0)
            Train_fraction = stof(argv[i + 1]);
        else if (strcmp(arg, "--batch") == 0)
            batch = stoi(argv[i + 1]);
        else if (strcmp(arg, "--rootInput") == 0)
            rootInput = argv[i + 1];
        else if (strcmp(arg, "--window") == 0)
            window = stof(argv[i + 1]);

        else if (strcmp(arg, "--N_classes") == 0)
            N_classes = stoi(argv[i + 1]);
        else if (strcmp(arg, "--N_ev_classes") == 0)
            N_ev_classes = stoi(argv[i + 1]);
        else if (strcmp(arg, "--NROOT") == 0)
            NROOT = stoi(argv[i + 1]);
        else if (strcmp(arg, "--file_id_GS") == 0)
            file_id_GS = stoi(argv[i + 1]);
        else if (strcmp(arg, "--exclude_L0") == 0)
            _exclude_L0 = stoi(argv[i + 1]);
        else if (strcmp(arg, "--exclude_L1") == 0)
            _exclude_L1 = stoi(argv[i + 1]);
        else if (strcmp(arg, "--help") == 0)
        {
            PrintHelp();
            exit(0);
        }
        
    }
    rootInput = SNN_PATH + rootInput;
    
    if (ReadPars != "none") {
        S = new SNN();
        cout << "Retrieving parameters from file" << endl;
        S->loadFromJson(ReadPars);
        S->Reset_Parameters(_alpha, _L1inhibitfactor, _K, _K1, _K2, _IPSP_dt_dilation, _MaxDelay, _tau_m, _tau_s, _tau_r, _tau_plus, _tau_minus, _a_plus, _a_minus, _taud_plus, _taud_minus, _taud_plus_2, _taud_minus_2, _d_plus, _d_minus, _Threshold0, _Threshold1);
    }
    else{
        S = new SNN(_NL0, _NL1,
          _alpha,
          _CFI0, _CFI1, _CF01,
         _L1inhibitfactor,
         _K,  _K1,  _K2,
         _IPSP_dt_dilation,
         _MaxDelay,

          _tau_m, _tau_s, _tau_r, _tau_plus, _tau_minus,
         _a_plus,  _a_minus,

         _taud_plus, _taud_minus,
         _taud_plus_2, _taud_minus_2,
         _d_plus, _d_minus,

          _N_InputStreams,
          _Threshold0, _Threshold1, _sparsity, _split_layer0);
    }

    if(_exclude_L0 > 0 || _exclude_L1 > 0){
        vector<int> excluded_neurons = generate_excluded_neurons(_exclude_L0, _exclude_L1, S->N_neuronsL[0], S->N_neurons);

        // Exclude neurons in SNN
        S->Exclude_neurons(excluded_neurons);
    }

    SNN_Tracking(*S, file_id_GS);
    cout << "Creating the file for the potentials plot" << endl;
    PlotPotentials(rootInput.c_str(), *S, 12);

    delete S;
    
    return 0;
}

