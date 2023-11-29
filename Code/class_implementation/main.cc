#include <iostream>
#include <random>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include "SNN.h"

static int N_part; // Number of generated particles in an event
static float First_angle;
static long int ievent;
static int N_events;

static vector<float> PreSpike_Time;
static vector<int> PreSpike_Stream;
static vector<int> PreSpike_Signal;

static vector<int> neurons_index;
static bool insert = true;

// clear hits vector
void Reset_hits()
{
    hit_pos.clear();
    return;
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

} // read Weights from root file

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

// plot neuron potentials as a function of time
void PlotPotentials(const char *rootWeight, const char *rootInput, SNN &P, int _N_events, bool read_weights=true, bool no_firing_mode=false)
{
    vector<int>neurons_index;
    // initialization of neurons_index vector
    for (int i = 0; i < P.N_neurons; i++)
        neurons_index.push_back(i);
    // vector used to keep track of the neurons firing order
    vector<int> Fire_ID;
    N_events = _N_events;

    // vectors to plot
    vector<float> Time[N_events];
    vector<float> Potential[N_events][P.N_neurons];
    int fire_count[P.N_neurons];
    for (int ic = 0; ic < P.N_neurons; ic++)
        fire_count[ic] = 0;

    cout << "Initializaing the plot SNN" << endl;
    SNN P_plot(P.N_neuronsL[0], P.N_neuronsL[1]);
    
    cout << "Opening the weight file" << endl;
    TFile *file_weight = TFile::Open(rootWeight, "READ");
    if (!file_weight || file_weight->IsZombie())
    {
        cerr << "Error: Cannot open file " << rootWeight << endl;
        return;
    }
    //Uncomment to read weights
    if(read_weights){
        ReadWeights(file_weight, P);
        ReadWeights(file_weight, P_plot);
    }
    else{
        P.Set_weights();
        P_plot.Set_weights();
    }
    

    // the network is ready
    // we need to fecth the events and compute the plots

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

    int ievent = 1;
    // Loop on events ----------------------------------------------
    do
    {
        float previous_firetime = 0;
        cout << "Event " << ievent << endl;
        PreSpike_Time.clear();
        PreSpike_Stream.clear();
        PreSpike_Signal.clear();
        Fire_ID.clear();
        
        ReadFromProcessed(IT, OT, ievent);

        float t_in = (ievent-1) * (max_angle + Empty_buffer) / omega; // Initial time -> every event adds 25 ns
        Encode(t_in);

        // Loop on spikes and modify neuron and synapse potentials
        // -------------------------------------------------------
        for (int ispike = 0; ispike < PreSpike_Time.size(); ispike++)
        {
            // By looping to size(), we can insert along the way and still make it to the end
            float t = PreSpike_Time[ispike];

            // Modify neuron potentials based on synapse weights
            // -------------------------------------------------
            float min_fire_time = largenumber - 1.; // if no fire, neuron_firetime returns largenumber
            int in_first = -1;

            // Loop on neurons, but not in order to not favor any neuron
            // ---------------------------------------------------------

            // Shuffle order
            auto rng = default_random_engine{};
            shuffle(neurons_index.begin(), neurons_index.end(), rng);

            for (auto in : neurons_index)
            {
                //  We implement a scheme where input streams produce an IE signal into L0, an EPS into L1, and L0 neurons EPS into L1
                //  Add to neuron history, masking out L1 spikes for L0 neurons
                int is = PreSpike_Stream[ispike];
                if (is < N_InputStreams || P.Neuron_layer[in] > 0)
                { // otherwise stream "is" does not lead to neuron "in"
                    if(insert){
                        P.History_time[in].push_back(t);
                        // All input spikes lead to EPSP
                        P.History_type[in].push_back(1);
                        // P.History_type[in].push_back(1);
                        P.History_ID[in].push_back(is);
                    }
                    
                    // Compute future fire times of neurons and their order
                    float fire_time = P.Neuron_firetime_past(in, t);

                    if (fire_time < min_fire_time)
                    {
                        in_first = in;
                        min_fire_time = fire_time;
                    }
                }
            }
            if (in_first == -1){
                //finish the plot for the previous spike
                float t_prime;
                float delta_t;
                if(!insert){
                    t_prime = previous_firetime;
                    delta_t = (t - t_prime)/11;    
                    
                    for(int inc = 0; inc < 11; inc++){
                        Time[ievent-1].push_back(t_prime+inc*delta_t);
                        for (auto in : neurons_index)
                        {
                            Potential[ievent-1][in].push_back(P.Neuron_Potential(in, t_prime+inc*delta_t, false));
                        }
                    }
                    
                }
                // plot the result of the incoming spike
                else
                {
                    if (ispike==0) t_prime=t_in;
                    else t_prime = PreSpike_Time[ispike-1];
                    
                    delta_t =(t - t_prime)/11;
                    for(int inc = 0; inc < 11; inc++){
                        Time[ievent-1].push_back(t_prime+inc*delta_t);
                        for (auto in : neurons_index)
                        {
                            Potential[ievent-1][in].push_back(P.Neuron_Potential(in, t_prime+inc*delta_t, false));
                        }
                    } 
                }

                insert = true;
                
                continue;
            }

            // Ok, neuron in_first is going to fire next.
            // Peek at next event in list, to see if it comes before in_first fires
            // --------------------------------------------------------------------

            if (ispike < PreSpike_Time.size() - 1)
            {
                if (PreSpike_Time[ispike + 1] >= min_fire_time)
                { // otherwise we go to next spike in list
                    float t_prime;
                    float delta_t;
                    // handle firing of neuron in_first
                    //That means that we are handling the first neruon activation bewtween two EPSP spike
                    
                    if(insert){
                        if(ispike==0) t_prime = t_in;
                        else t_prime =  PreSpike_Time[ispike-1];
                    }
                    else t_prime = previous_firetime;
                    delta_t = (min_fire_time - t_prime)/11;    
                    
                    for(int inc = 0; inc < 11; inc++){
                        Time[ievent-1].push_back(t_prime+inc*delta_t);
                        for (auto in : neurons_index)
                        {
                            Potential[ievent-1][in].push_back(P.Neuron_Potential(in, t_prime+inc*delta_t, false));
                        }
                    }
                    
                    
                    cout << "Neuron " << in_first << " is firing -> computing the consequences" << endl;

                    P.Fire_time[in_first].push_back(min_fire_time);
                    Fire_ID.push_back(in_first);
                    // Reset history of this neuron
                    P.History_time[in_first].clear();
                    P.History_type[in_first].clear();
                    P.History_ID[in_first].clear();
                    P.History_time[in_first].push_back(min_fire_time);
                    P.History_type[in_first].push_back(0);
                    P.History_ID[in_first].push_back(0); // ID is not used for type 0 history events

                    // IPSP for all others at relevant layer
                    for (int in2 = 0; in2 < P.N_neurons; in2++)
                    {
                        if (in2 != in_first)
                        {
                            if (P.Neuron_layer[in2] == P.Neuron_layer[in_first])
                            { // inhibitions within layer or across
                                P.History_time[in2].push_back(min_fire_time);
                                P.History_type[in2].push_back(2);
                                P.History_ID[in2].push_back(in_first);
                            }
                        }
                    }

                    // Create EPS signal in L0 neuron-originated streams
                    if (P.Neuron_layer[in_first] == 0)
                    { // this is a Layer-0 neuron
                        for(int in = P.N_neuronsL[0]; in < P.N_neurons; in++){
                            P.History_time[in].push_back(min_fire_time);
                            P.History_type[in].push_back(1);
                            P.History_ID[in].push_back(N_InputStreams + in_first);
                        }
                    }
                    ispike-=1;
                    insert=false;

                    //take a step back and search for another activation
                    previous_firetime = min_fire_time;
                }
            } // end if in_first fires
            else
            {
                cout << "Waiting" << endl;
            }
            
        }     // end ispike loop, ready to start over        
        ievent++; // only go to next event if we did a backward pass too
    } while (ievent <= N_events);
    
    cout << "Out";
    // dump the potentials inside a csv file
    ofstream outfile;
    outfile.open("potentials.csv");

    // Header
    outfile << "Event,Time";
    for (int in = 0; in < P_plot.N_neurons; in++)
    {
        outfile << ",V(t)_" << in;
    }
    outfile << endl;

    // content
    for (int ievent = 1; ievent <= N_events; ievent++)
    {
        for (int it = 0; it < Time[ievent-1].size(); it++)
        {
            outfile << ievent << "," << Time[ievent-1][it];

            for (int in = 0; in < P_plot.N_neurons; in++)
            {
                outfile << "," << Potential[ievent-1][in][it];
            }
            outfile << endl;
        }
    }

    // closing the input file
    delete IT;
    delete OT;
    delete dirIT;
    delete dirOT;

    file->Close();
    file_weight->Close();
    outfile.close();

    delete file_weight;
    delete file;
}

int main()
{
    // Creazione di un oggetto SNN con valori specificati solo per var1, var2 e var3
    SNN P(6,  6);
    // ReadWeights(TFile::Open("../MODE/SNNT/Histos13_NL0=6_NL1=6_NCl=6_CF01=1.00_CFI0=1.00_CFI1=1.00_alfa=0.25_0.root", "READ"), P);
    cout << "SNN initialized, let's plot the potentials" << endl;
    //command for Ema
    PlotPotentials("../MODE/SNNT/Histos13_NL0=6_NL1=6_NCl=6_CF01=1.00_CFI0=1.00_CFI1=1.00_alfa=0.25_0.root", "./ordered.root", P, 12, true, false);
    //command for Fabio
    //PlotPotentials("../MODE/SNNT/Histos13_NL0=6_NL1=6_NCl=6_CF01=1.00_CFI0=1.00_CFI1=1.00_alfa=2.00_0.root", "/Users/Fabio/Desktop/CMS-SpikingNeuralNetwork/Code/6ev_6cl_100bkg.root", P, 7);

    return 0;
}
