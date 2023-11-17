#ifndef SNN_H
#define SNN_H

#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TF1.h"

#include "Snnt_constants.h"

using namespace std;

class SNN
{
private:
public:


    float alpha;              // 0.25; // factor tuning inhibition strength

    float CFI0;
    float CFI1;
    float CF01;

    float L1inhibitfactor;      // multiplier for L1 inhibits

    float K;
    float K1;                   // constants to tune post-synaptic shape
    float K2;                   // see above

    float IE_Pot_const;        // Constant for IE modeling
    float IPSP_dt_dilation;     // shape factor of exponential IPSP signal

    float MaxDelay;         // Determines shape of IE signal[s]

    float tau_m;             // membrane time constant[s]
    float tau_s;          // synapse time constant[s]
    float tau_plus;       // [s]
    float tau_minus;      // [s]
    float a_plus;          // for model of EPSP
    float a_minus;       // 0.85*a_plus;

    int N_streams;
    int N_classes;

    float MaxFactor;           // Initial factor of excursion of parameters for optimization



    //Variables that depend on the upper ones
    int N_neurons;
    float Threshold[2];     // Neuron threshold in arbitrary units; in paper it is 550V but they have 1000 channels, 100x ours
    float tmax;
    float MaxDeltaT;    // time window wherein pre-synaptic, inhibition, and post-synaptic kernels affect neuron potential


    /* data */
    //------- variables ---------
    float Weight[MaxNeurons][MaxStreams];         // Weight of synapse-neuron strength
    bool check_LTD[MaxNeurons][MaxStreams];       // checks to generate LTD after neuron discharge
    bool Void_weight[MaxNeurons][MaxStreams];     // These may be used to model disconnections
    float Weight_initial[MaxNeurons][MaxStreams]; // store to be able to return to initial conditions when optimizing
    float OldWeight[MaxNeurons][MaxStreams];      // for renorm
    float Delay[MaxNeurons][MaxStreams];          // Delay in incoming signals
    vector<float> History_time[MaxNeurons];       // Time of signal events per each neuron
    vector<int> History_type[MaxNeurons];         // Type of signal
    vector<int> History_ID[MaxNeurons];           // ID of generating signal stream or neuron
    vector<float> Fire_time[MaxNeurons];          // Times of firing of each neuron
    int Neuron_layer[MaxNeurons];
    float sumweight[MaxNeurons]; // summed weights of streams for each neurons for the purpose of normalization
    int N_neuronsL[2];           // Number of neurons in layers 0 and 1

    SNN(int NL0, int NL1);
    ~SNN();


    //------- Functions ---------
    void Init_neurons(float t_in);
    void Init_neurons();
    void Init_weights();
    void Init_connection_map();
    float EPS_potential(float delta_t);
    float Spike_potential(float delta_t, int ilayer);
    float Inhibitory_potential(float delta_t, int ilayer);
    float Neuron_firetime(int in, float t);
    float Neuron_Potential(int in, float t);
    float IE_potential(float delta_t, int in, int is);

};
#endif