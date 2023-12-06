#ifndef SNN_H
#define SNN_H

#include <iostream>
#include <vector>

#include "TF1.h"
#include "TRandom3.h"

using namespace std;

class SNN
{
private:
    float bisectionMethod(float a, float b, int in, float epsilon, std::function<float(int, float, bool)> func);

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
    float tau_r;          // refractory time[s]
    float tau_plus;       // [s]
    float tau_minus;      // [s]

    float a_plus;          // for model of EPSP
    float a_minus;       // 0.85*a_plus;

    int N_InputStreams;
    int N_streams;

    float fire_granularity;    //it defines how much close we will look for neuron's activation.
    float fire_precision;      // [V] it defines the precision of the neuron firetime detection.


    //Variables that depend on the upper ones
    int N_neurons;
    float Threshold[2];     // Neuron threshold in arbitrary units; in paper it is 550V but they have 1000 channels, 100x ours
    float tmax;
    float MaxDeltaT;    // time window wherein pre-synaptic, inhibition, and post-synaptic kernels affect neuron potential


    /* data */
    //------- variables ---------
    float **Weight;         // Weight of synapse-neuron strength
    bool **check_LTD;       // checks to generate LTD after neuron discharge
    bool **Void_weight;     // These may be used to model disconnections
    float **Weight_initial; // store to be able to return to initial conditions when optimizing
    float **OldWeight;      // for renorm
    float **Delay;          // Delay in incoming signals
    vector<float> *History_time;       // Time of signal events per each 1neuron
    vector<int> *History_type;         // Type of signal
    vector<int> *History_ID;           // ID of generating signal stream or neuron
    vector<float> *Fire_time;          // Times of firing of each neuron
    int *Neuron_layer;
    float *sumweight; // summed weights of streams for each neurons for the purpose of normalization
    int N_neuronsL[2];           // Number of neurons in layers 0 and 1
    TRandom3 *myRNG;
    float largenumber;
    float epsilon;

    SNN(int _NL0, int _NL1,
         float _alpha,
         float _CFI0, float _CFI1, float _CF01,
         float _L1inhibitfactor,
         float _K, float _K1, float _K2,
         float _IE_Pot_const, float _IPSP_dt_dilation,
         float _MaxDelay,

         float _tau_m, float _tau_s, float _tau_r, float _tau_plus, float _tau_minus,
         float _a_plus, float _a_minus,

         int _N_InputStreams,
         float _Threshold0, float _Threshold1);
         
    ~SNN();


    //------- Functions ---------
    void Init_neurons();
    void Init_weights();
    void Set_weights();

    void Init_delays();
    void Reset_weights();

    void Init_connection_map();
    float EPS_potential(float delta_t);
    float Spike_potential(float delta_t, int ilayer);
    float Inhibitory_potential(float delta_t, int ilayer);
    float Neuron_firetime(int in, float t);
    float Neuron_Potential(int in, float t, bool delete_history);
    float IE_potential(float delta_t, int in, int is);

};
#endif