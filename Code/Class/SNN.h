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
    double bisectionMethod(double a, double b, int in, double epsilon, std::function<float(int, double, bool)> func);

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
    double IPSP_dt_dilation;     // shape factor of exponential IPSP signal

    double MaxDelay;         // Determines shape of IE signal[s]

    double tau_m;             // membrane time constant[s]
    double tau_s;          // synapse time constant[s]
    double tau_r;          // refractory time[s]
    double tau_plus;       // [s]
    double tau_minus;      // [s]

    double a_plus;          // for model of EPSP
    double a_minus;       // 0.85*a_plus;

    int N_InputStreams;
    int N_streams;

    double fire_granularity;    //it defines how much close we will look for neuron's activation.
    float fire_precision;      // [V] it defines the precision of the neuron firetime detection.


    //Variables that depend on the upper ones
    int N_neurons;
    float Threshold[2];     // Neuron threshold in arbitrary units; in paper it is 550V but they have 1000 channels, 100x ours
    double tmax;
    double MaxDeltaT;    // time window wherein pre-synaptic, inhibition, and post-synaptic kernels affect neuron potential


    /* data */
    //------- variables ---------
    float **Weight;         // Weight of synapse-neuron strength
    bool **check_LTD;       // checks to generate LTD after neuron discharge
    bool **Void_weight;     // These may be used to model disconnections
    double **Delay;          // Delay in incoming signals
    vector<double> *History_time;       // Time of signal events per each 1neuron
    vector<int> *History_type;         // Type of signal
    vector<int> *History_ID;           // ID of generating signal stream or neuron
    vector<double> *Fire_time;          // Times of firing of each neuron
    int *Neuron_layer;
    float *sumweight; // summed weights of streams for each neurons for the purpose of normalization
    int N_neuronsL[2];           // Number of neurons in layers 0 and 1
    TRandom3 *myRNG;
    double largenumber;
    double epsilon;

    SNN(int _NL0, int _NL1,
         float _alpha,
         float _CFI0, float _CFI1, float _CF01,
         float _L1inhibitfactor,
         float _K, float _K1, float _K2,
         float _IE_Pot_const, double _IPSP_dt_dilation,
         double _MaxDelay,

         double _tau_m, double _tau_s, double _tau_r, double _tau_plus, double _tau_minus,
         double _a_plus, double _a_minus,

         int _N_InputStreams,
         float _Threshold0, float _Threshold1);
         
    ~SNN();


    //------- Functions ---------
    void Init_neurons();
    void Init_weights();
    void Set_weights();

    void Init_delays();
    //void Reset_weights(float &Weight_initial);

    void Init_connection_map();
    float EPS_potential(double delta_t);
    float Spike_potential(double delta_t, int ilayer);
    float Inhibitory_potential(double delta_t, int ilayer);
    float Neuron_firetime(int in, double t);
    float Neuron_Potential(int in, double t, bool delete_history);
    float IE_potential(double delta_t, int in, int is);
    void LTP(int in, double fire_time, bool nearest_spike_approx, SNN &old);  
    void LTD(int in, int is, double spike_time,bool nearest_spike_approx, SNN &old);
    void Renorm(int in, SNN &old);
    void Renorm_Opt(int in, float delta_weight, SNN &old);
    void PrintWeights();
    void PrintSNN();

};
#endif