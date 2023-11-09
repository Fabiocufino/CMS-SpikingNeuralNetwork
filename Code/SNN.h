#ifndef SNN_H
#define SNN_H

#include <iostream>
#include <vector>

#include "Snnt_constants.h"


class SNN
{
private:
    /* data */
public:
    SNN();
    ~SNN();

    //------- Functions ---------
    void Init_neurons();
    void Init_weights();
    void Init_connection_map();
    float EPS_potential(float delta_t);
    float Spike_potential(float delta_t, int ilayer);
    float Inhibitory_potential(float delta_t, int ilayer);
    float Neuron_firetime(int in, float t);

    //------- variables ---------
    float Weight[MaxNeurons][MaxStreams];    // Weight of synapse-neuron strength
    bool check_LTD[MaxNeurons][MaxStreams];   // checks to generate LTD after neuron discharge
    bool Void_weight[MaxNeurons][MaxStreams]; // These may be used to model disconnections
    float Weight_initial[MaxNeurons][MaxStreams]; // store to be able to return to initial conditions when optimizing
    float OldWeight[MaxNeurons][MaxStreams];      //for renorm
    float Delay[MaxNeurons][MaxStreams];          // Delay in incoming signals
    vector<float> History_time[MaxNeurons];       // Time of signal events per each neuron
    vector<int> History_type[MaxNeurons];          // Type of signal
    vector<int> History_ID[MaxNeurons];            // ID of generating signal stream or neuron
    vector<float> Fire_time[MaxNeurons];          // Times of firing of each neuron
    int Neuron_layer[MaxNeurons];
    float sumweight[MaxNeurons];                      //summed weights of streams for each neurons for the purpose of normalization     
    int N_neuronsL[2]; // Number of neurons in layers 0 and 1
    int N_streams;
    int N_neurons;
    int N_classes;
    float CFI0;
    float CFI1;
    float CF01;
};


SNN::SNN()
{
}

SNN::~SNN()
{
}

#endif