#ifndef default_parameters_H
#define default_parameters_H

int NL0 = 6;
int NL1 = 6;
float _alpha = 2;
float _CFI0 = 1; 
float _CFI1 = 1; 
float _CF01 = 1;
float _L1inhibitfactor = 1;
float _K = 1; 
float _K1 = 2; 
float _K2 = 4;
float _IE_Pot_const = 1; 
float _IPSP_dt_dilation = 1;
float _MaxDelay =  0.1e-9;

float _tau_m = 1e-09 / 2;
float _tau_s =  0.25e-09 / 2;
float _tau_r = 0.5e-09 / 2;
float _tau_plus = 1.68e-09 / 2;
float _tau_minus = 3.37e-09 / 2;
float _a_plus = 0.00003125;
float _a_minus = 0.00002656;

float _Threshold0 = 0.1;
float _Threshold1 = 0.1;

int _N_InputStreams = ;
int _N_streams = ;
#endif