#ifndef SNNT_CONSTANTS_H
#define SNNT_CONSTANTS_H
#include<math.h>
#include<vector>

using namespace std;

// hard coded constants and data used throughout the code

// -------------------------------------------
// Time encoding constants
// -------------------------------------------

static int Empty_buffer = 0;                  // (old) We process events every 300 TimeSteps, leaving time for L0 neurons to pass delayed signal to L1 ones
static double delta = 0.7;                    // max delta for 1Gev is 0.66rad
static double max_angle = 2.0 * M_PI + delta; // max angle to scan
static double frequency = 40e6;               // CMS tracker reading frequency [Hz]
static double omega = max_angle * frequency;  // reading angular velocity

static const short int N_bin_r = 10;
static const int N_bin_z = 1;
static float Left_Layers[10]  = {25, 55, 100, 140, 210, 335, 490, 655, 835, 1050}; //mm 
static float Right_Layers[10] = {40, 75, 116, 158, 290, 400, 540, 720, 895, 1110}; //mm
static float z_range = 1200.;                 // bin z in [-z_range/2, z_range/2]          
static float max_R = 1200;                    // bin r in [0, max_R]
static short int N_TrackingLayers = 10;

static float z_bin_length = z_range / N_bin_z;    
static float r_bin_length = max_R / N_bin_r;     


// -------------------------------------------
// Hits and signals related constants
// -------------------------------------------
struct Hit //holds info about one hit (position of a cluster)
{
    float z;
    float r, phi;
    short int id;
    short int pclass;

    Hit(float r_, float z_, float phi_, short int id_, short int pclass_)
        : r(r_), z(z_), phi(phi_), id(id_), pclass(pclass_)
    {
    }
};

static vector<Hit> hit_pos;

static long int last_row_event_IT = 1;         // last row associated to the previous event read in IT
static long int last_row_event_OT = 1;         // last row associated to the previous event read in OT
static const int MaxClasses = 20;
static int pclass;

static const short int BGR = 1;
static const short int SIG = 2;

// -------------------------------------------
// tracking constants
// -------------------------------------------

static int N_events = 20000;
static int N_epochs = 1;
static int NevPerEpoch = N_events / N_epochs;
static bool batch = false;
static string rootInput = "/Code/Data/muons_100k_100br_new.root";
static int N_classes = 3;
static int N_ev_classes = 3;
static int TrainingCode = 0;
static string ReadPars = "none";
static long int NROOT = 100000;                //number of events inside the root file
static bool update9 = false;                // controls whether to optimize 7 network parameters
static bool updateDelays = false;           // controls whether to optimize neuron delays
static bool updateConnections = false;      // controls whether to optimize connections between streams and neurons

// -------------------------------------------
// class constants
// -------------------------------------------

static int _NL0 = 10;
static int _NL1 = 10;
static float _alpha = 0.5;
static float _CFI0 = 1; 
static float _CFI1 = 1; 
static float _CF01 = 1;
static float _L1inhibitfactor = 1;
static float _K = 1; 
static float _K1 = 2; 
static float _K2 = 4;
static float _IE_Pot_const = 1; 
static double _IPSP_dt_dilation = 1;

static double _tau_m = 1.e-09 / 2;
static double _tau_s =  0.25e-09 / 2;
static double _tau_r = 0.5e-09 / 2;
static double _tau_plus = 1.68e-09 /2;
static double _tau_minus = 3.37e-09 /2;
static double _MaxDelay =  1.e-09;
static double _a_plus = 0.00003125;
static double _a_minus = 0.00002656; 

static double _d_plus = 1.e-11;
static double _d_minus = 1.e-11; 
static double _taud_plus = _tau_minus;
static double _taud_minus = _tau_minus;

static float _Threshold0 = 0.45;
static float _Threshold1 = 0.45;

static float _sparsity = 2;
static bool _split_layer0 = true;

static int _N_InputStreams = N_bin_r*N_bin_z;
static string SNN_PATH = "";

// -------------------------------------------
// neural network constants
// -------------------------------------------

static const int MaxEvents = 10000000;
static const double largenumber = 999999999.;
static const double epsilon = 1. / largenumber;
static const int MaxNeurons = 100;
static float ProbWSwitchUp = 0.5;
static float ProbWSwitchDown = 0.05;
static float MaxFactor = 0.2;           // Initial factor of excursion of parameters for optimization
static float eff_target = 0.9;
static float acc_target = 0.05;
static bool learnDelays = false;
static const bool nearest_spike_approx_weights = false; // Used to turn on the nearest spike approximation inside LTD and LTP functions
static const bool nearest_spike_approx_delays = true; // Used to turn on the nearest spike approximation inside LTD and LTP functions
static int N_display = 500;
static float Train_fraction = 0.9;
static double window = _tau_m*7.;

#endif // SNNT_CONSTANTS_H
