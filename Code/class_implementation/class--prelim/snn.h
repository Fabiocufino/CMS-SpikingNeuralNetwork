

class SNN
{
private:
public:









static int Empty_buffer = 0;                  // (old) We process events every 300 TimeSteps, leaving time for L0 neurons to pass delayed signal to L1 ones
static double delta = 0.7;                    // max delta for 1Gev is 0.66rad
static double max_angle = 2.0 * M_PI + delta; // max angle to scan
static double frequency = 40e6;               // CMS tracker reading frequency [Hz]
static double omega = max_angle * frequency;  // reading angular velocity

static double z_range = 1200.;                 // bin z in [-z_range/2, z_range/2]          
static double max_R = 1200.;                  // bin r in [0, max_R]
static const short int N_bin_r = 50;
static const int N_bin_z = 0;
static const int N_InputStreams = (N_bin_r)*N_bin_z;
static double z_bin_length = z_range / N_bin_z;    
static double r_bin_length = max_R / N_bin_r;       

// -------------------------------------------
// Hits and signals related constants
// -------------------------------------------
struct Hit //holds info about one hit (position of a cluster)
{
    float z;
    float r, phi;
    short int id;

    Hit(float r_, float z_, float phi_, short int id_)
        : r(r_), z(z_), phi(phi_), id(id_)
    {
    }
};

static vector<Hit> hit_pos;

static long int last_row_event_IT = 1;         // last row associated to the previous event read in IT
static long int last_row_event_OT = 1;         // last row associated to the previous event read in OT
static const int MaxClasses = 20;
static int pclass;

static const short int BGR; /// = 1;
static const short int SIG; // = 2;

static long int NROOT;// = 100000;                //number of events inside the root file

// -------------------------------------------
// neural network constants
// -------------------------------------------
static const int MaxEvents;// = 10000000;
static const double largenumber;// = 999999999.;
static const double epsilon = 1.;// / largenumber;
static const int MaxNeurons ;///= 12;
static double ProbWSwitchUp;// = 0.5;
static double ProbWSwitchDown;// = 0.05;
static double MaxDelay;// = 0.1e-9;         // Determines shape of IE signal[s]
static double tau_m;// = 1e-09;             // membrane time constant[s]
static double tau_s;// = 0.25e-09;          // synapse time constant[s]
static double K1;// = 2.;                   // constants to tune post-synaptic shape
static double K2;// = 4.;                   // see above
static double IE_Pot_const;// = 2.5;        // Constant for IE modeling
static double Threshold[2];// = {15., 10.}; // Neuron threshold in arbitrary units; in paper it is 550V but they have 1000 channels, 100x ours
static double alpha; //= 0.25;              // 0.25; // factor tuning inhibition strength
static double L1inhibitfactor;// = 1.;      // multiplier for L1 inhibits
static double MaxDeltaT;// = 7. * tau_m;    // time window wherein pre-synaptic, inhibition, and post-synaptic kernels affect neuron potential
static double tau_plus;// = 1.68e-09;       // [s]
static double tau_minus;// = 3.37e-09;      // [s]
static double IPSP_dt_dilation;// = 1.;     // shape factor of exponential IPSP signal
static double a_plus;// = 0.03125;          // for model of EPSP
static double a_minus;// = 0.0265625;       // 0.85*a_plus;
static double MaxFactor;// = 0.2;           // Initial factor of excursion of parameters for optimization
static double eff_target;// = 0.9;
static double acc_target;// = 0.05;
static bool learnDelays;// = false;
static const int MaxStreams;// = MaxNeurons + N_bin_r * N_bin_z;



static int N_part; // Number of generated particles in an event
static double First_angle;
static double Eff[MaxNeurons * MaxClasses];   //   // Efficiency of each neuron to signals of different classes
static double Efftot[MaxClasses];                // Global efficiency to a different class
static double Weight[MaxNeurons][MaxStreams];    // Weight of synapse-neuron strength
static bool check_LTD[MaxNeurons][MaxStreams];   // checks to generate LTD after neuron discharge
static bool Void_weight[MaxNeurons][MaxStreams]; // These may be used to model disconnections
static bool bestVoid_weight[MaxNeurons][MaxStreams];
static double Weight_initial[MaxNeurons][MaxStreams]; // store to be able to return to initial conditions when optimizing
static double OldWeight[MaxNeurons][MaxStreams];      //for renorm
static double DeltaWeight[MaxNeurons][MaxStreams];
static double Delay[MaxNeurons][MaxStreams];          // Delay in incoming signals
static double bestDelay[MaxNeurons][MaxStreams];      // Opt delay in incoming signals
static vector<double> History_time[MaxNeurons];       // Time of signal events per each neuron
static vector<int> History_type[MaxNeurons];          // Type of signal
static vector<int> History_ID[MaxNeurons];            // ID of generating signal stream or neuron
static vector<double> Fire_time[MaxNeurons];          // Times of firing of each neuron
static int Neuron_layer[MaxNeurons];
double sumweight[MaxNeurons]={0};                      //summed weights of streams for each neurons for the purpose of normalization     
double SumofSquaresofWeight[MaxNeurons]={0};                    //sum of squares synaptic weights for each neuron for RMS calc
double MeanofSquaresofWeight[MaxNeurons]={0};                  //mean of squares of synaptic weights for each neuron for RMS calc
double MaxWeight[MaxNeurons];
double MinWeight[MaxNeurons];
double RMSWeight[MaxNeurons];
//double checksumweight[MaxNeurons]={0};
static int N_neuronsL[2]; // Number of neurons in layers 0 and 1
static int N_streams;
static int N_neurons;
static int N_classes;
static int N_events;
static int N_epochs;
static int NevPerEpoch;
static double ConnectedFraction_Input_L0;
static double ConnectedFraction_Input_L1;
static double ConnectedFraction_L0_L1;
static vector<double> PreSpike_Time;
static vector<int> PreSpike_Stream;
static vector<int> PreSpike_Signal; // 0 for background hit, 1 for signal hit, 2 for L1 neuron spike
static vector<int> neurons_index;   // contains neurons identifiers in random positions
static double tmax;                 // t of max value for EPSP
static double Pmax_EPSP;            // maximum EPSP spike height
static double K;                    // constant computed such that it sets the max of excitation spike at 1V
static bool update9;                // controls whether to optimize 7 network parameters
static bool updateDelays;           // controls whether to optimize neuron delays
static bool updateConnections;      // controls whether to optimize connections between streams and neurons
static bool anyHits = true;         // Whether to accept tracks with any number of hits <8 or not
static double Q_best;
static double SelL1_best;
static double Eff_best;
static double Acc_best;
static double T0_best;
static double T1_best;
static double A_best;
static double L1if_best;
static double K_best;
static double K1_best;
static double K2_best;
static double IEPC_best;
static double IPSPdf_best;
static int indfile;
static char progress[53] = "[--10%--20%--30%--40%--50%--60%--70%--80%--90%-100%]"; // Progress bar
static long int ievent;

    SNN(int NL0, int NL1);
    ~SNN();
    
//void Write_Parameters();
void Init_weights();
void Init_neurons();
 void Renorm(int in);
 void Reset_weights();
 void Init_delays();
 void Init_connection_map();
 //void Reset_hits();
 //void Encode(double t_in);
 double EPS_potential(double delta_t);
 double Spike_potential(double delta_t, int ilayer);
 double IE_potential(double delta_t, int in, int is);
 double Inhibitory_potential(double delta_t, int ilayer);
 // ---------------------------------------------------------------------------------------
void LTP(int in, int is, int this_spike, double fire_time);
void LTD(int in, int is, double spike_time);
double Neuron_firetime(int in, double t);
double LR_Scheduler(double LR0, int epoch, int Nepochs);
double Compute_Selectivity(int level, int mode);
double Compute_Q(double eff, double acc, double sel);
//void ReadFromProcessed(TTree *IT, TTree *OT, long int id_event_value);

 























}
