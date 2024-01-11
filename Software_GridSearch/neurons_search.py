import itertools

#---------------------------- STEPS ----------------------------
#create a main folder containing SNNT13.C and this file
#create a subfolder starting_parameters
#edit this file to choose the combinations you're intrested in
#launch this code
#it will generate a .sh file for each configuration in order to launch the neural network
#edit and launch run.sh script
#----------------------------------------------------------------

#------------------------- PARAMETERS --------------------------
#edit with the parameters value you're interested in
N_ev = 100000                     
N_ep = 1                        #1 to turn off the parameters optimization
# rootInput = "/lustre/cmsdata/SNN_group/100k_100br.root"   #name of the root file
# batch = 1                       #True -> Batch mode
NL0_list =          [10, 15, 20, 25, 30]       
NL1_list =          [10]
tau_m_list =        [1e-09]
tau_s_list =        [0.25e-09]
tau_plus_list =     [1.68e-09]
tau_minus_list =    [1.68e-09]
a_plus_list =       [0.03125]
a_minus_list =      [0.02656]
CFI0_list =         [1, 0.8, 0.6] 
CFI1_list =         [1]
CF01_list =         [1]
a_list =            [1, 5, 10, 50, 100]      #alfa parameter -> inhibition strength
Thresh0 =           35
Thresh1 =           35

parameters_list = (NL0_list, NL1_list, tau_m_list, tau_s_list, tau_plus_list, tau_minus_list, a_plus_list, a_minus_list, CFI0_list, CFI1_list, CF01_list, a_list)
#----------------------------------------------------------------
tot_comb = 1
for param in parameters_list:
    tot_comb*=len(param)

print(f"Total combinations: {tot_comb}")

input("Press enter to continue")

#File creation
#create an iterator with all the possible combinations
all_combinations = itertools.product(*parameters_list)

for i, combination in enumerate(all_combinations):
    NL0, NL1, tau_m, tau_s, tau_plus, tau_minus, a_plus, a_minus, CFI0, CFI1, CF01, a = combination
    file_content = f''' --N_ev {N_ev} --N_ep {N_ep} --NL0 {NL0} --NL1 {NL0} --tau_m {tau_m} --tau_s {tau_s} --tau_plus {tau_plus} --tau_minus {tau_minus} --a_plus {a_plus} --a_minus {a_minus} --CFI0 {CFI0} --CF01 {CFI0} --CFI1 {CFI0} --alpha {a} --TH0 {Thresh0} --TH1 {Thresh1} --file_id_GS {i}'''

    print("File creation")
    filename = f"starting_parameters/start_{i}.cmd"
    i+=1
    with open(filename, "w") as file:
        file.write(file_content)
    print(f"File {filename} successfully created.")