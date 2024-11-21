import numpy as np
import random
import operator
import sys
import os
import subprocess
import re
import pygad
import time
import json
from scipy.optimize import minimize_scalar

# Define the spike potential function in Python
def spike_potential(delta_t, K1, K2, tau_m, tau_s):
    if delta_t < 0:
        return float('inf')  # Handle out-of-range cases
    return K1 * np.exp(-delta_t / tau_m) - K2 * (np.exp(-delta_t / tau_m) - np.exp(-delta_t / tau_s))

def convert_taus(A, t_max):
    smaller_tau = A / (A + 1) / np.log(A + 1) * t_max
    greater_tau = (A + 1) * smaller_tau
    return greater_tau, smaller_tau

# Paths
PATH_JSON = '/home/ubuntu/CMS-SpikingNeuralNetwork/Code/MODE/JSON/'

# Default values
NL0, NL1 = 10, 10
#tau_m = 9.e-10
#tau_s = 1.5e-10
K1, K2 = 3.45, 5
CFI0, CF01, CFI1, K, tau_plus, tau_minus, a_plus, a_minus = 0, 0, 0, 0, 0, 0, 0, 0
exclude_L0, exclude_L1 = 0,0
# Genetic algorithm parameters
N_ev = 30000
num_generations = 100
num_parents_mating = 32
sol_per_pop = 64

# Function to read JSON files
def read_json_file(file_path):
    try:
        with open(file_path, 'r') as file:
            return json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        return None

# Function to run SNN with the given parameters and file ID
def run_SNN(N_ev, NL0, NL1, tau_m, tau_s, tau_r, tau_plus, tau_minus, a_plus, a_minus, CFI0, CF01, CFI1, alpha, TH0, TH1, K, K1, K2, IPSP_dt_dilation, d_plus, d_minus, taud_plus, taud_minus, taud_plus_2, taud_minus_2, exclude_L0, exclude_L1, file_id_GS):
    try:
        file_id_GS = str(file_id_GS).zfill(5)  # Ensure file ID is 5 digits
        command = (
            f'../SNNT13.out --NL0 {NL0} --NL1 {NL1} --N_ev {N_ev} --tau_m {tau_m} --tau_s {tau_s} --tau_r {tau_r} '
            f'--TH0 {TH0} --TH1 {TH1} --file_id_GS {file_id_GS} --Train_fraction 0.65 --d_plus {d_plus} --d_minus {d_minus} --taud_plus {taud_plus} --taud_minus {taud_minus} '
            f'--taud_plus_2 {taud_plus_2} --taud_minus_2 {taud_minus_2} --a_plus 0. --a_minus 0. '
            f' --K1 {K1} --K2 {K2} --IPSP_dt_dilation {IPSP_dt_dilation} --alpha {alpha} '
            f'--exclude_L0 {exclude_L0} --exclude_L1 {exclude_L1}'

        )
        print('Command:', command)

        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        try:
            process.wait(timeout=20 * 60)  # Timeout after 45 minutes
        except subprocess.TimeoutExpired:
            print(f"Process {file_id_GS} exceeded the time limit and was killed.")
            process.kill()
            return {'Eff': 0, 'Fr': 1, 'Q': 0, 'Selectivity': 0}

        data_js = read_json_file(PATH_JSON + f'Parameters_{file_id_GS}.json')
        if data_js is None:
            raise ValueError(f"Failed to read JSON file for file_id_GS {file_id_GS}")

        return {
            'Eff': data_js['Efficiency'],
            'Fr': data_js['Fake_rate'],
            'Q': data_js['Q'],
            'Selectivity': data_js['Selectivity']
        }
    except Exception as e:
        print(f"Error running SNN: {e}")
        return None

# Function to generate a solution index
def generate_solution_idx(population_size, generation_number, individual_idx):
    generation_digits = str(generation_number).zfill(3)
    individual_digits = str(individual_idx).zfill(2)
    return int('9' + generation_digits + individual_digits)

# Fitness function (updated to handle dynamic parameters)
def fitness_func(ga_instance, solution, solution_idx):
    global NL0, NL1, tau_m, tau_s
    current_generation = ga_instance.generations_completed
    current_id_GS = solution_idx

    # Extract solution parameters based on dynamic parameter names
    solution_dict = {parameter_names[i]: solution[i] for i in range(len(solution))}
    
    # Extract individual parameters
    TH0 = solution_dict.get('TH0', 0)
    TH1 = solution_dict.get('TH1', 0)

    d_plus = solution_dict.get('d_plus', 0)
    d_minus = solution_dict.get('d_minus', 0)
    A_plus = solution_dict.get('A_plus', 0)
    A_minus = solution_dict.get('A_minus', 0)
    t_max_plus = solution_dict.get('t_max_plus', 0)
    t_max_minus = solution_dict.get('t_max_minus', 0)
    A_EPSP = solution_dict.get('A_EPSP', 0)
    t_max_EPSP = solution_dict.get('t_max_EPSP', 0)
    alpha = solution_dict.get('alpha', 0)
    IPSP_dt_dilation = solution_dict.get('IPSP_dt_dilation', 0) 
    '''
    exclude_L0 = solution_dict.get('exclude_L0', 0)
    exclude_L1 = solution_dict.get('exclude_L1', 0)
    '''    
    

    taud_plus, taud_plus_2 = convert_taus(A_plus, t_max_plus)
    taud_minus, taud_minus_2 = convert_taus(A_minus, t_max_minus)
    tau_m, tau_s = convert_taus(A_EPSP, t_max_EPSP)
    '''
    d_plus, d_minus, A_plus, A_minus, t_max_plus, t_max_minus, alpha, IPSP_dt_dilation, exclude_L0, exclude_L1 = 0,0,0,0,0,0,0,0,0,0
    taud_plus, taud_plus_2, taud_minus, taud_minus_2 = 0, 0, 0, 0
    '''
    tau_r = minimize_scalar(spike_potential, bounds=(0, tau_m * 10), args=(K1, K2, tau_m, tau_s), method='bounded', options={'xatol': 1e-13}).x * 1.1

    output_values = run_SNN(N_ev, NL0, NL1, tau_m, tau_s, tau_r, tau_plus, tau_minus, a_plus, a_minus, CFI0, CF01, CFI1, alpha, TH0, TH1, K, K1, K2, IPSP_dt_dilation, d_plus, d_minus, taud_plus, taud_minus, taud_plus_2, taud_minus_2, exclude_L0, exclude_L1, solution_idx)

    if output_values is None:
        return [1, 1]  # Indicates failure

    selectivity = output_values['Selectivity']

    # Log solution details dynamically
    with open('values_ga_30k.csv', 'a') as file:
        solution_str = ','.join(map(str, solution))
        output_str = f"{solution_str},{output_values['Eff']},{output_values['Fr']},{selectivity},{solution_idx}\n"
        file.write(output_str)
        
    return [selectivity, 1/(output_values['Fr']+1.e-4)]
    return [selectivity]
    
    
# Fitness function wrapper
def fitness_func_wrapper(ga_instance, solution, solution_idx):
    generation_number = ga_instance.generations_completed
    population_size = ga_instance.population.shape[0]
    formatted_solution_idx = str(generate_solution_idx(population_size, generation_number, solution_idx)).zfill(5)
    return fitness_func(ga_instance, solution, formatted_solution_idx)

# Function to define the gene space dynamically
def define_gene_space():
    parameter_config = {
        'TH0': {'low': 1.4, 'high': 2.},
        'TH1': {'low': 0.65, 'high': 1.2},
        'd_plus': {'low': 0.5e-12, 'high': 1.e-11},
        'd_minus': {'low': 0.5e-12, 'high': 1.e-11},
        'A_plus': {'low': 0.5, 'high': 5},
        'A_minus': {'low': 1,  'high': 10},
        't_max_plus': {'low': 5e-10, 'high': 1.5e-9},
        't_max_minus': {'low': 2e-10, 'high': 7e-10},
        'alpha': {'low': 0.5, 'high': 3},
        'IPSP_dt_dilation': {'low': 0.1, 'high': 0.5},
        'A_EPSP': {'low': 0.1, 'high': 10},
        't_max_EPSP': {'low': 0.1e-10, 'high': 2.e-10},
        #'exclude_L0': {'low': 0, 'high': 7, 'step': 1},
        #'exclude_L1': {'low': 0, 'high': 7, 'step': 1}
    }
    return parameter_config

# Generate gene space and get parameter names
parameter_config = define_gene_space()
gene_space = [parameter_config[key] for key in parameter_config]
parameter_names = list(parameter_config.keys())

# Initialize CSV file with dynamic header
with open('values_ga_30k.csv', 'w') as file:
    header = ','.join(parameter_names) + ',Efficiency,FakeRate,Selectivity,#ID\n'
    file.write(header)

# GA configuration
ga_instance = pygad.GA(
    num_generations=num_generations,
    num_parents_mating=num_parents_mating,
    sol_per_pop=sol_per_pop,
    num_genes=len(gene_space),
    fitness_func=fitness_func_wrapper,
    gene_space=gene_space,
    parent_selection_type='nsga2',
    crossover_type="single_point",
    parallel_processing=["thread", 15],
    save_best_solutions=True,
    mutation_percent_genes=[25, 75],
    mutation_probability=[0.2, 0.75],
    mutation_type="adaptive",
    on_generation=lambda ga: print(
        f"Generation: {ga.generations_completed}, Population Size: {len(ga.population)}"
    )
)

# Run the genetic algorithm
ga_instance.run()

# Plot the fitness values
ga_instance.plot_fitness()

# Get the best solution
solution, solution_fitness, solution_idx = ga_instance.best_solution()
solution_idx = generate_solution_idx(population_size, ga_instance.generations_completed)

print(f"Parameters of the best solution: {solution}")
print(f"Fitness value of the best solution: {solution_fitness}")

# Save the best solution to the CSV file
with open('values_ga_30k.csv', 'a') as file:
    file.write(','.join(map(str, solution)) + ',' + ','.join(map(str, solution_fitness)) + '\n')
