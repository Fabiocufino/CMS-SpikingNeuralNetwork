import numpy as np
import random
from random import randint
import operator
from array import *
import sys 
import os
import subprocess
import re
import pygad
import time
import json

# PATHS
PATH_JSON = '/home/centos/CMS-SpikingNeuralNetwork/Code/MODE/JSON/'

# Function to read JSON files
def read_json_file(file_path):
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
            return data
    except Exception as e:
        print(f"Error during the reading of the JSON file: {e}")
        return None

# Function to run SNN with the given parameters and file ID
def run_SNN(N_ev, NL0, NL1, tau_m, tau_s, tau_r, tau_plus, tau_minus, a_plus, a_minus, CFI0, CF01, CFI1, alpha, TH0, TH1, K, K1, K2, IPSP_dt_dilation, file_id_GS):
    try:
        command = f'../SNNT13.out --NL0 {NL0} --NL1 {NL1} --N_ev {N_ev} --tau_m {tau_m} --tau_s {tau_s} --tau_r {tau_r} --tau_plus {tau_plus} --tau_minus {tau_minus} --a_plus {a_plus} --a_minus {a_minus} --CFI0 {CFI0} --CF01 {CF01} --CFI1 {CFI1} --alpha {alpha} --TH0 {TH0} --TH1 {TH1} --K {K} --K1 {K1} --K2 {K2} --IPSP_dt_dilation {IPSP_dt_dilation} --file_id_GS {file_id_GS}'

        print(f"Command: {command}")
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        
        values = {
            'Eff': 0,
            'Fr': 0,
            'Q': 0,
            'Selectivity': 0,
        }

        #wait for the process to terminate
        process.wait()

        # Read the output of the SNN
        data_js = read_json_file(PATH_JSON + 'Parameters_' + file_id_GS + '.json')

        values['Eff'] = data_js['Efficiency']
        values['Fr'] = data_js['Fake_rate']
        values['Q'] = data_js['Q']
        values['Selectivity'] = data_js['Selectivity']

        return values
    
    except Exception as e:
        print(f"Error during the execution of SNN: {e}")
        return None

# Function to generate solution index
def generate_solution_idx(population_size, generation_number):
    """
    Generates a solution index based on the population size and generation number.

    Args:
        population_size (int): The total size of the population.
        generation_number (int): The current generation number.

    Returns:
        int: The solution index as a four-digit integer.
    """
    population_digits = str(population_size)[:2]  # Get the first two digits of the population size
    generation_digits = str(generation_number).zfill(2)  # Get the last two digits of the generation number
    solution_idx = int(population_digits + generation_digits)  # Combine and convert to integer
    return solution_idx

# Your original fitness function
def fitness_func(ga_instance, solution, solution_idx):
    # Extract the current generation and population size from the GA instance
    current_generation = ga_instance.generations_completed
    population_size = len(ga_instance.population)

    N_ev = 1000
    print("---------------------------", solution_idx)

    NL0, NL1, tau_m, tau_s, tau_r, tau_plus, tau_minus, a_plus, a_minus, CFI0, CF01, CFI1, alpha, TH0, TH1, K, K1, K2, IPSP_dt_dilation = solution

    output_values = run_SNN(N_ev, NL0, NL1, tau_m, tau_s, tau_r, tau_plus, tau_minus, a_plus, a_minus, CFI0, CF01, CFI1, alpha, TH0, TH1, K, K1, K2, IPSP_dt_dilation, str(solution_idx))
    
    if output_values is None:
        return [1000, 1000, 1000]  # Large values to indicate failure
    
    efficiency = output_values['Eff']
    fake_rate = output_values['Fr']
    selectivity = output_values['Selectivity']

    # Minimize fake rate (thus using 1/(fake_rate + 1e-6) in fitness)
    fitness = [efficiency, 1/(fake_rate + 1e-6), selectivity]

    # Append the solution and its fitness to the CSV file
    with open('values_ga_30k.csv', 'a') as file:
        file.write(','.join(map(str, solution)) + ',' + str(efficiency) + ',' + str(fake_rate) + ',' + str(selectivity) + ',' + str(current_generation) + ',' + str(population_size) + '\n')

    return fitness

# Wrapper function to ensure the correct parameters are passed
def fitness_func_wrapper(ga_instance, solution, solution_idx):
    solution_idx = generate_solution_idx(ga_instance.population.shape[0], ga_instance.generations_completed)
    return fitness_func(ga_instance, solution, solution_idx)

# Create the values.csv and write the header
with open('values_ga_30k.csv', 'w') as file:
    file.write('NL0, NL1,tau_m,tau_s,tau_r,tau_plus,tau_minus,a_plus,a_minus,CFI0,CF01,CFI1,alpha,TH0,TH1,K,K1,K2,IPSP_dt_dilation,Efficiency,FakeRate,Selectivity\n')

# Gene space --------------------------
NLO_MIN, NLO_MAX = 5, 12
NL1_MIN, NL1_MAX = 5, 12
tau_m_MIN, tau_m_MAX = 1e-10, 1e-8
tau_s_MIN, tau_s_MAX = 1e-10, 1e-8
tau_r_MIN, tau_r_MAX = 1e-10, 1e-8
tau_plus_MIN, tau_plus_MAX = 1e-10, 1e-8
tau_minus_MIN, tau_minus_MAX = 1e-10, 1e-8
a_minus_MIN, a_minus_MAX = 0.00000656, 0.00009656
a_plus_MIN, a_plus_MAX = 0.00000125, 0.00009125
CFI0_MIN, CFI0_MAX = 0.5, 0.9
CF01_MIN, CF01_MAX = 0.5, 0.9
CFI1_MIN, CFI1_MAX = 0.5, 0.9
alpha_MIN, alpha_MAX = 0.1, 1
TH0_MIN, TH0_MAX = 0.6, 0.95
TH1_MIN, TH1_MAX = 0.6, 0.95
K_MIN, K_MAX = 1, 10
K1_MIN, K1_MAX = 1, 10
K2_MIN, K2_MAX = 1, 10
IPSP_dt_dilation_MIN, IPSP_dt_dilation_MAX = 0.001, 1.0

gene_space = [
    {'low': NLO_MIN, 'high': NLO_MAX, 'step': 1},
    {'low': NL1_MIN, 'high': NL1_MAX, 'step': 1},
    {'low': tau_m_MIN, 'high': tau_m_MAX},
    {'low': tau_s_MIN, 'high': tau_s_MAX},
    {'low': tau_r_MIN, 'high': tau_r_MAX},
    {'low': tau_plus_MIN, 'high': tau_plus_MAX},
    {'low': tau_minus_MIN, 'high': tau_minus_MAX},
    {'low': a_minus_MIN, 'high': a_minus_MAX},
    {'low': a_plus_MIN, 'high': a_plus_MAX},
    {'low': CFI0_MIN, 'high': CFI0_MAX},
    {'low': CF01_MIN, 'high': CF01_MAX},
    {'low': CFI1_MIN, 'high': CFI1_MAX},
    {'low': alpha_MIN, 'high': alpha_MAX},
    {'low': TH0_MIN, 'high': TH0_MAX},
    {'low': TH1_MIN, 'high': TH1_MAX},
    {'low': K_MIN, 'high': K_MAX},
    {'low': K1_MIN, 'high': K1_MAX},
    {'low': K2_MIN, 'high': K2_MAX},
    {'low': IPSP_dt_dilation_MIN, 'high': IPSP_dt_dilation_MAX}
]
# --------------------------

## Algorithm

# GA Configuration
num_generations = 10000
num_parents_mating = 4
sol_per_pop = 4
num_genes = len(gene_space)

# Initialize the GA instance with the wrapper function
ga_instance = pygad.GA(num_generations=num_generations,
                       num_parents_mating=num_parents_mating,
                       sol_per_pop=sol_per_pop,
                       num_genes=num_genes,
                       fitness_func=fitness_func_wrapper,  # Use the wrapper
                       gene_space=gene_space,
                       parent_selection_type='nsga2',
                       crossover_type="single_point",
                       mutation_type="random",
                       parallel_processing=["thread", 28],
                       save_best_solutions=True,
                       on_generation=lambda ga_instance: print(
                           f"Generation: {ga_instance.generations_completed}, "
                           f"Population Size: {len(ga_instance.population)}"
                       ))

# Run the GA
ga_instance.run()

# Plot the fitness values
ga_instance.plot_fitness()

# Get the best solution
solution, solution_fitness, solution_idx = ga_instance.best_solution()
population_size = 10  # Example population size
solution_idx = generate_solution_idx(population_size, ga_instance.generations_completed)
print(f"Parameters of the best solution: {solution}")
print(f"Fitness value of the best solution: {solution_fitness}")

# Save the best solution to the CSV file
with open('values_ga_30k.csv', 'a') as file:
    file.write(','.join(map(str, solution)) + ',' + ','.join(map(str, solution_fitness)) + '\n')
