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

import subprocess
import re

def run_SNN(N_ev,NL0, NL1, tau_m, tau_s, tau_r, tau_plus, tau_minus, a_plus, a_minus, CFI0, CF01, CFI1, alpha, TH0, TH1, K, K1, K2,IPSP_dt_dilation ):
    try:
        command = f'../SNNT13.out --NL0 {NL0} --NL1 {NL1} --N_ev {N_ev} --tau_m {tau_m} --tau_s {tau_s} --tau_r {tau_r} --tau_plus {tau_plus} --tau_minus {tau_minus} --a_plus {a_plus} --a_minus {a_minus} --CFI0 {CFI0} --CF01 {CF01} --CFI1 {CFI1} --alpha {alpha} --TH0 {TH0} --TH1 {TH1} --{K} --{K1} --{K2} --{IPSP_dt_dilation}'
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        

        values = {
            'Eff': 0,
            'Fr': 0,
            'Q': 0,
            'Selectivity': 0,
        }

        for line in iter(process.stdout.readline, ''):
            line = line.strip()

            print(line)
            if 'Average efficiency:' in line:
                match = re.search(r'Average efficiency: (\d+\.?\d*)', line)
                if match:
                    values['Eff'] = float(match.group(1))
                else:
                    values['Eff'] = 0
                    
            elif 'Average fake rate:' in line:
                match = re.search(r'Average fake rate: (\d+\.?\d*)', line)
                if match:
                    values['Fr'] = float(match.group(1))
            elif 'Maximum Q value:' in line:
                match = re.search(r'Maximum Q value: (\d+\.?\d*)', line)
                if match:
                    values['Q'] = float(match.group(1))
            elif 'L1 selectivity:' in line:
                match = re.search(r'L1 selectivity: (\d+\.?\d*)', line)
                if match:
                    values['Selectivity'] = float(match.group(1))
        
        process.stdout.close()
        process.wait()
        
        return values
    
    except Exception as e:
        print(f"Error during the execution of SNN: {e}")
        return None

# Fitness function for multi-objective optimization
def fitness_func(ga_instance, solution, solution_idx):
    N_ev = 1000
    NL0, NL1, tau_m, tau_s, tau_r, tau_plus, tau_minus, a_plus, a_minus, CFI0, CF01, CFI1, alpha, TH0, TH1, K, K1, K2, IPSP_dt_dilation = solution

    output_values = run_SNN(N_ev, NL0, NL1, tau_m, tau_s, tau_r, tau_plus, tau_minus, a_plus, a_minus, CFI0, CF01, CFI1, alpha, TH0, TH1, K, K1, K2, IPSP_dt_dilation)
    
    if output_values is None:
        return [1000, 1000, 1000]  # Large values to indicate failure
    
    efficiency = output_values['Eff']
    fake_rate = output_values['Fr']
    selectivity = output_values['Selectivity']

    # Minimize fake rate (thus using 1/(fake_rate + 1e-6) in fitness)
    fitness = [efficiency, 1/(fake_rate + 1e-6), selectivity]

    # Append the solution and its fitness to the CSV file, but save the real fake_rate
    with open('values_ga_30k.csv', 'a') as file:
        file.write(','.join(map(str, solution)) + ',' + str(efficiency) + ',' + str(fake_rate) + ',' + str(selectivity) + '\n')

    return fitness


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
num_parents_mating = 2
sol_per_pop = 2
num_genes = len(gene_space)


ga_instance = pygad.GA(num_generations=num_generations,
                       num_parents_mating=num_parents_mating,
                       sol_per_pop=sol_per_pop,
                       num_genes=num_genes,
                       fitness_func=fitness_func,
                       gene_space=gene_space,
                       parent_selection_type='nsga2',
                       crossover_type="single_point",
                       mutation_type="random",
                       parallel_processing=["thread", 28],
                       save_best_solutions=True,
                       on_generation=lambda ga_instance: print(f"Generation: {ga_instance.generations_completed}"))

# Run the GA
ga_instance.run()

# Plot the fitness values
ga_instance.plot_fitness()

# Get the best solution
solution, solution_fitness, solution_idx = ga_instance.best_solution()
print(f"Parameters of the best solution: {solution}")
print(f"Fitness value of the best solution: {solution_fitness}")

# Note: Saving the best solution separately is optional since all solutions are being saved during the run.
# Save the best solution to the CSV file
with open('values_ga_30k.csv', 'a') as file:
    file.write(','.join(map(str, solution)) + ',' + ','.join(map(str, solution_fitness)) + '\n')
