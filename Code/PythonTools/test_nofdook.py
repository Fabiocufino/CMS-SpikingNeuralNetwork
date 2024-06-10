# Import necessary libraries
import pandas as pd
import numpy as np

# Step 1: Organize Data
data = {
    'Neuron': [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
    'Class': [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
    'Fires': [1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1]
}

df = pd.DataFrame(data)

# Step 2: Calculate Probabilities
def calculate_probabilities(df, neuron, class_label):
    p_x1_y1 = len(df[(df['Neuron'] == neuron) & (df['Class'] == class_label) & (df['Fires'] == 1)]) / len(df)
    p_x1_y0 = len(df[(df['Neuron'] == neuron) & (df['Class'] != class_label) & (df['Fires'] == 1)]) / len(df)
    p_x0_y1 = len(df[(df['Neuron'] == neuron) & (df['Class'] == class_label) & (df['Fires'] == 0)]) / len(df)
    p_x0_y0 = len(df[(df['Neuron'] == neuron) & (df['Class'] != class_label) & (df['Fires'] == 0)]) / len(df)
    
    p_x1 = p_x1_y1 + p_x1_y0
    p_x0 = p_x0_y1 + p_x0_y0
    p_y1 = p_x1_y1 + p_x0_y1
    p_y0 = p_x1_y0 + p_x0_y0
    
    return p_x1_y1, p_x1_y0, p_x0_y1, p_x0_y0, p_x1, p_x0, p_y1, p_y0

# Step 3: Compute Metrics
def mutual_information(p_x1_y1, p_x1_y0, p_x0_y1, p_x0_y0, p_x1, p_x0, p_y1, p_y0):
    mi = 0
    if p_x1_y1 > 0:
        mi += p_x1_y1 * np.log2(p_x1_y1 / (p_x1 * p_y1))
    if p_x1_y0 > 0:
        mi += p_x1_y0 * np.log2(p_x1_y0 / (p_x1 * p_y0))
    if p_x0_y1 > 0:
        mi += p_x0_y1 * np.log2(p_x0_y1 / (p_x0 * p_y1))
    if p_x0_y0 > 0:
        mi += p_x0_y0 * np.log2(p_x0_y0 / (p_x0 * p_y0))
    return mi

def selectivity_index(p_x1_y1, p_x1_y0):
    if p_x1_y1 + p_x1_y0 > 0:
        si = (p_x1_y1 - p_x1_y0) / (p_x1_y1 + p_x1_y0)
    else:
        si = 0
    return si

def sensitivity(p_x1_y1, p_y1):
    return p_x1_y1 / p_y1 if p_y1 > 0 else 0

def specificity(p_x0_y0, p_y0):
    return p_x0_y0 / p_y0 if p_y0 > 0 else 0

# Aggregate measures for all neurons and classes
neurons = df['Neuron'].unique()
classes = df['Class'].unique()

results = []

for neuron in neurons:
    for class_label in classes:
        probabilities = calculate_probabilities(df, neuron, class_label)
        p_x1_y1, p_x1_y0, p_x0_y1, p_x0_y0, p_x1, p_x0, p_y1, p_y0 = probabilities
        
        mi = mutual_information(p_x1_y1, p_x1_y0, p_x0_y1, p_x0_y0, p_x1, p_x0, p_y1, p_y0)
        si = selectivity_index(p_x1_y1, p_x1_y0)
        sen = sensitivity(p_x1_y1, p_y1)
        spec = specificity(p_x0_y0, p_y0)
        
        results.append({
            'Neuron': neuron,
            'Class': class_label,
            'Mutual Information': mi,
            'Selectivity Index': si,
            'Sensitivity': sen,
            'Specificity': spec
        })

results_df = pd.DataFrame(results)
#import ace_tools as tools; tools.display_dataframe_to_user(name="Neural Network Specialization Metrics", dataframe=results_df)
print(results_df)

