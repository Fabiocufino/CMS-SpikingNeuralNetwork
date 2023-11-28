# %%
#read potential.csv

#%matplotlib ipympl
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv

# read data from csv file
data = pd.read_csv('potentials.csv')
data.describe()


# %%
data['Time'][data['Event'] == 5].shape

# %%
#fig, ax = plt.subplots(6, 2)
#fig size
#fig.set_size_inches(11, 9)

# plot potential vs. time
for EV in range(1, 13):
    #plot in a matrix
    plt.subplot(6, 2, EV)
    #plt.scatter(data['Time'][data['Event'] == EV], data['V(t)_6'][data['Event'] == EV], c='b', marker='.',s=2)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_6'][data['Event'] == EV], c='red', marker='None', alpha=0.5)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_7'][data['Event'] == EV], c='orange', marker='None', alpha=0.5)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_8'][data['Event'] == EV], c='yellow', marker='None', alpha=0.5)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_9'][data['Event'] == EV], c='green', marker='None', alpha=0.5)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_10'][data['Event'] == EV], c='blue', marker='None', alpha=0.5)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_11'][data['Event'] == EV], c='purple', marker='None', alpha=0.5)

    #if EV == 1:
    #    plt.xlim(0.48e-8, 0.59e-8)

    #plot trashhold
    plt.axhline(y=0.1, color='r', linestyle='-')
    plt.xlabel('Time (s)')
    plt.ylabel('Potential (V)')
    plt.title('Event ' + str(EV))
plt.tight_layout()
#plt.savefig("V6.pdf")
plt.show()

# %%
fig, ax = plt.subplots(6, 2)
#fig size
#fig.set_size_inches(12, 11)

# plot potential vs. time
for EV in range(1, 13):
    #plot in a matrix
    plt.subplot(6, 2, EV)
    #plt.scatter(data['Time'][data['Event'] == EV], data['V(t)_6'][data['Event'] == EV], c='b', marker='.',s=2)

    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_0'][data['Event'] == EV], c='purple', marker='None', alpha=0.5)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_1'][data['Event'] == EV], c='red', marker='None', alpha=0.5)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_2'][data['Event'] == EV], c='orange', marker='None', alpha=0.5)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_3'][data['Event'] == EV], c='black', marker='None', alpha=0.5)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_4'][data['Event'] == EV], c='green', marker='None', alpha=0.5)
    plt.plot(data['Time'][data['Event'] == EV], data['V(t)_5'][data['Event'] == EV], c='blue', marker='None', alpha=0.5)

    # if EV == 4:
    #     plt.xlim(0.97e-7, 1e-7)

    #plot threshold
    plt.axhline(y=0.1, color='r', linestyle='-')
    plt.xlabel('Time (s)')
    plt.ylabel('Potential (V)')
    plt.title('Event ' + str(EV))
plt.tight_layout()
#plt.savefig("V6.pdf")
plt.show()

# %%


# %%



