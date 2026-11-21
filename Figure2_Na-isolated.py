# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 18:14:23 2023

@author: P70073624
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


segments =  ['PT', 'DCT', 'mTAL', 'cTAL', 'CCD']
inhib = ['inhibition: none','NaKATPase', 'NaKATPase, \n ENaC', 'NCC, NKCC,\n ENaC', 'All']

# percent_seg = [0.23,0.02,0.73,0.02]
percent_seg = [1,1,1,1]

filename = '../01_baseline_uncalibrated/outlets_isolated_cell/C_outlet_isolated_cell_'
control_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            control_Na.append(float(words_in_line[i+1][1]))
            
print('c_na',control_Na)
            
            
control_Na[2] += control_Na[3]
control_Na[3] = control_Na[4]

control_Na = control_Na[0:4]        
control_Na = [control_Na[i]*percent_seg[i] for i in range(len(percent_seg))]
            
filename = '../01_baseline_uncalibrated/outlets_isolated_cell/O_outlet_isolated_cell_'

O_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            O_Na.append(float(words_in_line[i+1][1]))
            
print('O_na',O_Na)
            
O_Na[2] += O_Na[3]
O_Na[3] = O_Na[4]

O_Na = O_Na[0:4]        
O_Na = [O_Na[i]*percent_seg[i] for i in range(len(percent_seg))]
            
filename = '../01_baseline_uncalibrated/outlets_isolated_cell/BTAO_outlet_isolated_cell_'
            
BTAO_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            BTAO_Na.append(float(words_in_line[i+1][1]))   
            
print('BTAO_na',BTAO_Na)
            
BTAO_Na[2] += BTAO_Na[3]
BTAO_Na[3] = BTAO_Na[4]

BTAO_Na = BTAO_Na[0:4]        
BTAO_Na = [BTAO_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

# filename = '../01_baseline_uncalibrated/outlets_isolated_cell/BTA_outlet_isolated_cell_'

# BTA_Na = []

# for seg in segments:
#     fname = filename + seg + '_male_hum.txt'
    
#     with open(fname, 'r') as file:
            
#         # Split the line into individual words
#         words_in_line = [line.strip().split() for line in file]
    
#     for i,line in enumerate(words_in_line):   
#         if 'Na' in line:
#             BTA_Na.append(float(words_in_line[i+1][1]))   
            
# BTA_Na[2] += BTA_Na[3]
# BTA_Na[3] = BTA_Na[4]
# BTA_Na = BTA_Na[0:4]        
# BTA_Na = [BTA_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

filename = '../01_baseline_uncalibrated/outlets_isolated_cell/AO_outlet_isolated_cell_'

OA_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            OA_Na.append(float(words_in_line[i+1][1]))   
            
print('OA_na',OA_Na)
OA_Na[2] += OA_Na[3]
OA_Na[3] = OA_Na[4]

OA_Na = OA_Na[0:4]        
OA_Na = [OA_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

filename = '../01_baseline_uncalibrated/outlets_isolated_cell/BO_outlet_isolated_cell_'

BO_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            BO_Na.append(float(words_in_line[i+1][1]))   
            

print('BO_na',BO_Na)
BO_Na[2] += BO_Na[3]
BO_Na[3] = BO_Na[4]

BO_Na = BO_Na[0:4]        
BO_Na = [BO_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

filename = '../01_baseline_uncalibrated/outlets_isolated_cell/TO_outlet_isolated_cell_'

TO_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            TO_Na.append(float(words_in_line[i+1][1]))   
print('TO_na',TO_Na)            
TO_Na[2] += TO_Na[3]
TO_Na[3] = TO_Na[4]

TO_Na = TO_Na[0:4]        
TO_Na = [TO_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

y1 = np.array([sum(O_Na),  sum(BO_Na), sum(TO_Na), sum(OA_Na), sum(BTAO_Na)])/ sum(O_Na) * 100        


filename = '../01_baseline/outlets_isolated_cell/C_outlet_isolated_cell_'
control_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            control_Na.append(float(words_in_line[i+1][1]))
            
print('c_na',control_Na)
            
            
control_Na[2] += control_Na[3]
control_Na[3] = control_Na[4]

control_Na = control_Na[0:4]        
control_Na = [control_Na[i]*percent_seg[i] for i in range(len(percent_seg))]
            
filename = '../01_baseline/outlets_isolated_cell/O_outlet_isolated_cell_'

O_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            O_Na.append(float(words_in_line[i+1][1]))
            
print('O_na',O_Na)
            
O_Na[2] += O_Na[3]
O_Na[3] = O_Na[4]

O_Na = O_Na[0:4]        
O_Na = [O_Na[i]*percent_seg[i] for i in range(len(percent_seg))]
            
filename = '../01_baseline/outlets_isolated_cell/BTAO_outlet_isolated_cell_'
            
BTAO_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            BTAO_Na.append(float(words_in_line[i+1][1]))   
            
print('BTAO_na',BTAO_Na)
            
BTAO_Na[2] += BTAO_Na[3]
BTAO_Na[3] = BTAO_Na[4]

BTAO_Na = BTAO_Na[0:4]        
BTAO_Na = [BTAO_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

# filename = '../01_baseline_uncalibrated/outlets_isolated_cell/BTA_outlet_isolated_cell_'

# BTA_Na = []

# for seg in segments:
#     fname = filename + seg + '_male_hum.txt'
    
#     with open(fname, 'r') as file:
            
#         # Split the line into individual words
#         words_in_line = [line.strip().split() for line in file]
    
#     for i,line in enumerate(words_in_line):   
#         if 'Na' in line:
#             BTA_Na.append(float(words_in_line[i+1][1]))   
            
# BTA_Na[2] += BTA_Na[3]
# BTA_Na[3] = BTA_Na[4]
# BTA_Na = BTA_Na[0:4]        
# BTA_Na = [BTA_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

filename = '../01_baseline/outlets_isolated_cell/AO_outlet_isolated_cell_'

OA_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            OA_Na.append(float(words_in_line[i+1][1]))   
            
print('OA_na',OA_Na)
OA_Na[2] += OA_Na[3]
OA_Na[3] = OA_Na[4]

OA_Na = OA_Na[0:4]        
OA_Na = [OA_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

filename = '../01_baseline/outlets_isolated_cell/BO_outlet_isolated_cell_'

BO_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            BO_Na.append(float(words_in_line[i+1][1]))   
            

print('BO_na',BO_Na)
BO_Na[2] += BO_Na[3]
BO_Na[3] = BO_Na[4]

BO_Na = BO_Na[0:4]        
BO_Na = [BO_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

filename = '../01_baseline/outlets_isolated_cell/TO_outlet_isolated_cell_'

TO_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            TO_Na.append(float(words_in_line[i+1][1]))   
print('TO_na',TO_Na)            
TO_Na[2] += TO_Na[3]
TO_Na[3] = TO_Na[4]

TO_Na = TO_Na[0:4]        
TO_Na = [TO_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

y2 = np.array([sum(O_Na),  sum(BO_Na), sum(TO_Na), sum(OA_Na), sum(BTAO_Na)])/ sum(O_Na) * 100        


plt.figure(1, figsize=(10,10))
# # Combine the data arrays into a single array
data_arrays = [O_Na, BO_Na, TO_Na, OA_Na, BTAO_Na]
inhib_labels = ['O', 'BO', 'TO', 'AO', 'BTAO']
# # plot bars in stack manner
# plt.bar(inhib_labels, [arr[0] for arr in data_arrays], color='#36454F', hatch = '-')
# plt.bar(inhib_labels, [arr[1] for arr in data_arrays], bottom=[arr[0] for arr in data_arrays], color='#A9A9A9', hatch = '/')
# plt.bar(inhib_labels, [arr[2] for arr in data_arrays], bottom=[arr[0]+arr[1] for arr in data_arrays], color='#6082B6', hatch = '.')
# plt.bar(inhib_labels, [arr[3] for arr in data_arrays], bottom=[arr[0]+arr[1]+arr[2] for arr in data_arrays], color='#818589', hatch = '\\')
# plt.xlabel("")
# ax = plt.gca()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.ylabel("Contribution to Na uptake", fontsize = 14)
# plt.text(-1.2, 190, '(b)', fontsize = 22)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)

# # Legend settings
# legend_labels = ["PT", "DCT", "TAL", 'CCD']
# plt.legend(legend_labels, fontsize=14, loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=4, frameon = False)
# plt.show()
# plt.savefig('../01_baseline_uncalibrated/percent-na-uptake-relative-ouabain-per-segment.png', dpi=300)
# plt.figure(2, figsize=(8, 6))
colors = ['#777777', '#777777', '#777777', '#777777', '#777777', '#777777']  # Grayscale colors
# print(y*sum(O_Na)) 

data_exp = pd.read_excel('Tub4p6_extracted_data.xlsx')
# Calculate mean
y_e = data_exp.mean()

# Calculate sample standard deviation
ye_std = data_exp.std(ddof=1)

f2data = pd.DataFrame(
    [
        [f"{mean:.2f} ± {std:.2f}" for mean, std in zip(y_e.values, ye_std.values)],
        y1,
        y2
    ],
    columns=y_e.index,
    index=['Experimental', 'Uncalibrated', 'Calibrated']
)
f2data.to_excel('figure2_data.xlsx')
# y_e = np.array([np.nan, 100, 91.1, 87.55, 75.27, 76.49])
bar_width = 0.2

ind = np.arange(len(y1[:4]))
r1 = [x - bar_width for x in ind]
r2 = [x for x in ind]
r3 = [x + bar_width for x in ind]
plots = plt.bar(r2, y1[:4], color = 'w',edgecolor = 'k', width = bar_width, label = 'Simulation-uncalibrated', hatch = '-')
plots = plt.bar(r3, y2[:4], color = 'w',edgecolor = 'k', width = bar_width, label = 'Simulation-calibrated', hatch = '/')

plots_e = plt.bar(r1, y_e[:4], yerr = ye_std[1:5], color = 'w',edgecolor = 'k', width = bar_width, label = 'Experiment' , capsize = 5)
plt.ylim(60, 120)
# Set y-axis ticks to 50, 100, and 150
plt.yticks([60, 80, 100, 120], fontsize = 24)
plt.xticks([r  for r in range(len(r1))], [ 'O', 'BO', 'TO', 'AO'], fontsize = 24)
plt.legend(fontsize = 20, frameon = False)
# plt.text(-1.3, 150, '(a)', fontsize = 22)
# for bar in plots.patches:
   
#     # Using Matplotlib's annotate function and
#     # passing the coordinates where the annotation shall be done
#     plt.annotate(format(bar.get_height(), '.2f'),
#                     (bar.get_x() + bar.get_width() / 2,
#                     bar.get_height()), ha='center', va='bottom',
#                     size=14, xytext=(0, 5), 
#                     textcoords='offset points', rotation = 90)
    
# for bar in plots_e.patches:
   
#     # Using Matplotlib's annotate function and
#     # passing the coordinates where the annotation shall be done
#     plt.annotate(format(bar.get_height(), '.2f'),
#                     (bar.get_x() + bar.get_width() / 2,
#                     bar.get_height()), ha='center', va='bottom',
#                     size=14, xytext=(0, 30),
#                     textcoords='offset points', rotation = 90)
plt.ylabel('Na uptake (%) normalised to  ouabain', fontsize = 24)

ax = plt.gca()
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
plt.savefig('Figure_2.png', dpi=300)










# for i in range(6):
#     r1 = np.arange(len(y))
#     r3 = [x - bar_width for x in r1]
#     fig, axes = plt.subplots(1,1)
#     # Plot experimental data in all subplots
#     axes.bar(r3, y_e, color='r', width=bar_width, label='Experiment')
    
#     # Incrementally add simulated data bars
#     if i > 0:
#         for j in range(i):
#             axes.bar(j, y[j], color=colors[j], width=bar_width, label='Simulation' if j == 0 else None)
    
#     axes.set_ylim(60, 120)
#     axes.set_yticks([0, 50, 100, 150])
#     axes.set_xticks([r - bar_width / 2 for r in range(len(y))])
#     axes.set_xticklabels(['O', 'BO', 'TO', 'AO', 'BTAO'])
#     axes.legend(fontsize=12)
    
#     axes.spines['top'].set_visible(False)
#     axes.spines['right'].set_visible(False)

#     # Add labels to the bars in each plot
#     for bar in axes.patches:
#         height = bar.get_height()
#         axes.annotate(f'{height:.2f}', (bar.get_x() + bar.get_width() / 2, height),
#                          ha='center', va='bottom', fontsize=10, rotation=90,
#                          textcoords="offset points", xytext=(0, 5))

#     plt.tight_layout()
#     plt.show()