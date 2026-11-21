# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:12:09 2023

@author: P70073624
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


segments =  ['PT', 'DCT', 'mTAL', 'cTAL', 'CCD']
inhib = ['inhibition: none','NaKATPase', 'NaKATPase, \n ENaC', 'NCC, NKCC,\n ENaC', 'All']

# percent_seg = [0.23,0.02,0.73,0.02]
percent_seg = [1,1,1,1]

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
            
print('control_Na',control_Na)
            
            
control_Na[2] += control_Na[3]
control_Na[3] = control_Na[4]

control_Na = control_Na[0:4]        
control_Na = [control_Na[i]*percent_seg[i] for i in range(len(percent_seg))]
            
filename = '../03_hypotonic-Cl-Copy/outlets_isolated_cell/O_outlet_isolated_cell_'

O_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            O_Na.append(float(words_in_line[i+1][1]))
            
print('O_na_hyp',O_Na)
            
O_Na[2] += O_Na[3]
O_Na[3] = O_Na[4]

O_Na = O_Na[0:4]        
O_Na_hyp = [O_Na[i]*percent_seg[i] for i in range(len(percent_seg))]
            
filename = '../03_hypotonic-Cl-Copy/outlets_isolated_cell/BTAO_outlet_isolated_cell_'
            
BTAO_Na = []

for seg in segments:
    fname = filename + seg + '_male_hum.txt'
    
    with open(fname, 'r') as file:
            
        # Split the line into individual words
        words_in_line = [line.strip().split() for line in file]
    
    for i,line in enumerate(words_in_line):   
        if 'Na' in line:
            BTAO_Na.append(float(words_in_line[i+1][1]))   
            
print('BTAO_na_hyp',BTAO_Na)
            
BTAO_Na[2] += BTAO_Na[3]
BTAO_Na[3] = BTAO_Na[4]

BTAO_Na = BTAO_Na[0:4]        
BTAO_Na_hyp = [BTAO_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

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
O_Na_iso = [O_Na[i]*percent_seg[i] for i in range(len(percent_seg))]
            
filename = '../02_isotonic-Cl/outlets_isolated_cell/BTAO_outlet_isolated_cell_'
            
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
BTAO_Na_iso = [BTAO_Na[i]*percent_seg[i] for i in range(len(percent_seg))]

# Data for experiment
data_exp = pd.read_excel('isotonic_hypotonic_exp.xlsx', sheet_name = 'Sheet4')
mean_exp = data_exp.mean()
std_exp = data_exp.std()

fig, ax = plt.subplots(figsize = (8,6))
bar_width = 0.3
# colors = ['#777777', '#777777', '#777777', '#777777', '#777777', '#777777']  # Grayscale colors
hatches = ['+', '//', '+', '//']
y = np.array([sum(O_Na_iso)/140,  sum(BTAO_Na_iso)/140, np.nansum(O_Na_hyp)/70, sum(BTAO_Na_hyp)/70])/ (sum(O_Na_iso)/140) * 100        
print(y)
plots = plt.bar(np.arange(4)+bar_width, y, bar_width, color = 'w', ec = 'k', hatch = '//', label = 'Simulated')
plots_e = plt.bar(np.arange(4), mean_exp,  bar_width, yerr = std_exp, color = 'w', edgecolor = 'k', label = 'In vitro', capsize = 5.0)
legend1 = plt.legend(frameon=False, fontsize=20, loc='upper left')
ax.add_artist(legend1)  # Add the first legend back to the plot
plt.ylim(50, 200)
# Set y-axis ticks to 50, 100, and 150
plt.yticks([60, 100, 140, 180], fontsize = 16)
# Iterating over the bars one-by-one
plt.text(0.3, 40, 'isotonic', fontsize = 16)
plt.text(2.3, 40, 'hypotonic', fontsize = 16)
plt.xticks([0.15, 1.15, 2.15, 3.15], labels = ['O', 'BTAO', 'O', 'BTAO'])

dummy_O = plt.bar([0], [0], color='w', hatch='+', edgecolor='k')
dummy_BTAO = plt.bar([0], [0], color='w', hatch='//', edgecolor='k')

# Add a second legend for the hatch patterns
# legend_hatches = plt.legend([dummy_O, dummy_BTAO], ['O', 'BTAO'], loc='upper right', frameon=False, fontsize=20)
# ax.add_artist(legend_hatches)

# for bar in plots.patches:
   
#     # Using Matplotlib's annotate function and
#     # passing the coordinates where the annotation shall be done
#     plt.annotate(format(bar.get_height(), '.2f'),
#                     (bar.get_x() + bar.get_width() / 2,
#                     bar.get_height()), ha='center', va='bottom',
#                     size=14, xytext=(0, 5), rotation = 90,
#                     textcoords='offset points')
    
# for bar in plots_e.patches:
#     plt.annotate(format(bar.get_height(), '.2f'),
#                     (bar.get_x() + bar.get_width() / 2,
#                     bar.get_height()), ha='center', va='bottom',
#                     size=14, xytext=(0, 5), rotation = 90,
#                     textcoords='offset points')
plt.ylabel('Na uptake (%) \n normalised to  ouabain (isotonic)', fontsize = 16)



ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('../03_hypotonic-Cl-Copy/Na-isotonic-hypotonic.png', dpi=300)



###############################################
#NEW PLOTS
segments =  ['PT', 'DCT', 'mTAL', 'cTAL', 'CCD']
percent_seg = [1,1,1,1]
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

y2 = np.array([sum(O_Na),  sum(BTAO_Na)])/ sum(O_Na) * 100        


data_e = pd.read_excel('Tub4p6_extracted_data.xlsx')
# Calculate mean
mean_e = data_e.mean()

# Calculate sample standard deviation
std_e = data_e.std(ddof=1)

fig, (axA, axB) = plt.subplots(
    1, 2,
    figsize=(11, 6),
    gridspec_kw={'width_ratios': [1, 3]}
)

bar_width = 0.3

# ============================================================
# PANEL A
# ============================================================

x = 0

# Experimental
axA.bar(
    x - bar_width/2,
    mean_e.iloc[4],
    bar_width,
    yerr=std_e.iloc[4],
    color='w',
    edgecolor='k',
    capsize=5,
    label='In vitro'
)

# Simulated
axA.bar(
    x + bar_width/2,
    y2[1],
    bar_width,
    color='w',
    edgecolor='k',
    hatch='//',
    label='Simulated'
)

axA.set_xticks([x])
axA.set_xticklabels(['BTAO'], fontsize=16)

axA.set_ylabel(
    'Na uptake (%) \n normalised to ouabain',
    fontsize=16
)


axA.tick_params(axis='y', labelsize=16)

axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)

# Panel label
axA.text(
    -0.1, 1.05, 'a',
    transform=axA.transAxes,
    fontsize=20,
    fontweight='bold',
    va='top'
)


# ============================================================
# PANEL B — your existing isotonic/hypotonic figure
# ============================================================
# mean_exp["O-iso"] = mean_e["O"]
# mean_exp["BTAO-Iso"] = mean_e["BTAO"]

y = np.array([
    sum(O_Na_iso)/140,
    sum(BTAO_Na_iso)/140,
    np.nansum(O_Na_hyp)/70,
    sum(BTAO_Na_hyp)/70
]) / (sum(O_Na_iso)/140) * 100

plots = axB.bar(
    np.arange(4) + bar_width,
    y,
    bar_width,
    color='w',
    ec='k',
    hatch='//',
    label='Simulated'
)

plots_e = axB.bar(
    np.arange(4),
    mean_exp,
    bar_width,
    yerr=std_exp,
    color='w',
    edgecolor='k',
    label='In vitro',
    capsize=5
)



data = pd.DataFrame(
    [
        y,
        [f"{m:.2f} ± {s:.2f}" for m, s in zip(mean_exp, std_exp)]
    ],
    columns=mean_exp.index,
    index=['Simulated', 'Experiment']
)

data.to_excel('Fig4b_data.xlsx', sheet_name = 'uncorrected')

axB.set_ylim(50, 200)
axB.set_yticks([60, 100, 140, 180])

axB.set_xticks([0.15, 1.15, 2.15, 3.15])
axB.set_xticklabels(['O', 'BTAO', 'O', 'BTAO'], fontsize=16)

axB.tick_params(axis='y', labelsize=16)

axB.text(0.3, 40, 'isotonic', fontsize=16)
axB.text(2.3, 40, 'hypotonic', fontsize=16)

axB.set_ylabel(
    'Na uptake (%) \n normalised to ouabain (isotonic)',
    fontsize=16
)

axB.legend(frameon=False, fontsize=16, loc='upper left')

axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)

axB.text(
    -0.07, 1.05, 'b',
    transform=axB.transAxes,
    fontsize=20,
    fontweight='bold',
    va='top'
)

for ax in [axA, axB]:
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
plt.tight_layout()
plt.savefig(
    '../03_hypotonic-Cl-Copy/Na-isotonic-hypotonic_subplots_uncorrected.png',
    dpi=300,
    bbox_inches='tight'
)


# for i in range(1,5):
#     r1 = np.arange(len(y))
#     r3 = [x - bar_width for x in r1]
#     fig, axes = plt.subplots(1,1)
#     # Plot experimental data in all subplots
#     axes.bar(r3, bar_data[:, 1], color='r', width=bar_width, label='Experiment')
#     axes.bar(r1[:2], y[:2], color ='#777777', width = bar_width, label = 'Simulated' )
    
#     # Incrementally add simulated data bars
#     if i > 2:
#         for j in range(i):
#             axes.bar(j, y[j], color=colors[j], width=bar_width)
    
#     axes.set_ylim(50, 180)
#     axes.set_yticks([60, 100, 140, 180])
#     axes.set_xticks([r - bar_width / 2 for r in range(len(y))])
#     axes.set_xticklabels(['O', 'BTAO', 'O', 'BTAO'])
#     axes.legend(fontsize=12)
    
#     axes.spines['top'].set_visible(False)
#     axes.spines['right'].set_visible(False)

#     # Add labels to the bars in each plot
#     for bar in axes.patches:
#         height = bar.get_height()
#         axes.annotate(f'{height:.2f}', (bar.get_x() + bar.get_width() / 2, height),
#                          ha='center', va='bottom', fontsize=10, rotation=90,
#                          textcoords="offset points", xytext=(0, 5))
#     plt.text(0.1, 38.5, 'isotonic', fontsize = 16)
#     plt.text(2.1, 38.5, 'hypotonic', fontsize = 16)
#     plt.tight_layout()
#     plt.show()

# plt.figure(7)

# # Number of bars
# num_segments = len(control_Na)

# # Set the positions for the bars
# ind = np.arange(num_segments)

# # Width of each bar
# width = 0.15

# # Create the grouped bar plot
# plt.bar(ind - 2 * width, control_Na, width, label='control_Na')
# plt.bar(ind - width, O_Na_iso, width, label='O_Na_iso')
# plt.bar(ind , O_Na_hyp, width, label='O_Na_hyp')
# plt.bar(ind + width, BTAO_Na_iso, width, label='BTAO_Na_iso')
# plt.bar(ind + 2 * width, BTAO_Na_hyp, width, label='BTAO_Na_hyp')


# # Set the x-axis labels
# plt.xticks(ind, ['PT', 'DCT', 'TAL', 'CNT', 'CCD'])

# # Set labels and title
# plt.xlabel('Segments')
# plt.ylabel('Na Values')
# plt.title('Na Values for Different Segments')
# plt.legend()

# # Show the plot
# plt.show()








