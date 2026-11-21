# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 15:34:45 2026

@author: swapn
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 13:21:51 2025

@author: P70073624
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


sns.set(style="white")


transporters = ['NaKATPase', 'NCC', 'NKCC2', 'ENaC']

folder_path = "../09_sensitivity/"

folders = [f for f in os.listdir(folder_path)
           if os.path.isdir(os.path.join(folder_path, f))]
sensitivity_df = pd.DataFrame(columns = folders, index = [ 'BO', 'TO', 'AO', 'BTAO'])

# left panel: Transporter contribution
segments =  ['PT', 'DCT', 'mTAL', 'cTAL', 'CCD']
inhib = ['inhibition: none','NaKATPase', 'NaKATPase, \n ENaC', 'NCC, NKCC,\n ENaC', 'All']

# percent_seg = [0.23,0.02,0.73,0.02]
percent_seg = [1,1,1,1]

for folder in folders:
    filename = f'../09_sensitivity/{folder}/outlets_isolated_cell/C_outlet_isolated_cell_'
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
                
    filename = f'../09_sensitivity/{folder}/outlets_isolated_cell/O_outlet_isolated_cell_'
    
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
                
    filename = f'../09_sensitivity/{folder}/outlets_isolated_cell/BTAO_outlet_isolated_cell_'
                
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
    
    
    
    filename = f'../09_sensitivity/{folder}/outlets_isolated_cell/AO_outlet_isolated_cell_'
    
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
    
    filename = f'../09_sensitivity/{folder}/outlets_isolated_cell/BO_outlet_isolated_cell_'
    
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
    
    filename = f'../09_sensitivity/{folder}/outlets_isolated_cell/TO_outlet_isolated_cell_'
    
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
    
    # plt.figure(1, figsize=(8, 6))
    # # Combine the data arrays into a single array
    y = np.array([ sum(BO_Na), sum(TO_Na), sum(OA_Na), sum(BTAO_Na)]) / sum(O_Na) * 100  
    
    sensitivity_df[folder] = y
fitted = ['NaKATPase-70 (fitted)', 'NCC-90 (fitted)', 'NKCC2-90 (fitted)', 'ENaC-90 (fitted)']
exp_data = pd.read_excel('exp_data_extraction.xlsx', index_col=0)
y_e = exp_data.loc['Mean'].iloc[1:]
y_e_std = exp_data.loc['Std'].iloc[1:]
fig, ax = plt.subplots(2,2,figsize=(17,16), sharex = True)
ax = ax.flatten()
labels = ['(a)', '(b)', '(c)', '(d)']
for t,transporter in enumerate(transporters):

    y_b = sensitivity_df['baseline']
    
    inhib_labels = ['BO', 'TO', 'AO', 'BTAO']
    cols = sensitivity_df.columns[sensitivity_df.columns.str.contains(transporter)]
    df_sub = sensitivity_df[cols]
    
    # Create DataFrame for export
    df = pd.DataFrame({'Measured': y_e, 'Baseline': y_b,
                       f'{cols[0]}':df_sub[cols[0]],f'{cols[1]}':df_sub[cols[1]]},
                       index=inhib_labels)
    
    
    # Bar parameters
    bar_width = 0.18
    x = np.arange(len(inhib_labels))
    
    colors = sns.color_palette("Greys", n_colors=4)  # Grayscale colors
    
    # Create a new figure and axis
    
    
    # Plot simulated values
    ax[t].bar(x - 1.5*bar_width, df.iloc[:,0], color='white', lw = 2, width=bar_width, label='In vitro', edgecolor='k')
    ax[t].errorbar(x - 1.5*bar_width, df.iloc[:,0], yerr=y_e_std, fmt='none', ecolor='black', capsize=4)
    
    # Plot experimental values
    ax[t].bar(x - 0.5*bar_width, df.iloc[:,1], color='w', width=bar_width, edgecolor='k', lw = 2, label='baseline', hatch = '/')
    ax[t].bar(x + 0.5*bar_width, df.iloc[:,2], color='w', width=bar_width, edgecolor='k', lw = 2, label=f'{cols[0]}', hatch ='.')
    ax[t].bar(x + 1.5*bar_width, df.iloc[:,3], color='w', width=bar_width, edgecolor='k', lw = 2, label=f'{cols[1]}', hatch = '-')
    # Axis formatting
    ax[t].set_ylim(0, 120)
    ax[t].set_xticks(x)
    ax[t].set_xticklabels(inhib_labels, fontsize=16)
    ax[t].tick_params(axis='y', labelsize=18)
    ax[t].tick_params(axis='x', labelsize=18)
    
    ax[t].legend(fontsize=14, frameon=False, ncol = 2)
    ax[t].text(
        0.05, 0.95, labels[t],                # position (x, y)
        transform=ax[t].transAxes,          # use axes coordinates
        fontsize=14,
        fontweight='bold',
        va='top', ha='left'
    )
    # ax[t].set_title(f'Sensitivity of {transporter} inhibition')
    
    # Layout and save
ax[0].set_ylabel('Na uptake (%) normalized to ouabain', fontsize=16)
ax[2].set_ylabel('Na uptake (%) normalized to ouabain', fontsize=16)
plt.tight_layout()
plt.savefig('../09_sensitivity/sensitivity.png', dpi=300)
plt.show()

