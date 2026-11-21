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
import matplotlib.gridspec as gridspec

sns.set(style="white")

transp = {
    "PT": ["NaKATPase", "NHE3"],
    "DCT": ["NCC", "NaKATPase", "ENaC", "NHE3"],
    "cTAL": ["NaKATPase", "NKCC2A", "NKCC2B", "NHE3"],
    "mTAL": ["NaKATPase", "NKCC2A", "NKCC2F", "NHE3"],
    "CNT": ["NaKATPase", "ENaC"],
    "CCD": ["NaKATPase", "ENaC"]
}

df = pd.DataFrame(columns=['C', 'O','BO', 'TO', 'AO', 'BTAO'])
for i in ['C', 'O','BO', 'TO', 'AO', 'BTAO']:
    folder = f'../01_baseline/{i}_Male_hum_isolated_cell/'
    for key, values in transp.items():
        for value in values:
            k = f'{key}-{value}'
            df.loc[k, i] = 0
            for root, dirs, files in os.walk(folder):
                for file in files:
                    if key in file and value in file and 'Na' in file.split('_')[4]:
                        with open(os.path.join(root, file), 'r') as f:
                            df.loc[k, i] += float(f.read())

df_normalized = df.copy()
mask = df['O'] != 0
df_normalized[df.columns] = df.div(np.where(mask, df['O'], 1), axis=0)
df_normalized[df.columns] = df_normalized[df.columns].replace(-0.0, 0.0)
df_hue = (df - df.min()) / (df.max() - df.min())
melted_hue = df_hue.reset_index().melt(id_vars='index', var_name='Inhibition', value_name='normalised to Na-K-ATPase')
df_normalized.drop(['DCT-NHE3'], inplace=True)
melted_df = df_normalized.reset_index().melt(id_vars='index', var_name='Inhibition', value_name='Normalized to Ouabain')
melted_df = melted_df.rename(columns={'index': 'Flux'})
melted_df['normalised to Na-K-ATPase'] = melted_hue['normalised to Na-K-ATPase']

fig = plt.figure(figsize=(8, 10))
gs = gridspec.GridSpec(2, 1)  # 1/3rd for panel A, 2/3rd for panel B

# Define axes
ax1 = fig.add_subplot(gs[0])  # Panel A (1/3rd)
ax2 = fig.add_subplot(gs[1:]) # Panel B (2/3rd)

# Mark rows where values exceed 100
melted_df['Capped'] = melted_df['Normalized to Ouabain'] > 100

# Cap values above 100 to 5
melted_df.loc[melted_df['Capped'], 'Normalized to Ouabain'] = 5

# Plot the clustered bar chart with grayscale color scheme

sns.barplot(data=melted_df, x='Flux', y='Normalized to Ouabain', hue='Inhibition', palette='Greys', ax = ax2)

# Highlight capped bars in dark red
for bar, (_, row) in zip(ax2.patches, melted_df.iterrows()):
    if row['Capped']:
        bar.set_color('darkred')
        

    bar.set_width(bar.get_width() * 0.7)  # Reduce bar width
    bar.set_edgecolor('black')  # Add black border
    bar.set_linewidth(1.5)  # Set border thickness

# Rotate x-axis labels by 45 degrees
for label in ax2.get_xticklabels():
    label.set_rotation(45)
    label.set_horizontalalignment('right')

ax2.set_ylabel('Predicted flux (mmol/s) normalized \n to ouabain', fontsize = 12)
ax2.set_xlabel('')
ax1.text(0.05, 110, '(a)', fontsize = 16)
ax2.text(-1.5, 5.5, '(b)', fontsize = 16)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
# Adjust layout
plt.tight_layout()
handles, _ = ax2.get_legend_handles_labels()
ax2.legend(handles, [ 'None', 'O','BO', 'TO', 'AO', 'BTAO'], title="Inhibition", frameon=False, ncol=3)

# Show the plot
plt.show()

# left panel: Transporter contribution
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

# plt.figure(1, figsize=(8, 6))
# # Combine the data arrays into a single array
data_arrays = [O_Na, BO_Na, TO_Na, OA_Na, BTAO_Na]
inhib_labels = ['None','O', 'BO', 'TO', 'AO', 'BTAO']



y = np.array([sum(control_Na), sum(O_Na),  sum(BO_Na), sum(TO_Na), sum(OA_Na), sum(BTAO_Na)])/ sum(O_Na) * 100    
colors = sns.color_palette("Greys", n_colors=len(y))  # Generate `n_colors` grayscale shades

y_e = np.array([np.nan, 100, 81.60, 93.61, 82.43, 72.38])

# Create a DataFrame
df = pd.DataFrame({'Measured': y_e,'Predicted': y}, index = inhib_labels)
df.to_csv('comparison_with_expt.csv')
bar_width = 0.3
r1 = np.arange(len(y))
r2 = [x - bar_width - 0.05 for x in r1]

ax1.bar(r1, y, color=colors, width=bar_width, label='Simulation', edgecolor='k', hatch = '//')
ax1.bar(r2, y_e, color='white', width=bar_width, edgecolor='k', label='Experiment')

ax1.set_ylim(0, 120)
ax1.set_xticks([r - bar_width/2 for r in range(len(y))])
ax1.set_xticklabels(['None', 'O', 'BO', 'TO', 'AO', 'BTAO'], fontsize = 16)
ax1.tick_params(axis = 'y', labelsize = 18)
ax1.tick_params(axis = 'x', labelsize = 18)
# for label in ax1.get_xticklabels():
#     label.set_rotation(45)
#     label.set_horizontalalignment('center')
ax1.legend(fontsize=10)
ax1.set_ylabel('Na uptake (%) normalized to ouabain', fontsize = 16)


plt.tight_layout()
plt.savefig('../01_baseline/two_panel_plot.png', dpi=300)
plt.show()


# Mark rows where values exceed 100
melted_df['Capped'] = melted_df['Normalized to Ouabain'] > 100

# Cap values above 100 to 5
melted_df.loc[melted_df['Capped'], 'Normalized to Ouabain'] = 5

# Create the plot
fig, ax = plt.subplots(figsize=(8, 6))

# Plot the clustered bar chart with grayscale color scheme
sns.barplot(data=melted_df, x='Flux', y='Normalized to Ouabain', hue='Inhibition', palette='Greys', ax=ax)

# Highlight capped bars in dark red
for bar, (_, row) in zip(ax.patches, melted_df.iterrows()):
    if row['Capped']:
        bar.set_color('darkred')

    bar.set_width(bar.get_width() * 0.7)  # Reduce bar width
    bar.set_edgecolor('black')           # Add black border
    bar.set_linewidth(1.5)               # Set border thickness

# Rotate x-axis labels by 45 degrees
for label in ax.get_xticklabels():
    label.set_rotation(45)
    label.set_horizontalalignment('right')

# Styling
ax.set_ylabel('Predicted flux normalized \n to ouabain', fontsize=16)
ax.set_xlabel('')
# ax.text(-1.5, 5.5, '(a)', fontsize=16)  # Single label now
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Adjust layout and legend
plt.tight_layout()
handles, _ = ax.get_legend_handles_labels()
ax.legend(handles, ['None', 'O', 'BO', 'TO', 'AO', 'BTAO'], title="Inhibition", frameon=False, ncol=3)
plt.savefig('../01_baseline/part_a.png', dpi=1200)
# Show the plot
plt.show()

y = np.array([sum(control_Na), sum(O_Na), sum(BO_Na), sum(TO_Na), sum(OA_Na), sum(BTAO_Na)]) / sum(O_Na) * 100    
y_e = np.array([np.nan, 100, 81.60, 93.61, 82.43, 72.38])
inhib_labels = ['None', 'O', 'BO', 'TO', 'AO', 'BTAO']

# Create DataFrame for export
df = pd.DataFrame({'Measured': y_e, 'Predicted': y}, index=inhib_labels)
df.to_csv('comparison_with_expt.csv')
#%%
# Bar parameters
bar_width = 0.3
r1 = np.arange(len(y))
r2 = r1 - bar_width - 0.05
colors = sns.color_palette("Greys", n_colors=len(y))  # Grayscale colors

# Create a new figure and axis
fig, ax = plt.subplots(figsize=(8, 6))

# Plot simulated values
ax.bar(r1, y, color=colors, width=bar_width, label='Simulation', edgecolor='k', hatch='//')

# Plot experimental values
ax.bar(r2, y_e, color='white', width=bar_width, edgecolor='k', label='Experiment')

# Axis formatting
ax.set_ylim(0, 120)
ax.set_xticks(r1 - bar_width / 2)
ax.set_xticklabels(inhib_labels, fontsize=16)
ax.tick_params(axis='y', labelsize=18)
ax.tick_params(axis='x', labelsize=18)
ax.set_ylabel('Na uptake (%) normalized to ouabain', fontsize=16)
ax.legend(fontsize=10, frameon=False)

# Layout and save
plt.tight_layout()
plt.savefig('../01_baseline/part_b.png', dpi=1200)
plt.show()

# plot in multiple subplots
#%%

df = melted_df[melted_df["Inhibition"] != "C"].copy()
flux_order = df["Flux"].unique()
conditions = ["AO", "BO", "TO", "BTAO"]
fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=True)
axes = axes.flatten()
letters = ['(a)', '(b)', '(c)', '(d)']

# Compute O reference once
O_values = df[df["Inhibition"] == "O"]["Normalized to Ouabain"]
O_mean = O_values.mean()

for i, (ax, cond) in enumerate(zip(axes, conditions)):

    # Subset ONLY this condition (exclude O)
    sub = df[df["Inhibition"] == cond].copy()

    # Plot bars
    sns.barplot(
        data=sub,
        x="Flux",
        y="Normalized to Ouabain",
        order=flux_order,
        color="white",
        ax=ax,
        width=0.7,
        errorbar=None
    )

    # Style bars
    for bar, (_, row) in zip(ax.patches, sub.iterrows()):
        if row.get("Capped", False):
            bar.set_color("darkred")
        bar.set_edgecolor("black")
        bar.set_linewidth(1.5)
        # bar.set_hatch('/')   # 👈 optional: keep hatch if you like

    #Add O horizontal reference
    ax.axhline(
        O_mean,
        color='black',
        linestyle='--',
        linewidth=2
    )

    # Labels
    ax.set_title(cond, fontsize=14)
    ax.set_xlabel("")
    ax.set_ylabel("Predicted flux normalized \n to ouabain" if ax in [axes[0], axes[2]] else "")
    ax.set_xticks(range(len(flux_order)))
    ax.set_xticklabels(flux_order, rotation=45, ha="right")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.text(-0.12, 1.05, letters[i],
            transform=ax.transAxes,
            fontsize=16,
            fontweight='bold',
            va='top')

# Optional: add ONE legend for the line
axes[0].plot([], [], linestyle='--', color='black', label='O reference')
axes[0].legend(frameon=False)



plt.tight_layout()
plt.savefig('../01_baseline/01_baseline_subplots.png', dpi=300)
plt.show()
