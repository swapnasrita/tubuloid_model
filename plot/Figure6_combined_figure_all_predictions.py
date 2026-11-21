# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:59:25 2026

@author: swapn
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# -----------------------------
# COMMON FUNCTIONS
# -----------------------------
segments = ['PT', 'DCT', 'mTAL', 'cTAL', 'CCD']
percent_seg = [1, 1, 1, 1]

def extract_Na(folder_prefix, suffix):
    data = []
    for seg in segments:
        fname = f"{folder_prefix}{seg}_{suffix}.txt"
        with open(fname, 'r') as file:
            words = [line.strip().split() for line in file]
        for i, line in enumerate(words):
            if 'Na' in line:
                data.append(float(words[i+1][1]))

    # merge TAL
    data[2] += data[3]
    data[3] = data[4]

    data = data[:4]
    data = [data[i]*percent_seg[i] for i in range(len(percent_seg))]
    return data

# -----------------------------
# CREATE FIGURE
# -----------------------------
fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.2])

ax1 = fig.add_subplot(gs[0, 0])  # top-left
ax2 = fig.add_subplot(gs[0, 1])  # top-right
ax3 = fig.add_subplot(gs[1, :])  # bottom (spans both columns)

axes = [ax1, ax2, ax3]

# ============================================================
# (a) CODE 1 → SOLO INHIBITION
# ============================================================
ax = axes[0]

base = '../08_solo_inhibition/outlets_isolated_cell/'

control = extract_Na(base + 'C_outlet_isolated_cell_', 'male_hum')
O = extract_Na(base + 'O_outlet_isolated_cell_', 'male_hum')
T = extract_Na(base + 'T_outlet_isolated_cell_', 'male_hum')
A = extract_Na(base + 'A_outlet_isolated_cell_', 'male_hum')
B = extract_Na(base + 'B_outlet_isolated_cell_', 'male_hum')

y = np.array([sum(control), sum(T), sum(A), sum(B)]) / sum(O) * 100

ax.bar(range(len(y)), y, color='w', width=0.5, ec = 'k')
ax.set_xticks(range(len(y)))
ax.set_xticklabels(['None', 'T', 'A', 'B'],  fontsize = 18)
ax.set_ylabel('Predicted Na uptake (%)', fontsize = 18)
# ax.set_title('Solo inhibition')

# ============================================================
# (b) CODE 2 → TRANSPORTER CONTRIBUTION
# ============================================================
ax = axes[2]

transp = {
    "cTAL": ["NaKATPase", "NKCC2A", "NKCC2B", "NHE3"],
    "mTAL": ["NaKATPase", "NKCC2A", "NKCC2F", "NHE3"],
}

columns = ['1', '0.5', '0.3', '0.1', '0.01']
df = pd.DataFrame(columns=columns)

for num, i in enumerate(['1', '0p5', '0p3', '0p1', '0p01']):
    folder = f'../07_upreg_NKCC2/{i}_C_Male_hum_isolated_cell/'

    for key in transp:
        for value in transp[key]:
            k = key + '-' + value
            df.loc[k, columns[num]] = 0

            for root, dirs, files in os.walk(folder):
                for file in files:
                    if all(word in file for word in [key, value]):
                        if 'Na' in file.split('_')[4]:
                            with open(folder + file, 'r') as f:
                                df.loc[k, columns[num]] += float(f.read())

df_norm = df.div(df['1'], axis=0)
plot_data = np.array(df_norm).T

ind = np.arange(len(df))
width = 0.15
offset = [-2, -1, 0, 1, 2]
colors = ['#000000', '#333333', '#666666', '#999999', '#BBBBBB']

for i in range(len(columns)):
    ax.bar(ind + offset[i]*width, plot_data[i], width, color=colors[i])
    
labels = [0, 50, 70, 90, 99]
hatches = ['', 'x', '-', '//', '\\']

for i in range(len(columns)):
    ax.bar(
        ind + offset[i]*width,
        plot_data[i],
        width,
        color='w',
        ec='k',
        label=labels[i],
        hatch = hatches[i]# 👈 important
    )
ax.legend(
    title="% inhibition",
    bbox_to_anchor=(0.001, 1),
    loc="upper left",
    fontsize=14,
    ncol=5,
    frameon=False
)

ax.set_xticks(ind)
ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize = 18)
ax.set_ylabel(r' $\frac{Predicted Na_{flux}}{Predicted Na_{ flux,baseline}}$', fontsize = 18)
# ax.set_title('Transporter contribution')

# ============================================================
# (c) CODE 3 → MALE vs FEMALE
# ============================================================
ax = axes[1]

base_m = '../01_baseline/outlets_isolated_cell/'
base_f = '../05_female_baseline/outlets_isolated_cell/'

control = extract_Na(base_m + 'C_outlet_isolated_cell_', 'male_hum')
O = extract_Na(base_m + 'O_outlet_isolated_cell_', 'male_hum')
BO = extract_Na(base_m + 'BO_outlet_isolated_cell_', 'male_hum')
TO = extract_Na(base_m + 'TO_outlet_isolated_cell_', 'male_hum')
AO = extract_Na(base_m + 'AO_outlet_isolated_cell_', 'male_hum')
BTAO = extract_Na(base_m + 'BTAO_outlet_isolated_cell_', 'male_hum')

control_f = extract_Na(base_f + 'C_outlet_isolated_cell_', 'Female_hum')
O_f = extract_Na(base_f + 'O_outlet_isolated_cell_', 'Female_hum')
BO_f = extract_Na(base_f + 'BO_outlet_isolated_cell_', 'Female_hum')
TO_f = extract_Na(base_f + 'TO_outlet_isolated_cell_', 'Female_hum')
AO_f = extract_Na(base_f + 'AO_outlet_isolated_cell_', 'Female_hum')
BTAO_f = extract_Na(base_f + 'BTAO_outlet_isolated_cell_', 'Female_hum')

y_m = np.array([sum(control), sum(O), sum(BO), sum(TO), sum(AO), sum(BTAO)]) / sum(O) * 100
y_f = np.array([sum(control_f), sum(O_f), sum(BO_f), sum(TO_f), sum(AO_f), sum(BTAO_f)]) / sum(O_f) * 100

r = np.arange(len(y_m))
width = 0.35

ax.bar(r - width/2, y_m, width, label='Male', color='w', ec = 'k')
ax.bar(r + width/2, y_f, width, label='Female', color='w', ec = 'k', hatch = '//')

ax.set_xticks(r)
ax.set_xticklabels(['C', 'O', 'BO', 'TO', 'AO', 'BTAO'], fontsize = 18)
ax.set_ylabel('Predicted Na uptake (%)', fontsize = 18)
# ax.set_title('Sex comparison')
ax.legend(frameon=False)

# -----------------------------
# COMMON STYLING
# -----------------------------
ax1.tick_params(axis='y', labelsize=18)
ax2.tick_params(axis='y', labelsize=18)
ax3.tick_params(axis='y', labelsize=18)
for i, ax in enumerate(axes):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # subplot labels (a), (b), (c)
ax1.text(-0.2, 0.95, '(a)',
            transform=ax1.transAxes,
            fontsize=18,
            va='top')
ax2.text(-0.24, 0.95,'(b)',
            transform=ax2.transAxes,
            fontsize=18,
            va='top')
ax3.text(-0.09, 0.95,'(c)',
            transform=ax3.transAxes,
            fontsize=18,
            va='top')

plt.tight_layout()
plt.savefig('combined_figure_all_predictions.png', dpi=300)
plt.show()
