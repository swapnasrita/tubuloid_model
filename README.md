# About
This is the latest version of isolated proximal tubule cell. Using this file we can specify the luminal and basolateral concentrations and find the cell and LIS concentrations in that specific cell. 

# Instructions
To run the code use command: 
python isolated_pt_cell.py --input_file [option] --inhibition [option] --sex [option]

The options here are:

input_file: **'PTparams_M_rat_first.dat', 'PTparams_M_rat_mid.dat', 'PTparams_M_rat_last.dat'**

sex: **'Male', 'Female'**

inhibition: **ACE, SGLT2, NHE3-50, NHE3-80, NKCC2-70, NKCC2-100, NCC-70, NCC-100, ENaC-70, ENaC-100, SNB-70, SNB-100** (optional, default: None).


### Understanding output

All the output files' names are in following structure: 'Male_rat_normal_isolated_cell/male_rat_PT_apical_K_sup.txt'. 

The unit of concentration from outputs is **mmol/L (mM)**.

The unit of volume is **nl/min**.

The unit of flow is **pmol/min**.

**/plot/** contains some example scripts for plotting output

## Related Work
Please cite appropriate paper(s) when using this model.

Published papers:
* **multiple nephron model (male only):** [2016 Layton et al. "A computational model for simulating solute transport and oxygen consumption along the nephrons"](https://journals.physiology.org/doi/full/10.1152/ajprenal.00293.2016)
* **superficial nephron (sex-specific):** [2019 Hu et al. "Functional implications of the sex differences in transporter abundance along the rat nephron: modeling and analysis"](https://journals.physiology.org/doi/full/10.1152/ajprenal.00352.2019)
* **multiple nephron (sex-specific):** [2020 Hu et al. "Sex differences in solute transport along the nephrons: effects of Na+ transport inhibition"](https://journals.physiology.org/doi/abs/10.1152/ajprenal.00240.2020?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org)
* **human multiple nephron (male only):** [2019 Layton and Layton "A computational model of epithelial solute and water transport along a human nephron"](https://journals.plos.org/ploscompbiol/article?rev=2&id=10.1371/journal.pcbi.1006108)
* **human multiple nephron (sex-specific):** [2021 Hu et al. "Sex differences in solute and water handling in the human kidney: Modeling and functional implications"](https://www.sciencedirect.com/science/article/pii/S2589004221006350)
* **diabetic human (sex-specific):** [2021 Hu et al. "A Computational Model of Kidney Function in a Patient with Diabetes"](https://www.mdpi.com/1422-0067/22/11/5819)
* **pregnant rat (superficial nephron):** [2022 Stadt and Layton "Adaptive changes in single-nephron GFR, tubular morphology, and transport in a pregnant rat nephron: modeling and analysis"](https://journals.physiology.org/doi/abs/10.1152/ajprenal.00264.2021)

### Previous versions
Previous versions of this model code are available [here](https://github.com/uwrhu).
