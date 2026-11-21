# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 12:30:32 2023

@author: P70073624
"""

#%%
from defs import *
import copy
import set_params
import numpy as np
import equations
from values import *
import Newton
import timeit
import os


from set_params import set_torq_params
import output
import argparse
import pdb

#%%
solute = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
compart = ['Lumen','Cell','ICA','ICB','LIS','Bath']
cw=Vref*60e6

# parsing arguments:
# required input
parser = argparse.ArgumentParser()
parser.add_argument('--input_file', required = False,  default = 'cTALparams_M_hum.dat',type = str,help = 'cell number')

# optional input:
parser.add_argument('--sex',choices=['Male','Female'],default = 'Male',type = str,help = 'Sex')
parser.add_argument('--species',choices=['human','rat','mouse'],default = 'rat',type = str, help = 'Human, Rat, or Mouse model')
parser.add_argument('--type',choices = ['superficial','multiple'],default = 'superficial',type=str,help='superficial nephron or multiple nephrons?')
parser.add_argument('--inhibition',choices=['ACE','SGLT2','NHE3-50','NHE3-80','NKCC2-70','NKCC2-100','NCC-70','NCC-100','ENaC-70','ENaC-100','SNB-70','SNB-100'],default = None,type = str,help = 'any transporter inhibition?')
segments = ['PT', 'DCT', 'mTAL', 'cTAL', 'CNT', 'CCD']
inhibitions = ['C', 'O', 'BO', 'TO', 'AO', 'BTAO']

for inhib in inhibitions:
    for seg in segments:
        fname = seg + 'params_M_hum.dat'
        args = parser.parse_args()
        sex = 'Male' #args.sex
        species = 'human' #args.species
        sup_or_multi = args.type
        # inhib = args.inhibition
        input_filename = fname
        
        print("Input file is {}".format(input_filename))
        #%%
        # conditions check:
        if inhib != None:
            file_to_save = inhib+'_'+sex+'_'+species[0:3]
        else:
            file_to_save = sex + '_' + species[0:3] +'_normal'
            
        file_to_save = file_to_save + "_isolated_cell"
        
        if os.path.isdir(file_to_save) == False:
            os.makedirs(file_to_save)
            
        if sup_or_multi == 'superficial':
            parts = 'sup'
        else:
            parts = ['sup','jux1','jux2','jux3','jux4','jux5']
        
        
        
        #%%
        def compute_segment(
                    sup_or_jux,
                    sex,
                    species,
                    input_filename,
                    file_to_save):
        
        
            solute = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
            compart = ['Lumen','Cell','ICA','ICB','LIS','Bath']
            cw=Vref*60e6
            
            #========================================================
            # Proximal convolute tubule
            #========================================================
            # if sex == 'Male':
            #     filename = './datafiles/PTparams_M_'+species[0:3]+'_'+ input_filename +'.dat'
            # elif sex == 'Female':
            #     filename = './datafiles/PTparams_F_'+species[0:3]+'_'+ input_filename +'.dat'
            # else:
            #     filename ='./datafiles/PTparams_F_'+species[0:3]+'_'+ input_filename +'.dat'
        
            
            print("\n**** Finding Concentrations, Volumes and Membrane Potentials in {} Cell ****\n".format(
                input_filename[0:-16]
                ))
            
            filename = './datafiles/'+ input_filename
            
            
            pt=compute(
                    filename,  
                    'Newton',
                    sup_or_jux=None,
                    species = species,
                    inhibition=inhib,
                    unx = 'N')
        
            
            # for j in range(1):
            #     if pt.segment == 'PT':
            #         TS = 1.3
            #         scaleT = 1.0
            #     elif pt.segment == 'S3':
            #         TS = 1.3
            #         scaleT = 0.5
        
            #     #torque-modulated effects                
            #     PM=pt.pres[0]                
            #     Radref,torqR,torqvm,PbloodPT,torqL,torqd = set_torq_params(pt.species,pt.sex)
                
                                            
                # if pt.species == 'rat':
                #     fac1 = 8.0*visc*(pt.vol_init[0]*Vref)*torqL/(Radref**2)
                # elif pt.species == 'mou':
                #     fac1 = 8.0*visc*(pt.vol_init[0]*Vref)*torqL/(Radref**2)
                # elif pt.species == 'hum':
                #     fac1 = 8.0*visc*(pt.volref[0]*Vref)*torqL/(Radref**2)
                # else:
                #     print('pt.species: ' + str(pt.species))
                #     raise Exception('what is species?')
                # fac2 = 1.0 + (torqL+torqd)/Radref + 0.50*((torqL/Radref)**2)
                # TM0= fac1*fac2
                    
                # RMtorq = torqR*(1.0e0+torqvm*(PM - PbloodPT))
        
                # factor1 = 8.0*visc*(pt.vol[0]*Vref)*torqL/(RMtorq**2) 
                # factor2 = 1.0 + (torqL+torqd)/RMtorq + 0.50*((torqL/RMtorq)**2)
                # Torque = factor1*factor2
                    
            Scaletorq = 1.0 #SS+ TS*scaleT*(Torque/TM0-1.0)
        # for calculations of torque, https://www.notion.so/Layton2019_epithelial-solute-water-transport-model-254253e63a3347a8b12f90d09d9e217d?pvs=4
            
            if os.path.isdir(file_to_save) == False:
                    os.makedirs(file_to_save)
            output.output_segment_results(pt,sup_or_jux,Scaletorq,file_to_save,1)
        
            return pt
        
        #%%
        def compute(filename,  
                    method,
                    sup_or_jux=None,
                    species = 'rat',
                    inhibition=inhib,
                    unx = 'N'):
        
            start= timeit.default_timer()
            cell = membrane()
            
            # initialize the species of cell
            if species == 'human':
                    cell.species = 'hum'
            elif species == 'rat':
                    cell.species = 'rat'
            elif species == 'mouse':
                    cell.species = 'mou'
            # the diabetic status of cell.
            else:
                raise Exception('what is species', species)
            # superficial nephron or juxtamedullary nephron
            cell.type = sup_or_jux
            
            if inhibition == None:
                cell.inhib = ''
            else:
                cell.inhib = inhib
        
            cell.unx = unx
        
            # initialize intersitital concentration and read parameters.
            set_params.read_params_isolated(cell,filename)
            
            # initial guess of unknowns
        
            celln = copy.deepcopy(cell)
            dx = 1.0e-3
                    
            x = np.zeros(18)
        
        
            
            x[0:3] = cell.conc[0:3,0]
            x[3:6] = cell.conc[0:3,1]
            x[6:9] = cell.conc[0:3,4]
            x[9:12] = cell.conc[0:3,5]
            
            x[12] = cell.vol[0]
            x[13] = cell.vol[1]
            x[14] = cell.vol[4]
            
            
    
            x[15] = cell.ep[0]
            x[16] = cell.ep[1]
            x[17] = cell.ep[4]
    
            
        
            equations.conservation_init (cell,cell,celln,dx)
            
            fvec = equations.conservation_eqs (x)
            
            
        
            # solving the system
            if method == 'Newton':
                if species == 'human':
                    sol = Newton.newton_human(equations.conservation_eqs,x,0,cell)
                elif species == 'rat':
                        sol = Newton.newton_rat(equations.conservation_eqs,x,0,cell)
                elif species == 'mouse':
                    sol = Newton.newton_rat(equations.conservation_eqs,x,cell)
                else:
                    print('species:' + species)
                    raise Exception('human, rat or mouse?', species)
            elif method == 'Broyden':
                sol, funxk = Newton.broyden(equations.conservation_eqs,x,cell.segment)
                
            else:
                raise Exception('what is this method?', method)
            
               # print("Fvec is:{}".format(fvec))
                # set up next cell
                    
            
            cell.conc[0:3,0] = sol[0:3]
            cell.conc[0:3,1] = sol[3:6]
            cell.conc[0:3,4] = sol[6:9]
            cell.conc[0:3,5] = sol[9:12]
            
            cell.vol[0] = sol[12]
            cell.vol[1] = sol[13]
            cell.vol[4] = sol[14]
            
    
            cell.ep[0] = sol[15]
            cell.ep[1] = sol[16]
            cell.ep[4] = sol[17]
            
    
            
        
        
        #================================OUTPUT IN TO FILE================================        
                # directory to store outlet files
            if os.path.isdir('outlets_isolated_cell') == False:
                os.makedirs('outlets_isolated_cell')
        
            file = open('./outlets_isolated_cell/'+ inhib+ '_outlet_isolated_cell'+ '_' + filename[12:-16] + '_' +cell.sex +'_'+cell.species+ '.txt','w')
            file.write("\nConcentrations\n")
        
            for j in range(3):
                file.write("{} (Isolated):\n".format(solute[j]))
        
                file.write('{} {} {} {} \n\n'.format(
                    cell.conc[j,0],cell.conc[j,1],cell.conc[j,4],cell.conc[j,5]
                    ))
                
        
                #file.write('{} {} {} {} \n'.format())
        
            file.write("Volumes (Isolated):\n")
            file.write('{} {} {} \n\n'.format(
                cell.vol[0],cell.vol[1],cell.vol[4])
                )
        
        
        
            file.write("Membrane Potentials (Isolated):\n")
            file.write('{} {} {} \n\n'.format(
                cell.ep[0],cell.ep[1],cell.ep[4])
                )
        
        
            file.write("Luminal Pressure (Isolated):\n")
            file.write(str(cell.pres[0])+'\n\n')
            file.close()
        
        
        
            stop=timeit.default_timer()
            ComputationTime=stop-start
        
        
        
            file=open('ComputationTime.txt','a')
            file.write('Computation Time for %s %s: %f \n' %(cell.type,cell.segment,ComputationTime))
            file.close()
        
            print("\n** This cell is done **")
            print("\n** Results file saved in {} **\n\n".format('outlets_isolated_cell'))
        
            return cell
        # define desired cell to be outputed:
        #%%
        # load original fluxes:
        which_cell = input_filename[0:-16]
        # pdb.set_trace()
        pt = compute_segment(
                            parts,
                            sex,
                            species,
                            input_filename,
                            file_to_save 
                            )
    

# %%
