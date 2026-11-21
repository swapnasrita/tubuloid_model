# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 18:13:30 2023

@author: P70073624
"""

#%%
from defs import *
from values import *
from set_params import set_torq_params
import math
import flux
import copy
import numpy as np
from numpy.linalg import det
#%%
# parameters
dx = 1e-3
PI = 3.14159265
cell  = membrane()
celln = membrane()

count=0

def conservation_init (my_cell0, my_cell1 ,my_celln, my_dx):

    # global parameters
    global dx
    global cell1,cell1,cellm,celln

    cell0  = copy.deepcopy(my_cell0)
    cell1  = copy.deepcopy(my_cell1)
    celln = copy.deepcopy(my_celln)
    dx    = my_dx


#---------------------------------------------------------------------------
# compute residuals of conservation equations
# inputs:
#  cell:  current cell values
#  celln: cell function values at previous step
#  dx:    step
# returns:
#  fvec:  residuals
#
# no. of unknowns = no. of solutes, in cell and in LIS, plus 4 membrane
#  potentials (lumen-cell/LIS, cell/LIS-bath)
#---------------------------------------------------------------------------

def conservation_eqs (x):


    Jvol1=np.zeros(NC*NC).reshape(NC,NC)
    Jsol1=np.zeros(NS*NC*NC).reshape(NS,NC,NC)    


    if cell1.segment:
    
        if np.isnan(x).any():
            print("nan")


        #cell1.conc[:,0] = x[0:NS] 
        
        cell1.conc[0:3,0] = x[0:3]
        cell1.conc[0:3,1] = x[3:6]
        cell1.conc[0:3,4] = x[6:9]
        cell1.conc[0:3,5] = x[9:12]
     
        cell1.vol[0] = x[12] 
        cell1.vol[1] = x[13]
        cell1.vol[4] = x[14]
    
        cell1.ep[0] = x[15] 
        cell1.ep[1] = x[16]
        cell1.ep[4] = x[17]
        
        #cell1.pres[0] = x[3*NS+6]
    
        ph = np.zeros(NC)
    
        for j in [0,1,4,5]:
            ph[j] = -np.log(abs(cell1.conc[11][j])/1e3)/np.log(10.0)
        
        Bm = PI*cell1.diam
        Am = PI*(cell1.diam**2)/4
        

        fsol0 = np.zeros(3)
        fsol1 = np.zeros(3)
    
        for j in range(3):
            fsol0[j] = cell1.conc[j][0]*cell1.vol[0]*Vref/href 
            fsol1[j] = cell1.conc[j][0]*cell1.vol[0]*Vref/href
    
               
        fvol0 = cell1.vol[0]*Vref
        fvol1 = cell1.vol[0]*Vref
    
        PM0 = cell1.pres[0] 
        PM1 = cell1.pres[0] #x[3*NS+6]

        Jvol1,Jsol1, _ = flux.compute_fluxes(cell1)

    #---------------------------------------------------------------------
    #   Initialize  source terms
    #---------------------------------------------------------------------
        fvec = np.zeros(18)          # conservation for individual solutes
        S    = np.zeros(18)          # reaction terms
                                    
    #---------------------------------------------------------------------
    #---------------------------------------------------------------------
    #   COMPUTE SOURCE TERMS FOR FLOW OF SOLUTE I IN THE LUMEN
    #   The non-dimensional volume flux needs to be multiplied by 
    #    (Vref/href)
    #---------------------------------------------------------------------
    #---------------------------------------------------------------------
        
        for j in range(3):
            sumJsol1 = Jsol1[j,0,1]+Jsol1[j,0,4]
            fsola = cell1.vol[0]*cell1.conc[j][0]*Vref/href # NOT USED IN CALCULATION
            fsolb = cell1.vol[0]*cell1.conc[j][0]*Vref/href # NOT USED IN CALCULATION
            diff = fsolb - fsola
            if diff !=0 :
                print("fsol diffs: {}".format(fsolb-fsola))
            S[4*j] = Bm*cell1.len*sumJsol1/cell1.total
        
    #---------------------------------------------------------------------
    #---------------------------------------------------------------------
    #    COMPUTE SOURCE TERMS FOR SOLUTE I IN EACH CELLULAR COMPARTMENT
    #---------------------------------------------------------------------
    #---------------------------------------------------------------------
             
        
            
        for j in range(3):
            S[4*j+1] = -Jsol1[j][0][1]+Jsol1[j][1][4]+Jsol1[j][1][5]
            S[4*j+2] = -Jsol1[j][1][4]-Jsol1[j][0][4]+Jsol1[j][4][5]
            #[2*j] =  Jsol1[j][0][1] +  Jsol1[j][0][4]
            #S[2*j+1] =  Jsol1[j][1][5]   + Jsol1[j][4][5]
        # print(Jsol1[2,1,4],Jsol1[2,1,5],Jsol1[2,0,1])
        
     #---------------------------------------------------------------------
     #---------------------------------------------------------------------
     #   COMPUTE SOURCE TERMS FOR FLOW OF SOLUTE I IN THE BATH
     #   The non-dimensional volume flux needs to be multiplied by 
     #    (Vref/href)
     #---------------------------------------------------------------------
     #---------------------------------------------------------------------
        for j in range(3):
             sumJsol1 = -Jsol1[j,1,5]-Jsol1[j,4,5]
             fsola = cell1.vol[0]*cell1.conc[j][5]*Vref/href # NOT USED IN CALCULATION
             fsolb = cell1.vol[0]*cell1.conc[j][5]*Vref/href # NOT USED IN CALCULATION
             diff = fsolb - fsola
             if diff !=0 :
                 print("fsol diffs: {}".format(fsolb-fsola))
             S[4*j+3] = Bm*cell1.len*sumJsol1/cell1.total
        
  
    #---------------------------------------------------------------------
    #---------------------------------------------------------------------
    #   COMPUTE SOURCE TERMS FOR VOLUME FLOW IN THE LUMEN
    #   The non-dimensional volume flux needs to be multiplied by 
    #   Vref/(Pfref.Vwbar.Cref)
    #---------------------------------------------------------------------
    #---------------------------------------------------------------------
        
        fvmult = Pfref*Vwbar*Cref
    #        fvmult = Pfref*RTosm*Cref
    #        print('fvmult = %e' %fvmult)
    #        print('Jvol1[0,1] = %e'%Jvol1[0,1])
    #        print('Jvol1[0,4] = %e' %Jvol1[0,4])
        sumJvol1 = (Jvol1[0,1]+Jvol1[0,4])*fvmult
    #        print('sumJvol1 = %e'%sumJvol1)
    #        print('cell1.vol[0] = %f'%cell1.vol[0])
    #        print('Vref = %e'%Vref)
        fvola = cell1.vol[0]*Vref
    #        print('fvola = %e \n' %fvo dla)
        fvolb = cell1.vol[0]*Vref
    #        print('equations.py: cell1.len = %f' %cell1.len)
    #        print('equations.py: cell1.total = %f' %cell1.total)
        diff_vol = fvolb - fvola  
        if diff_vol !=0 :
            print("Vol diff {}".format(diff_vol))

        S[12] =  Bm*cell1.len*sumJvol1/cell1.total
        
        
    #---------------------------------------------------------------------
    #---------------------------------------------------------------------
    #   COMPUTE SOURCE TERMS FOR VOLUME IN EACH CELLULAR COMPARTMENT
    #---------------------------------------------------------------------
    #---------------------------------------------------------------------

        S[13] =   Jvol1[1,4]+Jvol1[1,5]-Jvol1[0,1]
        S[14] = Jvol1[4,5]-Jvol1[0,4]-Jvol1[1,4]
        
        #S[2*NS] =   Jvol1[0,1]+Jvol1[0,4]
        #S[1+2*NS] = Jvol1[4,5]-Jvol1[0,4]

    #---------------------------------------------------------------------
    #---------------------------------------------------------------------
    #   Convert to equations of the form F(X) = 0
    #---------------------------------------------------------------------
    #---------------------------------------------------------------------

    #---------------------------------------------------------------------
    #    For I = 1-3 and 9 (non-reacting solutes)
    #---------------------------------------------------------------------
        for m in range(12):
            fvec[m]=S[m]

    #     fvec[32]=S[32]
    #     fvec[33]=S[33]
    #     fvec[34]=S[34]
    #     fvec[35]=S[35]

    # #HCO2-/H2CO2 equations modified

    #     for m in range(12,15):
    #         fvec[4*m]=S[4*m]
    #         fvec[1+4*m]=S[1+4*m]
    #         fvec[2+4*m]=S[2+4*m] 
    #         fvec[3+4*m]=S[3+4*m] 


    # #---------------------------------------------------------------------
    # #    For CO2/HCO3/H2CO3
    # #    The dimensional factor for the kinetic term in the lumen is 
    # #    Am/Bm = Pi R^2 / 2 PI R = R / 2 = D / 4
    # #---------------------------------------------------------------------

    #     # conservation of reacting solutes in C and P: 
    #     fvec[12] = S[12]+S[16]+S[20] 
    #     fvec[13] = S[13]+S[17]+S[21] 
    #     fvec[14] = S[14]+S[19]+S[22]
    #     fvec[15] = S[15]+S[20]+S[23] 

    #     fvec[16] = ph[0]-pKHCO3-np.log(abs(float(cell1.conc[3,0])/float(cell1.conc[4,0])))/np.log(10.0)
    #     fvec[17] = ph[1]-pKHCO3-np.log(abs(float(cell1.conc[3,1])/float(cell1.conc[4,1])))/np.log(10.0)
    #     fvec[18] = ph[4]-pKHCO3-np.log(abs(float(cell1.conc[3,4])/float(cell1.conc[4,4])))/np.log(10.0)
    #     fvec[19] = ph[5]-pKHCO3-np.log(abs(float(cell1.conc[3,5])/float(cell1.conc[4,5])))/np.log(10.0)



    #     fkin1=(cell1.dkh[0]*cell1.conc[5,0]-cell1.dkd[0]*cell1.conc[4,0])
    #     fkin2=(cell1.dkh[0]*cell1.conc[5,0]-cell1.dkd[0]*cell1.conc[4,0])
    #     fvec[20] = S[20] + Am*cell1.len*fkin2/cell1.total/href

    #     facnd=Vref/href
    #     fvec[21] = S[21]+cell1.vol[1]*(cell1.dkh[1]*cell1.conc[5,1]-cell1.dkd[1]*cell1.conc[4,1])*facnd
    #     fvec[22] = S[22]+max(cell1.vol[4],cell1.volref[4])*(cell1.dkh[4]*cell1.conc[5,4]-cell1.dkd[4]*cell1.conc[4,4])*facnd
        
    #     #since default reaction rates for HCO3, H2Co3 is zero
    #     fkin1 = 0
    #     fkin2 = 0
    #     fvec[23] = S[23]

    # #---------------------------------------------------------------------
    # #    For HPO4(2-)/H2PO4(-)
    # #---------------------------------------------------------------------

    #     fvec[24] = S[24]+S[28]
    #     fvec[25] = S[25]+S[29]
    #     fvec[26] = S[26]+S[30]
    #     fvec[27] = S[27]+S[31]
    #     #fvec[20] = S[20]+S[23]
    #     fvec[28] = ph[0]-pKHPO4-np.log(abs(float(cell1.conc[6,0])/float(cell1.conc[7,0])))/np.log(10.0)
    #     fvec[29] = ph[1]-pKHPO4-np.log(abs(float(cell1.conc[6,1])/float(cell1.conc[7,1])))/np.log(10.0)
    #     fvec[30] = ph[4]-pKHPO4-np.log(abs(float(cell1.conc[6,4])/float(cell1.conc[7,4])))/np.log(10.0)
    #     fvec[31] = ph[5]-pKHPO4-np.log(abs(float(cell1.conc[6,5])/float(cell1.conc[7,5])))/np.log(10.0)

    # #---------------------------------------------------------------------
    # #    For NH3/NH4
    # #---------------------------------------------------------------------
    #     # TS = 1.3
    #     # scaleT = 1.0

    #     #     #torque-modulated effects
        
    #     #     # Note: 
    #     #     # Using attributes of cell1 over the attributes of cell1 makes very little difference in the computation. 
    #     #     # However, from the nature of Newton's method, attributes of cell1 are the initial guesses that we started with
    #     #     # and so we'll be using the same attributes for computation throughout.
            
    #     # PM=cell1.pres[0]

    #     # Radref,torqR,torqvm,PbloodPT,torqL,torqd = set_torq_params(cell1.species,cell1.sex)

    #     # if cell1.species == 'hum':
    #     #     fac1 = 8.0*visc*(cell1.volref[0]*Vref)*torqL/(Radref**2) 
    #     # elif cell1.species == 'rat':
    #     #     fac1 = 8.0*visc*(cell1.vol_init[0]*Vref)*torqL/(Radref**2)
    #     # elif cell1.species == 'mou':
    #     #     fac1 = 8.0*visc*(cell1.vol_init[0]*Vref)*torqL/(Radref**2)
    #     # else:
    #     #     print('cell1.species: ' + str(cell1.species))
    #     #     raise Exception('what is species?')
    #     # fac2 = 1.0 + (torqL+torqd)/Radref + 0.50*((torqL/Radref)**2)
    #     # TM0= fac1*fac2

    #     # RMtorq = torqR*(1.0e0+torqvm*(PM - PbloodPT)) 
    #     # factor1 = 8.0*visc*(cell1.vol[0]*Vref)*torqL/(RMtorq**2) 
    #     # factor2 = 1.0 + (torqL+torqd)/RMtorq + 0.50*((torqL/RMtorq)**2)
    #     # Torque = factor1*factor2

    #     # Scaletorq = 1.0 + TS*scaleT*(Torque/TM0-1.0)

    #     # Qnh4 = 0.25e-6/(href*Cref)*Scaletorq

    #     #else:
    #     #    Qnh4 = 0.0

    #     fvec[36] = S[36]+S[40] #- Qnh4
    #     fvec[37] = S[37]+S[41] 
    #     fvec[38] = S[38]+S[42]
    #     fvec[39] = S[39]+S[43]
        
    #     #fvec[19] = S[19]+S[21]-Qnh4

    #     fvec[40] = ph[0]-pKNH3-np.log(abs(float(cell1.conc[9,0])/float(cell1.conc[10,0])))/np.log(10.0)
    #     fvec[41] = ph[1]-pKNH3-np.log(abs(float(cell1.conc[9,1])/float(cell1.conc[10,1])))/np.log(10.0)
    #     fvec[42] = ph[4]-pKNH3-np.log(abs(float(cell1.conc[9,4])/float(cell1.conc[10,4])))/np.log(10.0)
    #     fvec[43] = ph[5]-pKNH3-np.log(abs(float(cell1.conc[9,5])/float(cell1.conc[10,5])))/np.log(10.0)

    #     # print(ph[0],pKNH3,cell1.conc[9,0],cell1.conc[10,0])
    # #---------------------------------------------------------------------
    # #    For HCO2-/H2CO2
    # #---------------------------------------------------------------------

    #     fvec[48] = S[48]+S[52]
    #     fvec[49] = S[49]+S[53]  
    #     fvec[50] = S[50]+S[54]
    #     fvec[51] = S[51]+S[55]
        
    #     fvec[52] = ph[0]-pKHCO2-np.log(abs(cell1.conc[12,0]/cell1.conc[13,0]))/np.log(10.0)
    #     fvec[53] = ph[1]-pKHCO2-np.log(abs(cell1.conc[12,1]/cell1.conc[13,1]))/np.log(10.0)
    #     fvec[54] = ph[4]-pKHCO2-np.log(abs(cell1.conc[12,4]/cell1.conc[13,4]))/np.log(10.0)
    #     fvec[55] = ph[5]-pKHCO2-np.log(abs(cell1.conc[12,5]/cell1.conc[13,5]))/np.log(10.0)
    #     #print(ph[1],pKHCO2,cell1.conc[12,1],cell1.conc[13,1])

    # #---------------------------------------------------------------------
    # #    For pH
    # #---------------------------------------------------------------------
    #     fvec[44] = S[44]+S[40]-S[24]-S[12]-S[48]
    #     fvec[45] = S[45]+S[41]-S[25]-S[13]-S[49]    # H + NH4 - HPO4 - HCO3 - HCO2
    #     fvec[46] = S[46]+S[42]-S[26]-S[14]-S[50]   
    #     fvec[47] = S[47]+S[43]-S[27]-S[15]-S[51]

    #---------------------------------------------------------------------
    #    For volume
    #---------------------------------------------------------------------

        
        fvec[12]=S[12]
        fvec[13]=S[13]
        fvec[14]=S[14]
        
        


    #---------------------------------------------------------------------72
    #     For EP, need to satisfy electroneutrality
    #---------------------------------------------------------------------72

        currM = 0
        for j in range(3):
            currM = currM+zval[j]*(Jsol1[j,0,1]+Jsol1[j,0,4])
        fvec[15]=currM

        volP=cell1.volref[1]/cell1.vol[1]
        CimpP=cell1.cimpref[1]*volP
        facP=np.exp(np.log(10.0)*(ph[1]-pKbuf))
        
        # Changed by Dania:
        CbufP=cell1.cbuftot[1]*volP*facP/(facP+1)

        elecP= -CimpP-CbufP
        elecE=0

        for j in range(3):
            elecP=elecP+zval[j]*cell1.conc[j,1]
            elecE=elecE+zval[j]*cell1.conc[j,4]

        fvec[16]=elecP

        fvec[17]=elecE
    
#---------------------------------------------------------------------
#     For Pressure
#---------------------------------------------------------------------  
        '''
       # PT and S3 are the only compliant segments.
        if cell1.segment == 'PT' or cell1.segment == 'S3':       
            Radref,torqR,torqvm,PbloodPT,torqL,torqd = set_torq_params(cell1.species,cell1.sex,cell1.preg)
       
        if cell1.species == 'hum':
            if cell1.segment == 'PT' or cell1.segment == 'S3':
                RMcompl = torqR*(1.0e0+torqvm*(cell1.pres[0] - PbloodPT))
                Amcompl = PI*(RMcompl**2)
                factor1 = 8.0*PI*visc/(Amcompl**2)
            else:
                factor1 = 8.0*PI*visc/(Am**2)
        elif cell1.species == 'rat':
            if cell1.segment == 'PT' or cell1.segment == 'S3':
                RMcompl = torqR*(1.0e0+torqvm*(cell1.pres[0] - PbloodPT))
                Amcompl = PI*(RMcompl**2)
                factor1 = 8.0*PI*visc/(Amcompl**2)
            else:
                factor1 = 8.0*PI*visc/(Am**2)

        elif cell1.species == 'mou':
            if cell1.segment == 'PT' or cell1.segment == 'S3':
                RMcompl = torqR*(1.0e0+torqvm*(cell1.pres[0] - PbloodPT))
                Amcompl = PI*(RMcompl**2)
                factor1 = 8.0*PI*visc/(Amcompl**2)
            else:
                factor1 = 8.0*PI*visc/(Am**2)
        else:
            print('cell1.species: ' + str(cell1.species))
            raise Exception('what is species?')

        if cell1.segment == 'IMCD':
            if cell1.species == 'rat':
                factor1 = factor1*coalesce*2
            elif cell1.species == 'mou':
                factor1 = factor1*coalesce*2
            elif cell1.species == 'hum':
                factor1 = factor1*coalesce*1.0
            else:
                print('cell1.species: ' + str(cell1.species))
                raise Exception('what is species?')

        #fvec[6+3*NS] = PM1-PM0+factor1*cell1.vol[0]*Vref*cell1.len/cell1.total
        #print(PM1,PM0,factor1,cell1.vol[0])
        '''

    else:

        raise Exception('what is this segment? ' + cell1.segment)

    #print("Active Transport For NH3: {}".format([ electro_jsol[9][0,1], ))
    print(fvec)
    return fvec, Jsol1, Jvol1









# %%
