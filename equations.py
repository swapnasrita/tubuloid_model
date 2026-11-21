#%%
from defs import *
from values import *
from set_params import set_torq_params
import math
import flux
import copy
import numpy as np
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

    
    if np.isnan(x).any():
        print("nan")


    #cell1.conc[:,0] = x[0:NS] 
    cell1.conc[0:3,1] = x[0:3]
    cell1.conc[0:3,4] = x[3:6]
 

    cell1.ep[0] = x[6] 
    cell1.ep[1] = x[7]
    cell1.ep[4] = x[8]

    #cell1.pres[0] = x[3*NS+6]

    ph = np.zeros(NC)

    for j in [0,1,4]:
        ph[j] = -np.log(abs(cell1.conc[11][j])/1e3)/np.log(10.0)
    
    # Bm = PI*cell1.diam
    # Am = PI*(cell1.diam**2)/4
    

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
    fvec = np.zeros(9)          # conservation for individual solutes
    S    = np.zeros(9)          # reaction terms
                                    
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#   COMPUTE SOURCE TERMS FOR FLOW OF SOLUTE I IN THE LUMEN
#   The non-dimensional volume flux needs to be multiplied by 
#    (Vref/href)
#---------------------------------------------------------------------
#---------------------------------------------------------------------
    '''
    for j in range(NS):
        sumJsol1 = Jsol1[j,0,1]+Jsol1[j,0,4]
        fsola = cell1.vol[0]*cell1.conc[j][0]*Vref/href
        fsolb = cell1.vol[0]*cell1.conc[j][0]*Vref/href
        diff = fsolb - fsola
        if diff !=0 :
            print("fsol diffs: {}".format(fsolb-fsola))
        S[j] = Bm*cell1.len*sumJsol1/cell1.total
    '''
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#    COMPUTE SOURCE TERMS FOR SOLUTE I IN EACH CELLULAR COMPARTMENT
#---------------------------------------------------------------------
#---------------------------------------------------------------------
         
    
        
    for j in range(3):
        S[2*j] = -Jsol1[j][0][1]+Jsol1[j][1][4]+Jsol1[j][1][5]
        S[2*j+1] = -Jsol1[j][1][4]-Jsol1[j][0][4]+Jsol1[j][4][5]
        #[2*j] =  Jsol1[j][0][1] +  Jsol1[j][0][4]
        #S[2*j+1] =  Jsol1[j][1][5]   + Jsol1[j][4][5]
    # print(Jsol1[2,1,4],Jsol1[2,1,5],Jsol1[2,0,1])
  
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#   COMPUTE SOURCE TERMS FOR VOLUME FLOW IN THE LUMEN
#   The non-dimensional volume flux needs to be multiplied by 
#   Vref/(Pfref.Vwbar.Cref)
#---------------------------------------------------------------------
#---------------------------------------------------------------------
    '''
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

        S[3*NS] =  Bm*cell1.len*sumJvol1/cell1.total
    '''
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#   COMPUTE SOURCE TERMS FOR VOLUME IN EACH CELLULAR COMPARTMENT
#---------------------------------------------------------------------
#---------------------------------------------------------------------

    # S[2*NS] =   Jvol1[1,4]+Jvol1[1,5]-Jvol1[0,1]
    # S[1+2*NS] = Jvol1[4,5]-Jvol1[0,4]-Jvol1[1,4]

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
    for m in range(6):
        fvec[m]=S[m]

    # fvec[16]=S[16]
    # fvec[17]=S[17]

#HCO2-/H2CO2 equations modified

    # for m in range(12,15):
    #     fvec[2*m]=S[2*m]
    #     fvec[1+2*m]=S[1+2*m]
        #fvec[2+3*m]=S[2+3*m] 


#---------------------------------------------------------------------
#    For CO2/HCO3/H2CO3
#    The dimensional factor for the kinetic term in the lumen is 
#    Am/Bm = Pi R^2 / 2 PI R = R / 2 = D / 4
#---------------------------------------------------------------------

    # conservation of reacting solutes in C and P: 
    # fvec[6] = S[6]+S[8]+S[10]
    # fvec[7] = S[7]+S[9]+S[11]
    #fvec[8] = S[8]+S[10]+S[12]

    #fvec[9] = ph[0]-pKHCO3-np.log(abs(float(cell1.conc[3,0])/float(cell1.conc[4,0])))/np.log(10.0)
    # fvec[8] = ph[1]-pKHCO3-np.log(abs(float(cell1.conc[3,1])/float(cell1.conc[4,1])))/np.log(10.0)
    # fvec[9] = ph[4]-pKHCO3-np.log(abs(float(cell1.conc[3,4])/float(cell1.conc[4,4])))/np.log(10.0)



    #fkin1=(cell1.dkh[0]*cell1.conc[5,0]-cell1.dkd[0]*cell1.conc[4,0])
    #fkin2=(cell1.dkh[0]*cell1.conc[5,0]-cell1.dkd[0]*cell1.conc[4,0])
    #fvec[15] = S[15] + Am*cell1.len*fkin2/cell1.total/href

    # facnd=Vref/href
    # fvec[10] = S[10]+cell1.vol[1]*(cell1.dkh[1]*cell1.conc[5,1]-cell1.dkd[1]*cell1.conc[4,1])*facnd
    # fvec[11] = S[11]+max(cell1.vol[4],cell1.volref[4])*(cell1.dkh[4]*cell1.conc[5,4]-cell1.dkd[4]*cell1.conc[4,4])*facnd

#---------------------------------------------------------------------
#    For HPO4(2-)/H2PO4(-)
#---------------------------------------------------------------------

    # fvec[12] = S[12]+S[14]
    # fvec[13] = S[13]+S[15]
    # #fvec[20] = S[20]+S[23]
    # #fvec[21] = ph[0]-pKHPO4-np.log(abs(float(cell1.conc[6,0])/float(cell1.conc[7,0])))/np.log(10.0)
    # fvec[14] = ph[1]-pKHPO4-np.log(abs(float(cell1.conc[6,1])/float(cell1.conc[7,1])))/np.log(10.0)
    # fvec[15] = ph[4]-pKHPO4-np.log(abs(float(cell1.conc[6,4])/float(cell1.conc[7,4])))/np.log(10.0)

#---------------------------------------------------------------------
#    For NH3/NH4
#---------------------------------------------------------------------
    # TS = 1.3
    # scaleT = 1.0

    #     #torque-modulated effects
    
    #     # Note: 
    #     # Using attributes of cell1 over the attributes of cell1 makes very little difference in the computation. 
    #     # However, from the nature of Newton's method, attributes of cell1 are the initial guesses that we started with
    #     # and so we'll be using the same attributes for computation throughout.
        
    # PM=cell1.pres[0]

    # Radref,torqR,torqvm,PbloodPT,torqL,torqd = set_torq_params(cell1.species,cell1.sex)

    # if cell1.species == 'hum':
    #     fac1 = 8.0*visc*(cell1.volref[0]*Vref)*torqL/(Radref**2) 
    # elif cell1.species == 'rat':
    #     fac1 = 8.0*visc*(cell1.vol_init[0]*Vref)*torqL/(Radref**2)
    # elif cell1.species == 'mou':
    #     fac1 = 8.0*visc*(cell1.vol_init[0]*Vref)*torqL/(Radref**2)
    # else:
    #     print('cell1.species: ' + str(cell1.species))
    #     raise Exception('what is species?')
    # fac2 = 1.0 + (torqL+torqd)/Radref + 0.50*((torqL/Radref)**2)
    # TM0= fac1*fac2

    # RMtorq = torqR*(1.0e0+torqvm*(PM - PbloodPT)) 
    # factor1 = 8.0*visc*(cell1.vol[0]*Vref)*torqL/(RMtorq**2) 
    # factor2 = 1.0 + (torqL+torqd)/RMtorq + 0.50*((torqL/RMtorq)**2)
    # Torque = factor1*factor2

    # Scaletorq = 1.0 + TS*scaleT*(Torque/TM0-1.0)

    # Qnh4 = 0.25e-6/(href*Cref)*Scaletorq

    #else:
    #    Qnh4 = 0.0

    # fvec[18] = S[18]+S[20] #- Qnh4
    # fvec[19] = S[19]+S[21] 
    
    # #fvec[19] = S[19]+S[21]-Qnh4

    # #fvec[30] = ph[0]-pKNH3-np.log(abs(float(cell1.conc[9,0])/float(cell1.conc[10,0])))/np.log(10.0)
    # fvec[20] = ph[1]-pKNH3-np.log(abs(float(cell1.conc[9,1])/float(cell1.conc[10,1])))/np.log(10.0)
    # fvec[21] = ph[4]-pKNH3-np.log(abs(float(cell1.conc[9,4])/float(cell1.conc[10,4])))/np.log(10.0)

    # print(ph[0],pKNH3,cell1.conc[9,0],cell1.conc[10,0])
#---------------------------------------------------------------------
#    For HCO2-/H2CO2
#---------------------------------------------------------------------

    # fvec[24] = S[24]+S[26]
    # fvec[25] = S[25]+S[27]  
    # # fvec[38] = S[38]+S[41]
    # # fvec[39] = ph[0]-pKHCO2-np.log(abs(cell1.conc[12,0]/cell1.conc[13,0]))/np.log(10.0)

    # fvec[26] = ph[1]-pKHCO2-np.log(abs(cell1.conc[12,1]/cell1.conc[13,1]))/np.log(10.0)
    # fvec[27] = ph[4]-pKHCO2-np.log(abs(cell1.conc[12,4]/cell1.conc[13,4]))/np.log(10.0)

    #print(ph[1],pKHCO2,cell1.conc[12,1],cell1.conc[13,1])

#---------------------------------------------------------------------
#    For pH
#---------------------------------------------------------------------
    # fvec[22] = S[22]+S[20]-S[12]-S[6]-S[24]   
    # fvec[23] = S[23]+S[21]-S[13]-S[7]-S[25]   
    #fvec[35] = S[35]+S[32]-S[20]-S[11]-S[38]

#---------------------------------------------------------------------
#    For volume
#---------------------------------------------------------------------

    # fvec[2*NS]=S[2*NS]
    # fvec[1+2*NS]=S[1+2*NS]
    #fvec[2+3*NS]=S[2+3*NS]


#---------------------------------------------------------------------72
#     For EP, need to satisfy electroneutrality
#---------------------------------------------------------------------72

    currM = 0
    for j in range(3):
        currM = currM+zval[j]*(Jsol1[j,0,1]+Jsol1[j,0,4])
    fvec[6]=currM

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

    fvec[7]=elecP

    fvec[8]=elecE

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

    # else:

    #     raise Exception('what is this segment? ' + cell1.segment)

    #print("Active Transport For NH3: {}".format([ electro_jsol[9][0,1], ))
        
    return fvec, Jsol1, Jvol1









# %%
