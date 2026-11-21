#%%
from values import *
import numpy as np
from defs import *


def female_conc(cell, i):
    # Concentrations of K, Cl are lower in female rat.
    cell.conc[1,5]=cell.conc[1,5]-1
    cell.conc[2,5]=cell.conc[2,5]-1
    if cell.segment == 'PT':
        if cell.preg == 'mid':
            # conc of Na, K, Cl, HCO3 different in MP rat plasma
            cell.conc[0,0] = cell.conc[0,0]*0.95
            cell.conc[1,0] = cell.conc[1,0]*1.0 #no change in K in MP
            cell.conc[2,0] = cell.conc[2,0]*0.95
            cell.conc[3,0] = cell.conc[3,0]*0.95

            # also change in starting for bath
            cell.conc[0,5] = cell.conc[0,0]
            cell.conc[1,5] = cell.conc[1,0] 
            cell.conc[2,5] = cell.conc[2,0]
            cell.conc[3,5] = cell.conc[3,0]
        elif cell.preg == 'late':
            # conc of Na, K, Cl, HCO3 different in LP rat plasma
            cell.conc[0,0] = cell.conc[0,0]*0.95
            cell.conc[1,0] = cell.conc[1,0]*1.15 # upregulated K in LP
            cell.conc[2,0] = cell.conc[2,0]*0.95
            cell.conc[3,0] = cell.conc[3,0]*0.9645

            # also change in starting for bath
            cell.conc[0,5] = cell.conc[0,0]
            cell.conc[1,5] = cell.conc[1,0]
            cell.conc[2,5] = cell.conc[2,0]
            cell.conc[3,5] = cell.conc[3,0]
    

def boundaryBath(cell,
                input_file = False
                ):

    if input_file == False:

        if cell.species == 'rat':
            if cell.inhib == 'NHE3-50':
                cell.oi[0] = TotSodOI_noinhib
                cell.oi[1] = TotPotOI_noinhib
                cell.oi[2] = TotCloOI_noinhib
                cell.oi[8] = TotureaOI_50inhib
                cell.pap[0] = TotSodPap_50inhib
                cell.pap[1] = TotPotPap_noinhib
                cell.pap[2] = TotCloPap_50inhib
                cell.pap[8] = TotureaPap_50inhib
            elif cell.inhib == 'NHE3-80':
                cell.oi[0] = TotSodOI_noinhib
                cell.oi[1] = TotPotOI_noinhib
                cell.oi[2] = TotCloOI_noinhib		
                cell.oi[8] = TotureaOI_80inhib
                cell.pap[0] = TotSodPap_80inhib
                cell.pap[1] = TotPotPap_80inhib
                cell.pap[2] = TotCloPap_80inhib
                cell.pap[8] = TotureaPap_80inhib
            elif cell.inhib == 'NKCC2-70':
                if cell.preg == 'non':
                    cell.oi[0] = TotSodOI_70NKCCinhib
                    cell.oi[1] = TotPotOI_70NKCCinhib
                    cell.oi[2] = TotCloOI_70NKCCinhib		
                    cell.oi[8] = TotureaOI_70NKCCinhib
                    cell.pap[0] = TotSodPap_70NKCCinhib
                    cell.pap[1] = TotPotPap_70NKCCinhib
                    cell.pap[2] = TotCloPap_70NKCCinhib
                    cell.pap[8] = TotureaPap_70NKCCinhib
                elif cell.preg == 'mid':
                    cell.oi[0] = TotSodOI_MP_70NKCCinhib
                    cell.oi[1] = TotPotOI_MP_70NKCCinhib
                    cell.oi[2] = TotCloOI_MP_70NKCCinhib		
                    cell.oi[8] = TotureaOI_70NKCCinhib
                    cell.pap[0] = TotSodPap_MP_70NKCCinhib
                    cell.pap[1] = TotPotPap_MP_70NKCCinhib
                    cell.pap[2] = TotCloPap_MP_70NKCCinhib
                    cell.pap[8] = TotureaPap_70NKCCinhib
                elif cell.preg == 'late':
                    cell.oi[0] = TotSodOI_LP_70NKCCinhib
                    cell.oi[1] = TotPotOI_LP_70NKCCinhib
                    cell.oi[2] = TotCloOI_LP_70NKCCinhib		
                    cell.oi[8] = TotureaOI_70NKCCinhib
                    cell.pap[0] = TotSodPap_LP_70NKCCinhib
                    cell.pap[1] = TotPotPap_LP_70NKCCinhib
                    cell.pap[2] = TotCloPap_LP_70NKCCinhib
                    cell.pap[8] = TotureaPap_70NKCCinhib
            elif cell.inhib == 'NKCC2-100':
                if cell.preg == 'non':
                    cell.oi[0] = TotSodOI_100NKCCinhib
                    cell.oi[1] = TotPotOI_100NKCCinhib
                    cell.oi[2] = TotCloOI_100NKCCinhib		
                    cell.oi[8] = TotureaOI_100NKCCinhib
                    cell.pap[0] = TotSodPap_100NKCCinhib
                    cell.pap[1] = TotPotPap_100NKCCinhib
                    cell.pap[2] = TotCloPap_100NKCCinhib
                    cell.pap[8] = TotureaPap_100NKCCinhib
                elif cell.preg == 'mid':
                    cell.oi[0] = TotSodOI_MP_100NKCCinhib
                    cell.oi[1] = TotPotOI_MP_100NKCCinhib
                    cell.oi[2] = TotCloOI_MP_100NKCCinhib		
                    cell.oi[8] = TotureaOI_100NKCCinhib
                    cell.pap[0] = TotSodPap_MP_100NKCCinhib
                    cell.pap[1] = TotPotPap_MP_100NKCCinhib
                    cell.pap[2] = TotCloPap_MP_100NKCCinhib
                    cell.pap[8] = TotureaPap_100NKCCinhib
                elif cell.preg == 'late':
                    cell.oi[0] = TotSodOI_LP_100NKCCinhib
                    cell.oi[1] = TotPotOI_LP_100NKCCinhib
                    cell.oi[2] = TotCloOI_LP_100NKCCinhib		
                    cell.oi[8] = TotureaOI_100NKCCinhib
                    cell.pap[0] = TotSodPap_LP_100NKCCinhib
                    cell.pap[1] = TotPotPap_LP_100NKCCinhib
                    cell.pap[2] = TotCloPap_LP_100NKCCinhib
                    cell.pap[8] = TotureaPap_100NKCCinhib
            elif cell.inhib == 'SNB-70':
                cell.oi[0] = TotSodOI_70NKCCinhib
                cell.oi[1] = TotPotOI_70NKCCinhib
                cell.oi[2] = TotCloOI_70NKCCinhib		
                cell.oi[8] = TotureaOI_70NKCCinhib
                cell.pap[0] = TotSodPap_70NKCCinhib
                cell.pap[1] = TotPotPap_70NKCCinhib
                cell.pap[2] = TotCloPap_70NKCCinhib
                cell.pap[8] = TotureaPap_70NKCCinhib
            elif cell.inhib == 'SNB-100':
                cell.oi[0] = TotSodOI_100NKCCinhib
                cell.oi[1] = TotPotOI_100NKCCinhib
                cell.oi[2] = TotCloOI_100NKCCinhib		
                cell.oi[8] = TotureaOI_100NKCCinhib
                cell.pap[0] = TotSodPap_100NKCCinhib
                cell.pap[1] = TotPotPap_100NKCCinhib
                cell.pap[2] = TotCloPap_100NKCCinhib
                cell.pap[8] = TotureaPap_100NKCCinhib
            elif cell.preg == 'mid':
                cell.oi[0] = TotSodOI_MP
                cell.oi[1] = TotPotOI_MP
                cell.oi[2] = TotCloOI_MP
                cell.oi[3] = TotHCO3OI_MP
                cell.pap[0] = TotSodPap_MP
                cell.pap[1] = TotPotPap_MP
                cell.pap[2] = TotCloPap_MP
                cell.pap[3] = TotHCO3Pap_MP
            elif cell.preg == 'late':
                cell.oi[0] = TotSodOI_LP
                cell.oi[1] = TotPotOI_LP
                cell.oi[2] = TotCloOI_LP
                cell.oi[3] = TotHCO3OI_LP
                cell.pap[0] = TotSodPap_LP
                cell.pap[1] = TotPotPap_LP
                cell.pap[2] = TotCloPap_LP
                cell.pap[3] = TotHCO3Pap_LP
                        
        if cell.segment=='cTAL' or cell.segment == 'MD' or cell.segment=='DCT' or cell.segment=='PT' or cell.segment == 'CNT' or cell.segment == 'CCD':
            if cell.species == 'rat':
                pHplasma = 7.323
                phpap = 7.0
            elif cell.species == 'mou':
                pHplasma = 7.323
                phpap = 7.0
            elif cell.species == 'hum':
                pHplasma = 7.4
                phpap = 7.3
            else:
                print('cell.species: ' + str(cell.species))
                raise Exception('what is species?')

            if cell.unx == 'Y':
                pHplasma = 7.4-0.1
                phpap = 7.3-0.1
            
            cell.pH[5] = pHplasma

            facbic = np.exp(np.log(10)*(cell.pH[5]-pKHCO3))
            facpho = np.exp(np.log(10)*(cell.pH[5]-pKHPO4))
            facamm = np.exp(np.log(10)*(cell.pH[5]-pKNH3))
            fachco2 = np.exp(np.log(10)*(cell.pH[5]-pKHCO2))

            cell.conc[0,5] = cell.cm[0]
            cell.conc[1,5] = cell.cm[1]
            cell.conc[2,5] = cell.cm[2]
            cell.conc[3,5] = cell.cm[3]
            cell.conc[4,5] = cell.cm[4]
            cell.conc[5,5] = cell.cm[5]
                
            cell.conc[4,5] = cell.conc[5,5]/342.07
            cell.conc[3,5] = cell.conc[4,5]*10**(cell.pH[5]-pKHCO3)
            cell.conc[6,5] = cell.cm[6]*facpho/(1+facpho)
            cell.conc[7,5] = cell.cm[6]/(1+facpho)

            cell.conc[8,5] = cell.cm[8]
            
        if  cell.segment =='PT':

            TotAmmCT = 0.1
            TotAmmCM = 1.5
            if cell.species == 'hum':
                total = cell.total*0.9
            elif cell.species == 'rat':
                total = cell.total*0.88
            elif cell.species == 'mou':
                total = cell.total*0.88
            else:
                print('cell.species: ' + str(cell.species))
                raise Exception('what is species?')
            pos = i/total


        Ammtotz = TotAmmCT+(TotAmmCM-TotAmmCT)*pos
                
        cell.conc[9,5] = Ammtotz*facamm/(1+facamm)
        cell.conc[10,5] = Ammtotz/(1+facamm)
        cell.conc[2,5] = cell.cm[2]+(Ammtotz/(1+facamm))

        cell.conc[11,5] = np.exp(-np.log(10)*cell.pH[5])*1000.0

        cell.conc[12,5] = cell.cm[12]*fachco2/(1+fachco2)
        cell.conc[13,5] = cell.cm[12]/(1+fachco2)
        cell.conc[14,5] = cell.cm[14]      

        elecS = 0
        for j in range(NS):
            elecS = elecS+zval[j]*cell.conc[j,5]

        cell.conc[2,5] = cell.conc[2,5]+elecS

        if cell.sex=='female' and cell.species == 'rat':
            female_conc(cell, i)
    else:
        return 0

    #cell.vol_init[:] = original_cell_vols[:]
    
    #cell.ep[:] = original_cell_ep[:]
    #cell.ep[:] = original_cell_ep[:]
   



        

#%%

'''

if cell.segment=='cTAL':
    TotAmmCT = 1.0
    TotAmmCM = 1.5
    pos = (1*(cell.total-i))/(1*cell.total)
elif cell.segment=='MD':
    TotAmmCT = 1.0
    TotAmmCM = 1.5
    pos = 1/200
elif cell.segment=='DCT':
    TotAmmCT = 0.1
    if cell.type !='sup':
        TotAmmCM = 2.9
    else:
        TotAmmCM = 1.0
    pos = (1.0*(cell.total-i))/(1.0*cell.total)
#elif cell.segment =='PT':
'''