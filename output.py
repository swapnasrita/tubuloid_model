from defs import *
from values import *
import flux
import electrochemical 
import glucose
import NHE3
import ATPase
import Pendrin
import AE1
import ENaC
import NCC, NKCC
import KCC

solute = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
compart = ['Lumen','Cell','ICA','ICB','LIS','Bath']
cw=Vref*60e6

def output_segment_results(cell,sup_or_jux,Scaletorq,file_to_save,N):
    
    if sup_or_jux != "":
        sup_or_jux = '_' + sup_or_jux

    # print as pregnant model if running pregnancy case
    sex_or_preg = cell.sex
    #========================================================
    # output concentrations in Lumen and Cell
    #========================================================
    for i in range(NS):
        file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_con_of_'+solute[i]+'_in_Lumen'+sup_or_jux+'.txt','w')
        file.write(str(cell.conc[i,0])+'\n')
        file.close()
    for i in range(NS):
        file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_con_of_'+solute[i]+'_in_Cell'+sup_or_jux+'.txt','w')
        file.write(str(cell.conc[i,1])+'\n')
        file.close()
    for i in range(NS):
        file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_con_of_'+solute[i]+'_in_Bath'+sup_or_jux+'.txt','w')
        file.write(str(cell.conc[i,5])+'\n')
        file.close()

    #========================================================
    # output water volume in Lumen and Cell
    #========================================================
    file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_water_volume_in_Lumen'+sup_or_jux+'.txt','w')
    file.write(str(cell.vol[0]*cw)+'\n')
    file.close()
    file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_water_volume_in_Cell'+sup_or_jux+'.txt','w')

    file.write(str(cell.vol[1]*cw)+'\n')
    file.close()
    
    #========================================================
    # output solute flows in Lumen and Cell
    #========================================================
    for i in range(NS):
        file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_flow_of_'+solute[i]+'_in_Lumen'+sup_or_jux+'.txt','w')
        file.write(str(cell.conc[i,0]*cell.vol[0]*cw)+'\n')
        file.close()
    for i in range(NS):
        file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_flow_of_'+solute[i]+'_in_Cell'+sup_or_jux+'.txt','w')
        file.write(str(cell.conc[i,1]*cell.vol[1]*cw)+'\n')
        file.close()
        
    #========================================================
    # output osmolality in Lumen, Cell, LIS, Bath
    #========================================================
    file_lumen = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_osmolality_in_Lumen'+sup_or_jux+'.txt','w')
    file_cell = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_osmolality_in_Cell'+sup_or_jux+'.txt','w')
    file_lis = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_osmolality_in_LIS'+sup_or_jux+'.txt','w')
    file_bath = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_osmolality_in_Bath'+sup_or_jux+'.txt','w')
    
    osm_l = 0
    osm_c = 0
    osm_lis = 0
    osm_b = 0

    for i in range(NS):
            osm_l = osm_l +cell.conc[i,0] 
            osm_c = osm_c +cell.conc[i,1] 
            osm_lis = osm_lis+cell.conc[i,4] 
            osm_b = osm_b +cell.conc[i,5]

    file_lumen.write(str(osm_l)+'\n')
    file_cell.write(str(osm_c)+'\n')
    file_lis.write(str(osm_lis)+'\n')
    file_bath.write(str(osm_b)+'\n')

    file_lumen.close()
    file_cell.close()
    file_lis.close()
    file_bath.close()

    #=======================================================
    # output pH
    #=======================================================
    file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_pH_in_Lumen'+sup_or_jux+'.txt','w')
    
    file.write(str(-np.log(abs(cell.conc[11,0])/1000)/np.log(10))+'\n')
    file.close()
    
    #========================================================
    # output luminal pressure
    #========================================================
    file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_pressure_in_Lumen'+sup_or_jux+'.txt','w')
    file.write(str(cell.pres[0])+'\n')
    file.close()

    #========================================================
    # output diameter 
    #========================================================
    file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_diameter'+sup_or_jux+'.txt', 'w')
    file.write(str(cell.diam)+'\n')
    file.close()

    #========================================================
    # output length 
    #========================================================
    file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_length'+sup_or_jux+'.txt', 'w')
    file.write(str(cell.len)+'\n')
    file.close()

    #==========================================================
    # potential gradient
    #===========================================================

    file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_potential_gradient_Lumen_Cell'+sup_or_jux+'.txt','w')
    file.write(str(cell.ep[0]-cell.ep[1])+'\n')
    file.close()

    file=open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_potential_gradient_Lumen_LIS'+sup_or_jux+'.txt','w')

    file.write(str(cell.ep[0]-cell.ep[4])+'\n')
    file.close()


    #========================================================
    # output transcellular and paracelluar Na fluxes 
    #========================================================


    # cell.area[4][5] = 0.02*max(cell.vol[4]/cell.volref[4],1.0)
    # cell.area[5][4] = cell.area[4][5]
    jvol,jsol, df = flux.compute_fluxes(cell) 

    df.to_csv('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_coupled_Na'+sup_or_jux+'.txt', sep = '\t')      

    file_Na_apical = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_apical_Na'+sup_or_jux+'.txt','a')
    file_Na_apical.write(str(jsol[0,0,1])+'\n')

    file_Na_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_paracellular_Na'+sup_or_jux+'.txt','a')
    file_Na_para.write(str(jsol[0,0,4])+'\n')

    file_K_apical = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_apical_K'+sup_or_jux+'.txt','a')
    file_K_apical.write(str(jsol[1,0,1])+'\n')

    file_K_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_paracellular_K'+sup_or_jux+'.txt','a')
    file_K_para.write(str(jsol[1,0,4])+'\n')

        #========================================================
        # output transporter-mediated fluxes 
        #========================================================
    for i in range(len(cell.trans)):
        transporter_type = cell.trans[i].type
        memb_id = cell.trans[i].membrane_id
        
        if transporter_type == 'SGLT1':
            solute_id,fluxs = glucose.sglt1(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')

        elif transporter_type == 'SGLT2':
            solute_id,fluxs = glucose.sglt2(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')

        elif transporter_type == 'GLUT1':
            solute_id,fluxs=glucose.glut1(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)
            for k in range(len([solute_id])):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+sup_or_jux+'.txt','a')
                file.write(str(fluxs*Scaletorq)+'\n')

        elif transporter_type == 'GLUT2':
            solute_id,fluxs=glucose.glut2(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)
            for k in range(len([solute_id])):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+sup_or_jux+'.txt','a')
                file.write(str(fluxs*Scaletorq)+'\n')

        elif transporter_type == 'NHE3':
            solute_id,fluxs=NHE3.nhe3(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
                
        elif transporter_type == 'NaKATPase':
            solute_id,fluxs=ATPase.nakatpase(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
                
        elif transporter_type == 'HATPase':
            solute_id,fluxs=ATPase.hatpase(cell,cell.ep,memb_id,cell.trans[i].act,cell.area)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
    
        elif transporter_type == 'HKATPase':
            solute_id,fluxs=ATPase.hkatpase(cell,memb_id,cell.trans[i].act,cell.area)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
                
        elif transporter_type == 'Pendrin':
            solute_id,fluxs=Pendrin.Pendrin(cell,memb_id,cell.trans[i].act,cell.area)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
                
        elif transporter_type == 'AE1':
            solute_id,fluxs=AE1.AE1(cell,memb_id,cell.trans[i].act,cell.area)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
                
        elif transporter_type == 'ENaC':
            solute_id,fluxs=ENaC.ENaC(cell,170,memb_id,cell.trans[i].act,cell.area, jvol)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
                
        elif transporter_type == 'NCC':
            solute_id,fluxs=NCC.NCC(cell,0,memb_id,cell.trans[i].act,cell.area)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
                
        elif transporter_type == 'NKCC2A':
            solute_id,fluxs=NKCC.nkcc2(cell,memb_id,cell.trans[i].act,cell.area, 'A')
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
                
        elif transporter_type == 'NKCC2B':
            solute_id,fluxs=NKCC.nkcc2(cell,memb_id,cell.trans[i].act,cell.area, 'B')
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
                
        elif transporter_type == 'NKCC2F':
            solute_id,fluxs=NKCC.nkcc2(cell,memb_id,cell.trans[i].act,cell.area, 'F')
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')
                
        elif transporter_type == 'KCC4':
            solute_id,fluxs=KCC.kcc4(cell.conc,memb_id,cell.trans[i].act,cell.area)
            for k in range(len(solute_id)):
                file = open('./'+file_to_save+'/'+sex_or_preg+'_'+cell.species+'_'+cell.segment+'_'+transporter_type+'_'+solute[solute_id[k]]+sup_or_jux+'.txt','a')
                file.write(str(fluxs[k]*Scaletorq)+'\n')

        else:
            print('transport: ' + transporter_type)
            raise Exception('What is this?',transporter_type)
        