from defs import *
from values import *
import numpy as np
import pandas as pd

def compute_cotransport (cell,delmu,jsol):

    numLA = len(cell.dLA)
    flux = np.zeros(numLA)
    df = pd.DataFrame(columns = ['flux'])
    
    for i in range(numLA):
        
        sid = list(cell.dLA[i].solute_id)
        mid = cell.dLA[i].membrane_id
        
        

        this_delmu = list(delmu[i][mid[0]][mid[1]] for i in sid)
        flux = cell.area[mid[0]][mid[1]] * cell.dLA[i].perm * sum(np.array(cell.dLA[i].coef) * this_delmu)

        ind = 0
        #overexpression of some basolateral transporters in case of NaKATPase inhibition
        # if 'O' in cell.inhib:
        #     if 0 in sid and 3 in sid:
        #         if 4 in mid:
        #             flux *= 2
                    
        
        for k in sid:
            jsol[k][mid[0]][mid[1]] += cell.dLA[i].coef[ind]*flux         
            ind += 1
            if k == 0:
                key = '_'.join(list(solute_id.keys())[list(solute_id.values()).index(ids)] for ids in sid)
                df.loc[key] = cell.dLA[i].coef[ind]*flux
                
    return jsol, df
