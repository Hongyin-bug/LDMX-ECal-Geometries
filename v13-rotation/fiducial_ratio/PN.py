import numpy as np
import glob
import uproot
import matplotlib.pyplot as plt
import concurrent.futures
import copy
import matplotlib
executor = concurrent.futures.ThreadPoolExecutor(12)

base = '/home/aminali/production/rotation_prod/v300_tskim/bkg_output/'
files = glob.glob(base+'*.root')

load_branches = ['EcalScoringPlaneHits_v3_v13.x_', 'EcalScoringPlaneHits_v3_v13.y_', 'EcalScoringPlaneHits_v3_v13.z_', 
'EcalScoringPlaneHits_v3_v13.px_', 'EcalScoringPlaneHits_v3_v13.py_', 'EcalScoringPlaneHits_v3_v13.pz_', 
'EcalScoringPlaneHits_v3_v13.pdgID_', 'EcalScoringPlaneHits_v3_v13.trackID_']

# Constants (mm)
EcalSP = 240.5005
EcalFace = 248.35
cell_radius = 5.0

# Projection Functions 
def projectionX(x,y,z,px,py,pz,zFinal):
    if (px == 0):
        return x + (zFinal - z)/99999
    else:
        return x + px/pz*(zFinal - z)

def projectionY(x,y,z,px,py,pz,zFinal):
    if (py == 0):
        return y + (zFinal - z)/99999
    else:
        return y + py/pz*(zFinal - z)

# Distance Function
def dist(p1, p2):
    return np.sqrt(np.sum( ( np.array(p1) - np.array(p2) )**2 ))

# Load the Cell Map
def loadCellMap():
    cellMap = {}
    for i, x, y in np.loadtxt('cellmodule_v13.txt'):
        cellMap[i] = (x, y)
    global cells
    cells = np.array(list(cellMap.values()))
    print("Loaded detector info")

# Function for getting the magnitudes of the recoil angles
def FiducialCounter(filelist):
      
    print("Reading files")
    
    f_events = 0 # number of fiducial events
    total_events = 0 # number of all events (fiducial and non-fiducial)
    
    file_number = 1
    for f in filelist:
        if (file_number > 20):
            print("    Finished.")
            break
        print("    File: {}".format(f))
        t = uproot.open(f)['LDMX_Events']
        if len(t.keys()) == 0:
            print("    File empty, skipping")
        table_temp = t.arrays(expressions=load_branches, interpretation_executor=executor)
        table = {}
        for k in load_branches:
            table[k] = table_temp[k]

        total_events += len(table['EcalScoringPlaneHits_v3_v13.pdgID_'])

        for event in range(len(table['EcalScoringPlaneHits_v3_v13.pdgID_'])):

            fiducial = False
                        
            for hit in range(len(table['EcalScoringPlaneHits_v3_v13.pdgID_'][event])):
                if ((table['EcalScoringPlaneHits_v3_v13.pdgID_'][event][hit] == 11) and
                   (table['EcalScoringPlaneHits_v3_v13.trackID_'][event][hit] == 1) and
                   (table['EcalScoringPlaneHits_v3_v13.z_'][event][hit] > 240.500) and
                   (table['EcalScoringPlaneHits_v3_v13.z_'][event][hit] < 240.501) and
                   (table['EcalScoringPlaneHits_v3_v13.pz_'][event][hit] > 0)): 

                    recoilX = table['EcalScoringPlaneHits_v3_v13.x_'][event][hit]
                    recoilY = table['EcalScoringPlaneHits_v3_v13.y_'][event][hit]
                    recoilZ = table['EcalScoringPlaneHits_v3_v13.z_'][event][hit]
                    recoilPx = table['EcalScoringPlaneHits_v3_v13.px_'][event][hit]
                    recoilPy = table['EcalScoringPlaneHits_v3_v13.py_'][event][hit]
                    recoilPz = table['EcalScoringPlaneHits_v3_v13.pz_'][event][hit]

                    # check if it's non-fiducial/fiducial
                    finalXY = (projectionX(recoilX,recoilY,recoilZ,recoilPx,recoilPy,recoilPz,EcalFace),projectionY(recoilX,recoilY,recoilZ,recoilPx,recoilPy,recoilPz,EcalFace))
                    if not recoilX == -9999 and not recoilY == -9999 and not recoilPx == -9999 and not recoilPy == -9999:
                        for cell in range(len(cells)):
                            celldis = dist(cells[cell], finalXY)
                            if celldis <= cell_radius:
                                fiducial = True
                                break
                    
            if fiducial == True: # filter for non-fiducial/fiducial
                f_events += 1
            
            if (event % 1000 == 0):
                print('    Finished Event ' + str(event))

        file_number += 1
        
    return f_events, total_events

if __name__ == '__main__':  
    print('--- Fiducial Fraction Program ---')
    loadCellMap() # Load Cell Map
    num_events, all_events = FiducialCounter(files)
    print()
    print("=== General Info ===")
    print("Total number of events: " + str(all_events))
    print("Number of fiducial events: " + str(num_events))
    print("Fiducial Fraction pn_v3_v13: " + str(num_events/all_events))
        

