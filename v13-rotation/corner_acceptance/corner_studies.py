from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import numpy as np
import matplotlib.pyplot as plt

import glob
import uproot
import math
import concurrent.futures
import copy
import matplotlib
executor = concurrent.futures.ThreadPoolExecutor(12)

base = '/home/aminali/production/rotation_prod/v300_tskim/sig_output/'
MeV10_files = glob.glob(base+'*0.01*.root')

load_branches = ['EcalRecHits_v3_v13.xpos_','EcalRecHits_v3_v13.ypos_']

r_= 85.0
R_= r_*2./np.sqrt(3)
gap_ = 1.5

#v12 unrotated
center_p = []
for i in range(1,7):
    x = (2. * r_ + gap_) * np.sin((i-1)*(np.pi/3.))
    y = (2. * r_ + gap_) * np.cos((i-1)*(np.pi/3.))
    center_p.append([x,y])

#v13 rotated
center_r = []
for i in range(1,7):
    x = (2. * r_ + gap_) * np.cos((i-1)*(np.pi/3.))
    y = - (2. * r_ + gap_) * np.sin((i-1)*(np.pi/3.))
    center_r.append([x,y])

p1 = Point(center_p[0][0]-R_/2., center_p[0][1]+r_)
p2 = Point(center_p[0][0]+R_/2., center_p[0][1]+r_)
p3 = Point(center_p[0][0]+R_, center_p[0][1])
p4 = Point(center_p[1][0]+R_/2., center_p[1][1]+r_)
p5 = Point(center_p[1][0]+R_, center_p[1][1])
p6 = Point(center_p[1][0]+R_/2., center_p[1][1]-r_)
p7 = Point(center_p[2][0]+R_, center_p[2][1])
p8 = Point(center_p[2][0]+R_/2., center_p[2][1]-r_)
p9 = Point(center_p[2][0]-R_/2., center_p[2][1]-r_)
p10 = Point(center_p[3][0]+R_/2., center_p[3][1]-r_)
p11 = Point(center_p[3][0]-R_/2., center_p[3][1]-r_)
p12 = Point(center_p[3][0]-R_, center_p[3][1])
p13 = Point(center_p[4][0]-R_/2., center_p[4][1]-r_)
p14 = Point(center_p[4][0]-R_, center_p[4][1])
p15 = Point(center_p[4][0]-R_/2., center_p[4][1]+r_)
p16 = Point(center_p[5][0]-R_, center_p[5][1])
p17 = Point(center_p[5][0]-R_/2., center_p[5][1]+r_)
p18 = Point(center_p[5][0]+R_/2., center_p[5][1]+r_)

r1 = Point(center_r[0][0]+r_,center_r[0][1]+R_/2.)
r2 = Point(center_r[0][0]+r_,center_r[0][1]-R_/2.)
r3 = Point(center_r[0][0],center_r[0][1]-R_)
r4 = Point(center_r[1][0]+r_,center_r[1][1]-R_/2.)
r5 = Point(center_r[1][0],center_r[1][1]-R_)
r6 = Point(center_r[1][0]-r_,center_r[1][1]-R_/2.)
r7 = Point(center_r[2][0],center_r[2][1]-R_)
r8 = Point(center_r[2][0]-r_,center_r[2][1]-R_/2.)
r9 = Point(center_r[2][0]-r_,center_r[2][1]+R_/2.)
r10 = Point(center_r[3][0]-r_,center_r[3][1]-R_/2.)
r11 = Point(center_r[3][0]-r_,center_r[3][1]+R_/2.)
r12 = Point(center_r[3][0],center_r[3][1]+R_)
r13 = Point(center_r[4][0]-r_,center_r[4][1]+R_/2.)
r14 = Point(center_r[4][0],center_r[4][1]+R_)
r15 = Point(center_r[4][0]+r_,center_r[4][1]+R_/2.)
r16 = Point(center_r[5][0],center_r[5][1]+R_)
r17 = Point(center_r[5][0]+r_,center_r[5][1]+R_/2.)
r18 = Point(center_r[5][0]+r_,center_r[5][1]-R_/2.)

unrotated = Polygon([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18])
rotated = Polygon([r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18])

#point = Point(center_r[0][0]+2*R_,center_r[0][1])

def CornerContainment(filelist):
    print("Reading files ... ")
    unrotated_total_hits = 0.
    rotated_total_hits = 0.
    both_total_hits = 0.
    neither_total_hits = 0.
    corner_rotated_total_ratio = 0.
    corner_unrotated_total_ratio = 0.
    corner_rotated_ratio_total = 0.
    corner_unrotated_ratio_total = 0.
    
    for f in filelist:
        unrotated_hits = 0.
        rotated_hits = 0.
        both_hits = 0.
        neither_hits = 0.
        print("    Loading File: {}".format(f))
        t = uproot.open(f)['LDMX_Events']
        if len(t.keys()) == 0:
            print("    File empty, skipping")
        table_temp = t.arrays(expressions=load_branches, interpretation_executor=executor)
        table = {}
        for k in load_branches:
            table[k] = table_temp[k]
            
        for event in range(len(table['EcalRecHits_v3_v13.xpos_'])):
            for hit in range(len(table['EcalRecHits_v3_v13.xpos_'][event])):
                x_pos = table['EcalRecHits_v3_v13.xpos_'][event][hit]
                y_pos = table['EcalRecHits_v3_v13.ypos_'][event][hit]
                point = Point(x_pos, y_pos)
                if unrotated.contains(point) or unrotated.touches(point):
                    unrotated_hits += 1.
                if rotated.contains(point) or rotated.touches(point):
                    rotated_hits += 1.
                if (unrotated.contains(point) or unrotated.touches(point)) and (rotated.contains(point) or rotated.touches(point)):
                    both_hits += 1.
                if (unrotated.contains(point) == False) and (unrotated.touches(point) == False) and (rotated.contains(point) == False) and (rotated.touches(point) == False):
                    neither_hits += 1.

        corner_rotated_ratio = (rotated_hits - both_hits)/both_hits
        corner_unrotated_ratio = (unrotated_hits - both_hits)/both_hits
        print("             rotated: {},   unrotated: {},   intersection: {},   neither: {} \n,             Rotated corner to intersection ratio: {},   Unrotated corner to intersection ratio: {}".format(rotated_hits, unrotated_hits, both_hits, neither_hits, corner_rotated_ratio, corner_unrotated_ratio))
        
        unrotated_total_hits += unrotated_hits
        rotated_total_hits += rotated_hits
        both_total_hits += both_hits
        neither_total_hits += neither_hits
        corner_rotated_ratio_total += corner_rotated_ratio
        corner_unrotated_ratio_total += corner_unrotated_ratio
        
    corner_rotated_averaged_ratio = corner_rotated_ratio_total/float(len(filelist))
    corner_unrotated_averaged_ratio = corner_unrotated_ratio_total/float(len(filelist))
    corner_rotated_total_ratio = (rotated_total_hits - both_total_hits)/both_total_hits
    corner_unrotated_total_ratio = (unrotated_total_hits - both_total_hits)/both_total_hits
    return rotated_total_hits, unrotated_total_hits, both_total_hits, neither_total_hits, corner_rotated_averaged_ratio, corner_unrotated_averaged_ratio, corner_rotated_total_ratio, corner_unrotated_total_ratio

#x,y = unrotated.exterior.xy
#plt.figure(1)
#plt.plot(x,y)

#xr,yr = rotated.exterior.xy
#plt.figure(2)
#plt.plot(xr,yr)

#plt.show()

if __name__ == '__main__':
    print('----------- Corner analysis of rotated vs unrotated: hit containment --------------')
    RotatedHit = CornerContainment(MeV10_files)[0]
    UnrotatedHit = CornerContainment(MeV10_files)[1]
    Intersection = CornerContainment(MeV10_files)[2]
    Neither = CornerContainment(MeV10_files)[3]
    RotatedCornerRatioAvg = CornerContainment(MeV10_files)[4]
    UnrotatedCornerRatioAvg = CornerContainment(MeV10_files)[5]
    RotatedCornerRatio = CornerContainment(MeV10_files)[6]
    UnrotatedCornerRatio = CornerContainment(MeV10_files)[7]
    print()
    print("============ 10 Mev: In summary... ============")
    print("Total Hits in v13 rotated area is: {} \n Total hits in v12 unrotated area is: {} \n Total hits in the intersection between the two geometries is: {} \n Total hits in neither areas is: {} \n Averaged ratio of rotated corner hits to intersection hits is: {} \n Averaged ratio of unrotated corner hits to intersection hits is: {} \n Total rotated corner hits to total intersection hits ratio is: {} \n Total unrotated corner hits to total intersection hits ratio is: {}".format(RotatedHit, UnrotatedHit, Intersection, Neither, RotatedCornerRatioAvg, UnrotatedCornerRatioAvg, RotatedCornerRatio, UnrotatedCornerRatio))
