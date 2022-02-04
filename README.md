# Ecal Rotation and Layer-shift Summary
## Contents

geometry_updates.pdf: Slides summarizing the changes made within ldmx-sw for Ecal's rotation and layer-shift geometries.
   
**I. v13-rotation**
LDMX software: https://github.com/LDMX-Software/ldmx-sw     
  Classifiers: v12-unrotated Ecal, v13-rotated Ecal   
  
**Contents in this repo**   
  **i. plots_validation**   
    - Folders v12-unrotated and v13-rotated: trigger-skimmed and no-trigger-skimmed plots of the unrotated and rotated geometrys, respectively   
    - Rotation v13 slides: a comparison of a few plots   
    
  **ii. cell_map**   
    - cellmodule.txt: map for v13 (rotated) cell centers on the Ecal face, for bdt training and particleNet   
    - cellModulePositions.png: a image of what cellmodule.txt looks like, drawing circles around cell centers   
    - map_cellmodule.cxx: outputs cellmodule.txt and cellModulePositions.png   
    
  **iii. fiducial_ratio**   
    - Trigger Skimmed Rotation slides: comparison of the Ecal electron fiducial ratios for the v12 vs. v13 geometries post-trigger-skim   
    - cellmodule_v12.txt and cellmodule_v13.txt: v12 and v13 cell center positions   
    - 10MV.py and PN.py: calculate the v13 fiducial ratio. This is a slight revision from the v12 fiducial scripts, with the original version provided by David Jiang (davidgjiang(at)ucsb.edu)   
    
  **iv. corner_acceptance**   
    - Corner Studies Rotation slides: Compare the number of events at the corners of the v12 and v13 goemetries. See that there are significantly more events at the corners for v13    
    - corner_studies.py: ouputs the ratio for the number of corner events to the number of events in the common region of v12 and v13 (overlapping region of v12 and v13 Ecal face), and an image for the contour lines of the rotated and unrotated goemetries, hexagons_image.png   

**II. layer-shift**   

In progress layer shift validation:
https://github.com/Hongyin-bug/ldmx-sw/tree/iss904-rotate-shift

