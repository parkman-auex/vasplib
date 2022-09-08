#!/Users/parkman/Documents/TOOLs/miniconda2/envs/my_pymatgen/bin/python
import numpy as np
import sys
sys.path.append('/Users/parkman/Desktop/晶体场劈裂/bin/PyCrystalField')
import PyCrystalField as cef

########### Import CIF file
CRSYTAL_Lig, ion = cef.importCIF('NaErS2.cif','Er')
########### print eigenvectors

ion.printLaTexEigenvectors()
