import os
import sys

sys.path.append(os.path.abspath('..'))
from md_setup.mae import MAEFile

mae = MAEFile("/lambda_stor/homes/heng.ma/Research/FoldingTraj/alpha3D/A3D-0-protein/A3D-0-protein.mae")
mae.struct.save('test.pdb')