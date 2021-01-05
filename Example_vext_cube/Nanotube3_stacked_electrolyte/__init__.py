import os

pdbfile = os.path.dirname(__file__) + "/equil_start.pdb"

try:
   from . import nanotube_data
except ImportError:
   pass
