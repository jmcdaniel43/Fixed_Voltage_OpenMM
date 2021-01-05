import os

pdbfile = os.path.dirname(__file__) + "/electrode_only_move_nanotube_1.5A.pdb"

try:
   from . import nanotube_data
except ImportError:
   pass
