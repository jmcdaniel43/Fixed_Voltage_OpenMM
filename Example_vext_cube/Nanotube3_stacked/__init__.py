import os

pdbfile = os.path.dirname(__file__) + "/electrode_shift.pdb"

try:
   from . import nanotube_data
except ImportError:
   pass
