from .subprocesses.polymer_structure import polymer_structure
from .subprocesses.polymer_topology import polymer_topology
from .subprocesses.polymer_clean import polymer_clean

def polymer_maker(number, angle):
	polymer_structure(angle)
	polymer_topology(number)
	polymer_clean()
