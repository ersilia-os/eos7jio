# imports
import os
import csv
import sys

from rdkit import Chem
import numpy as np

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

_MIN_PATH_LEN = 1
_MAX_PATH_LEN = 7
_N_BITS = 2048

class Descriptor:

    def __init__(self):
        self.minPathLen = _MIN_PATH_LEN
        self.maxPathLen = _MAX_PATH_LEN
        self.nbits = _N_BITS

    def _clip(self, v):
        if v > 127:
            v = 127
        return v

    def calc(self, mol):
        counts = Chem.RDKFingerprint(
            mol, minPath=self.minPathLen, maxPath=self.maxPathLen,
            fpSize=self.nbits)
        return [self._clip(c) for c in counts]

desc = Descriptor()

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# run model

outputs = []
empty_output = [None]*_N_BITS
for smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        outputs += [empty_output]
    else:
        fp = desc.calc(mol)
        fp = np.array(fp, dtype=int)
        outputs += [fp]

# check input and output have the same length
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

header = ["dimension_{0}".format(str(i).zfill(4)) for i in range(_N_BITS)]

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(header)  # header
    for o in outputs:
        writer.writerow(o)
