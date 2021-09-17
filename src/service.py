import numpy as np
import json
from rdkit import Chem
from typing import List

from bentoml import BentoService, api, artifacts
from bentoml.adapters import JsonInput
from bentoml.service.artifacts.common import JSONArtifact
from bentoml.types import JsonSerializable


_MIN_PATH_LEN = 1
_MAX_PATH_LEN = 7
_N_BITS = 2048


class Descriptor(object):

    def __init__(self):
        super().__init__()
        self.minPathLen = _MIN_PATH_LEN
        self.maxPathLen = _MAX_PATH_LEN
        self.nbits = _N_BITS

    def _clip(self, v):
        if v > 127:
            v = 127
        return v

    def _calc(self, mols):
        fingerprints = []
        for mol in mols:
            counts = Chem.RDKFingerprint(
                mol, minPath=self.minPathLen, maxPath=self.maxPathLen,
                fpSize=self.nbits)
            fingerprints += [[self._clip(c) for c in counts]]
        return np.array(fingerprints, dtype=np.int8)


@artifacts([JSONArtifact("model")])
class Service(BentoService):
    @api(input=JsonInput(), batch=True)
    def calculate(self, input: List[JsonSerializable]):
        desc = Descriptor()
        input = input[0]
        output = []
        for inp in input:
            mol = Chem.MolFromSmiles(inp["input"])
            fp = np.array(desc.calc(mol), dtype=np.int8)
            output += [{"fp": fp}]
        return [output]
