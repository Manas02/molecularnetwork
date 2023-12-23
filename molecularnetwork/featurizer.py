"""Molecular Featurization Pipeline"""

from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, FingerprintMols
from rdkit.Chem.AtomPairs import Pairs, Torsions


class MolecularDescriptors:
    """
    Class for managing molecular descriptors (featurizers).
    """

    def __init__(self):
        self.descriptors = {
            "atompairs": lambda m: Pairs.GetAtomPairFingerprint(m),
            "maccs": lambda m: MACCSkeys.GenMACCSKeys(m),
            "morgan2": lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048),
            "morgan3": lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, 2048),
            "rdkit": lambda m: FingerprintMols.FingerprintMol(m),
            "topo": lambda m: Torsions.GetTopologicalTorsionFingerprint(m),
        }

    def get_descriptor_function(self, descriptor_name):
        """
        Get the molecular descriptor function based on the descriptor name.

        Args:
            descriptor_name (str): Name of the molecular descriptor.

        Returns:
            callable: Molecular descriptor function.
        """
        return self.descriptors.get(descriptor_name)
