"""Molecular Featurization Pipeline

This module provides a class for calculating molecular fingerprints using various descriptors.

Overview:
    Molecular fingerprinting is a crucial step in the analysis of chemical compounds. This module offers
    a versatile class, FingerprintCalculator, capable of generating fingerprints using different descriptors.
    The fingerprints serve as compact representations of molecular structures, aiding in tasks such as
    similarity comparison, clustering, and machine learning model training.

Classes:
    FingerprintCalculator: A class for calculating molecular fingerprints.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs, Torsions

from molecularnetwork.utils import InvalidSMILESError


class FingerprintCalculator:
    """Class to calculate fingerprints for molecules."""

    def __init__(self, descriptor="morgan2"):
        """
        Initialize the FingerprintCalculator object with the specified descriptor.

        This class provides functionality to calculate molecular fingerprints using various descriptors.
        The default descriptor is 'morgan2', which generates Morgan fingerprints of radius 2.

        :param descriptor: The type of descriptor to use for fingerprint calculation (default: "morgan2").
                           Supported descriptors are:
                               - "atompairs": Atom pair fingerprints.
                               - "maccs": MACCS keys fingerprints.
                               - "morgan2": Morgan fingerprints with radius 2.
                               - "morgan3": Morgan fingerprints with radius 3.
                               - "rdkit": RDKit fingerprints.
                               - "topo": Topological torsion fingerprints.
        """

        self.descriptor = descriptor
        self.descriptors = {
            "atompairs": lambda m: Pairs.GetAtomPairFingerprint(m),
            "maccs": lambda m: MACCSkeys.GenMACCSKeys(m),
            "morgan2": lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048),
            "morgan3": lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, 2048),
            "rdkit": lambda m: FingerprintMols.FingerprintMol(m),
            "topo": lambda m: Torsions.GetTopologicalTorsionFingerprint(m),
        }

    def calculate_fingerprint(self, smiles):
        """
        Calculate the fingerprint for the given SMILES representation of a molecule.

        :param smiles: The SMILES representation of the molecule.
        :return: The calculated fingerprint.
        :raises InvalidSMILESError: If the input SMILES string is invalid.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fn = self.descriptors[self.descriptor]
            return fn(mol)
        raise InvalidSMILESError
