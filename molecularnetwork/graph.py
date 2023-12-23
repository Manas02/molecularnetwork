"""Generate Molecular Network"""


import numpy as np
import networkx
from rdkit import Chem
from pyvis.network import Network
from .utils import InvalidSMILESError

from .featurizer import MolecularDescriptors
from .similarity import SimilarityMetrics


class MolecularNetwork:
    """
    Class for creating a molecular network based on molecular similarities.
    """

    def __init__(self, descriptor="morgan2", sim_metric="tanimoto", sim_threshold=0.7):
        """
        Initialize the MolecularNetwork object.

        Args:
            descriptor (str): Molecular descriptor type.
            sim_metric (str): Similarity metric.
            sim_threshold (float): Similarity threshold for creating edges in the network.
        """
        self.sim_threshold = sim_threshold
        self.descriptor = descriptor
        self.descriptors_manager = MolecularDescriptors()
        self.sim_metric = sim_metric
        self.metrics_manager = SimilarityMetrics()
        self.graph = networkx.Graph()

        self._validate_function(
            self.descriptor,
            self.descriptors_manager.descriptors,
            "molecular descriptor",
        )
        self._validate_function(
            self.sim_metric, self.metrics_manager.metrics, "similarity metric"
        )

    def _validate_function(self, name, function_dict, function_type):
        if name not in function_dict:
            raise ValueError(f"Unsupported {function_type}: {name}")

    def calculate_fp(self, smi: str) -> np.ndarray:
        """
        Calculate the molecular fingerprint for a given SMILES string.

        Args:
            smi (str): SMILES string of the molecule.

        Returns:
            numpy.ndarray: Molecular fingerprint.
        """
        mol = Chem.MolFromSmiles(smi)
        if mol:
            return self.descriptors_manager.descriptors[self.descriptor](mol)
        raise InvalidSMILESError

    def _create_graph(self, smiles_list: list, classes: list):
        """
        Create a molecular network graph based on molecular similarities.

        Args:
            smiles_list (list): List of SMILES strings of the molecules.
            classes (list): List of class labels for the molecules.
        """
        fps = [self.calculate_fp(smi) for smi in smiles_list]
        unique_classes, categorical_labels = np.unique(classes, return_inverse=True)
        self.graph.add_nodes_from(
            enumerate(smiles_list, categorical_label=str(unique_classes[label]))
            for label in categorical_labels
        )
        self.graph.add_edges_from(
            (i, j)
            for i in range(len(fps))
            for j in range(i + 1, len(fps))
            if self._calculate_similarity(fps[i], fps[j]) > self.sim_threshold
        )

    def _calculate_similarity(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """
        Calculate similarity between two fingerprints.

        Args:
            fp1 (numpy.ndarray): Fingerprint of molecule 1.
            fp2 (numpy.ndarray): Fingerprint of molecule 2.

        Returns:
            float: Similarity value.
        """
        return max(
            self.metrics_manager.metrics[self.sim_metric](fp1, fp2),
            self.metrics_manager.metrics[self.sim_metric](fp2, fp1),
        )

    def create_graph(self, smiles_list: list, classes: list) -> networkx.Graph:
        """
        Create the molecular network graph based on the given SMILES list and classes.

        Args:
            smiles_list (list): List of SMILES strings of the molecules.
            classes (list): List of class labels for the molecules.

        Returns:
            networkx.Graph: Molecular network graph.
        """
        self._create_graph(smiles_list, classes)
        return self.graph

    def get_network(self) -> networkx.Graph:
        """
        Get the molecular network graph.

        Returns:
            networkx.Graph: Molecular network graph.
        """
        return self.graph

    def save_graph(self, file_path: str):
        """
        Save the molecular network graph to a file.

        Args:
            file_path (str): Path to the file.

        Raises:
            ValueError: If the file format is not supported.
        """
        supported_formats = [
            ".gml",
            ".graphml",
            ".gexf",
            ".gpickle",
            ".yaml",
            ".yml",
            ".json",
            ".adjlist",
        ]
        file_format = file_path.split(".")[-1].lower()

        if file_format not in supported_formats:
            raise ValueError(f"Unsupported file format: {file_format}")

        networkx.write_gml(self.graph, file_path)
