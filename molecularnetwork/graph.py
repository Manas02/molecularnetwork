"""Generate Molecular Network

This module provides a class for generating molecular networks based on molecular structures
represented by SMILES strings and their associated classes. The network is constructed using
molecular fingerprints and a specified similarity metric.

Classes:
    MolecularNetwork: A class for generating molecular networks.
"""


import numpy as np
import networkx
from joblib import dump, load
from molecularnetwork.featurizer import FingerprintCalculator
from molecularnetwork.similarity import SimilarityCalculator


class MolecularNetwork:
    """Generates a molecular network graph based on SMILES strings and associated classes.

    Attributes:
        descriptor (str): The type of molecular fingerprint to be used. Defaults to "morgan2".
        sim_metric (str): The similarity metric to be used for edge formation. Defaults to "tanimoto".
        sim_threshold (float): The threshold value for similarity above which edges are added between nodes. Defaults to 0.7.
        fingerprint_calculator (FingerprintCalculator): An instance of FingerprintCalculator for calculating fingerprints.
        similarity_calculator (SimilarityCalculator): An instance of SimilarityCalculator for calculating similarities.
        graph (networkx.Graph): The molecular network graph.

    Methods:
        create_graph(smiles_list, classes): Creates the molecular network graph.
        get_network(): Returns the generated molecular network graph.
        save_graph(graph_name: str): Saves the generated graph to a file.
        read_graph(graph_name: str): Reads a graph from a file.
    """

    def __init__(self, descriptor="morgan2", sim_metric="tanimoto", sim_threshold=0.7):
        """Initializes the MolecularNetwork.

        Args:
            descriptor (str, optional): The type of molecular fingerprint to be used. Defaults to "morgan2".
            sim_metric (str, optional): The similarity metric to be used for edge formation. Defaults to "tanimoto".
            sim_threshold (float, optional): The threshold value for similarity above which edges are added between nodes. Defaults to 0.7.
        """
        self.sim_threshold = sim_threshold
        self.fingerprint_calculator = FingerprintCalculator(descriptor)
        self.similarity_calculator = SimilarityCalculator(sim_metric)
        self.graph = networkx.Graph()

    def _create_graph(self, smiles_list, classes):
        """Creates the molecular network graph based on SMILES strings and associated classes.

        Args:
            smiles_list (list): List of SMILES strings representing molecular structures.
            classes (array-like): Array-like object containing the associated classes for each molecule.
        """
        fps = self._calculate_fingerprints(smiles_list)
        unique_classes, categorical_labels = self._convert_classes(classes)
        self._add_nodes(smiles_list, unique_classes, categorical_labels)
        self._add_edges(fps)

    def _calculate_fingerprints(self, smiles_list):
        """Calculates fingerprints for given SMILES strings.

        Args:
            smiles_list (list): List of SMILES strings.

        Returns:
            list: List of calculated fingerprints.
        """
        return [
            self.fingerprint_calculator.calculate_fingerprint(smi)
            for smi in smiles_list
        ]

    @staticmethod
    def _convert_classes(classes):
        """Converts classes to categorical labels.

        Args:
            classes (array-like): Array-like object containing the associated classes for each molecule.

        Returns:
            tuple: Tuple containing unique classes and corresponding categorical labels.
        """
        unique_classes = np.unique(classes)
        categorical_labels = np.arange(len(unique_classes))
        class_labels = np.array(
            [categorical_labels[np.where(unique_classes == c)[0][0]] for c in classes]
        )
        return unique_classes, class_labels

    def _add_nodes(self, smiles_list, unique_classes, categorical_labels):
        """Adds nodes to the graph.

        Args:
            smiles_list (list): List of SMILES strings representing molecular structures.
            unique_classes (numpy.ndarray): Array of unique classes.
            categorical_labels (numpy.ndarray): Array of categorical labels.
        """
        num_nodes = len(smiles_list)
        nodes = range(num_nodes)
        weighted_nodes = [
            (
                node,
                {
                    "smiles": smiles_list[node],
                    "categorical_label": str(unique_classes[value]),
                },
            )
            for node, value in zip(nodes, categorical_labels)
        ]
        self.graph.add_nodes_from(weighted_nodes)

    def _add_edges(self, fps):
        """Adds edges to the graph based on similarity threshold.

        Args:
            fps (list): List of molecular fingerprints.
        """
        num_nodes = len(fps)
        for i in range(num_nodes):
            for j in range(i + 1, num_nodes):
                sim_val = self._calculate_similarity(fps[i], fps[j])
                if sim_val > self.sim_threshold:
                    self.graph.add_edge(i, j)

    def _calculate_similarity(self, fp1, fp2):
        """Calculates similarity between two fingerprints.

        Args:
            fp1: Molecular fingerprint of the first molecule.
            fp2: Molecular fingerprint of the second molecule.

        Returns:
            float: Similarity value between the two fingerprints.
        """
        return self.similarity_calculator.calculate_similarity(fp1, fp2)

    def create_graph(self, smiles_list, classes):
        """Creates the molecular network graph.

        Args:
            smiles_list (list): List of SMILES strings representing molecular structures.
            classes (array-like): Array-like object containing the associated classes for each molecule.

        Returns:
            networkx.Graph: The generated molecular network graph.
        """
        self._create_graph(smiles_list, classes)
        return self.graph

    def get_network(self):
        """Returns the generated molecular network graph.

        Returns:
            networkx.Graph: The molecular network graph.
        """
        return self.graph

    def save_graph(self, graph_name: str):
        """Saves the generated graph to a file using joblib.

        Args:
            graph_name (str): The name of the file to save the graph.
        """
        dump(self.graph, graph_name)

    def read_graph(self, graph_name: str):
        """Reads a graph from a joblib networkx graph object file.

        Args:
            graph_name (str): The name of the file containing the graph.

        Returns:
            networkx.Graph: The read molecular network graph.
        """
        self.graph = load(graph_name)
        return self.graph
