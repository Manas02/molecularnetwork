"""Similarity Calculator

This class provides functionality to calculate molecular similarity using various metrics.

Attributes:
    sim_metric (str): The similarity metric to be used for calculation.
    metrics (dict): A dictionary mapping similarity metric names to corresponding functions.

Methods:
    __init__(self, sim_metric="tanimoto"): Initializes the SimilarityCalculator with a specified similarity metric.
    calculate_similarity(fp1, fp2): Calculates the similarity between two molecular fingerprints.
"""

from rdkit.Chem import DataStructs


class SimilarityCalculator:
    def __init__(self, sim_metric="tanimoto"):
        """Initialize the SimilarityCalculator.

        Args:
            sim_metric (str, optional): The similarity metric to be used for calculation.
                Defaults to "tanimoto". Supported metrics are:
                - "asymmetric": AsymmetricSimilarity
                - "braunblanquet": BraunBlanquetSimilarity
                - "cosine": CosineSimilarity
                - "dice": DiceSimilarity
                - "kulczynski": KulczynskiSimilarity
                - "mcconnaughey": McConnaugheySimilarity
                - "onbit": OnBitSimilarity
                - "rogotgoldberg": RogotGoldbergSimilarity
                - "russel": RusselSimilarity
                - "sokal": SokalSimilarity
                - "tanimoto": TanimotoSimilarity
                - "tversky": TverskySimilarity with parameters a=m1, b=m2
        """
        self.sim_metric = sim_metric
        self.metrics = {
            "asymmetric": DataStructs.AsymmetricSimilarity,
            "braunblanquet": DataStructs.BraunBlanquetSimilarity,
            "cosine": DataStructs.CosineSimilarity,
            "dice": DataStructs.DiceSimilarity,
            "kulczynski": DataStructs.KulczynskiSimilarity,
            "mcconnaughey": DataStructs.McConnaugheySimilarity,
            "onbit": DataStructs.OnBitSimilarity,
            "rogotgoldberg": DataStructs.RogotGoldbergSimilarity,
            "russel": DataStructs.RusselSimilarity,
            "sokal": DataStructs.SokalSimilarity,
            "tanimoto": DataStructs.TanimotoSimilarity,
            "tversky": lambda m1, m2, alpha, beta: DataStructs.TverskySimilarity(
                m1, m2, a=alpha, b=beta
            ),
        }

    def calculate_similarity(self, fp1, fp2):
        """Calculate the similarity between two molecular fingerprints.

        Args:
            fp1: The first molecular fingerprint.
            fp2: The second molecular fingerprint.

        Returns:
            float: The calculated similarity score.
        """
        return max(
            self.metrics[self.sim_metric](fp1, fp2),
            self.metrics[self.sim_metric](fp2, fp1),
        )
