""" Similarity Functions """

from rdkit.Chem import DataStructs


class SimilarityMetrics:
    """
    Class for managing molecular similarity metrics.
    """

    def __init__(self):
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
            "tversky": lambda m1, m2, a, b: DataStructs.TverskySimilarity(
                m1, m2, a=a, b=b
            ),
        }

    def get_metric_function(self, metric_name):
        """
        Get the similarity metric function based on the metric name.

        Args:
            metric_name (str): Name of the similarity metric.

        Returns:
            callable: Similarity metric function.
        """
        return self.metrics.get(metric_name)
