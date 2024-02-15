"""Utils for molecularnetwork

This module provides utility functions and classes for the molecularnetwork module.

Classes:
    InvalidSMILESError: An exception raised when an invalid SMILES string is encountered.
"""


class InvalidSMILESError(Exception):
    """
    Exception raised when an invalid SMILES string is encountered.

    Attributes:

        message (str): Explanation of the error.
    """
