# setup.py

from setuptools import setup, find_packages

setup(
    name="molecularnetwork",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "networkx",
        "rdkit",
    ],
    author="Manas Mahale",
    author_email="manas.m.mahale@gmail.com",
    description="A package for creating molecular networks based on molecular features and similarities.",
    url="https://github.com/Manas02/molecularnetwork",
)
