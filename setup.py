# setup.py

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name="molecularnetwork",
    version="0.3.2",
    packages=find_packages(),
    install_requires=[
        "setuptools~=68.2.2",
        "numpy~=1.26.3",
        "networkx~=3.1",
        "joblib~=1.2.0",
        "rdkit~=2023.9.5"
    ],
    extras_require={"test": ["pandas~=2.1.4"]},
    author="Manas Mahale",
    author_email="manas.m.mahale@gmail.com",
    description="A package for creating molecular networks based on molecular features and similarities.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Manas02/molecularnetwork",
)
