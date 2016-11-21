# denovoprediction
De Novo Protein Tertiary Structure Prediction Pipeline

# Pipeline
- Conformation Initializer - initializes the backbone chain
- Fragment Library Generator - creates the 3 and 9-mer fragment library
- Seef - computes the energy for a conformation
- Conformation Sampler - computes the next conformation and finds the minimum
- Residue Mapper - Map ROBETTA residues to fragments/fragment library/residues
- PDB Mapper - Map conformations to PDB file format

# To run
1) Install python 2
2) Initialize your pipeline in run.py
3) `python run.py`
