# denovoprediction
De Novo Protein Tertiary Structure Prediction Pipeline

# Pipeline
- Conformation Initializer - initializes the backbone chain
- Fragment Library Generator - creates the 3 and 9-mer fragment library
- Seef - computes the energy for a conformation (Rw potential)
- Conformation Sampler - computes the next conformation and finds the minimum
- Residue Mapper - Map ROBETTA residues to fragments/fragment library/residues
- PDB Mapper - Map conformations to PDB file format (CRANKITE)

# To run the example w/ Flask webserver
1. Install python 2
2. Install calRW (http://zhanglab.ccmb.med.umich.edu/RW/)
3. Install CRANKITE (https://sites.google.com/site/crankite/)
4. Run the webserver
    `python run_webserver.py`
5. Open http://localhost:5000 to launch the website

