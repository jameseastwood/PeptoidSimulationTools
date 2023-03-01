# PeptoidSimulationTools

This repository contains a collection of scripts and resources for setting up and performing Molecular Dynamics Simulations of peptoids (oligomers of N-substituted glycine). Most of the tools were developed over the course of my dissertation research at New York University.

The broadest workflow supported by these tools is as follows:
1. Look up peptoid residues from the Peptoid Data Bank (databank.peptoids.org) and save their associated SMILES string representation to individual files
2. Generate capped monomer models of these residues and save them as .pdb files
3. Generate Amber prepin files of these models for use with GAFF
4. Generate Amber input files (topology and coordinate files) from some starting models containing peptoid residues.

The python script `build_models.py` depends on the rdkit package.

Scripts are also included for performing crystal structure analysis relevant to peptoid structures using the CSD python API from Cambridge Crystallographic Data Centre.
