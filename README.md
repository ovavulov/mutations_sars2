# SARS Cov-2 mutations modelling and molecular docking

Materials for semester project in Bioinformatics Institute (fall 2020)

## Problem & Results

### Introduction

SARS-CoV2 virus caused an epidemic with more than a million people deaths and the number grows. Predicting the evolution of a new virus can help fight disease more successfully. The studies on virus receptor binding domain (RBD) revealed the increase in affinity to angiotensin converting enzyme 2 (ACE2) compared to SARS-CoV that caused 2002-2004 SARS outbreak. Thus we attempted to simulate further evolution of RBD in this direction and suggest the class of drugs that could possibly inhibit mutant RBD (RBD-mut).
We analysed RBD-ACE2 interface in the PyMOL visualization system and identified key amino acid residues of their interaction. Using the FoldX based pipeline we went through all possible missense mutations in these codons and selected combinations of them that: a) preserve the stability of RBD; b) increase the stability of the RBD/ACE2 complex. (The resulting RBD-mut variants showed increased hydrophobicity of interface.)
FDA-approved molecules were docked into (two) selected RBD-muts using AutoDock software, and the resulting RBD-mut/ligand complexes were ranked by interaction energy. Possible intermolecular interactions in the first 100 RBD-mut/ligand complexes (for each RBD-mut) were manually analyzed in PyMOL. The molecules that could successfully bind the interface region by the formation of polar contacts and hydrophobic interactions were suggested as potentially effective competitive inhibitors of each RBD-mut. For RBD-mut 1 molecules with heterocyclic core and polar groups on the sides were chosen. For RBD-mut 2 molecules with larger aromatic parts were suggested instead of heterocycles. (Thus, we identified the classes of substances that could act as the most effective inhibitors for each RBD-mut.)

### Objective

Prediction of mutations in SARS-CoV-2 S-protein RBD-domain, that increase its affinity for ACE2, and the search for potentially effective competitive inhibitors of the resulting proteins

### Plan

1) Introduction to the PDB database and the basics of visualization of three-dimensional structures in PyMOL, search for mutation variants that increase affinity for ACE2.
2) Selection of the most stable mutant variants of RBD (RBD-mut) with increased affinity using FoldX
Docking of FDA-approved drug molecules using AutoDock and search for potential competitive inhibitors of selected RBD-mut, detailed analysis of the interaction of 3) selected ligands.

### Results

...


## Methods

...

## Requirements

...

## User Guide

...

## Refereneces

...
