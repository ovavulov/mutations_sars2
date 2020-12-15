# SARS Cov-2 mutations modelling and molecular docking

Materials for semester project in Bioinformatics Institute (fall 2020)

__Students__:   Ivanova E.,  Varchenko K. (SPb team), Akimenkova M., Shemyakina A., Vavulov O. (Moscow team)

__Supervisor__: Zolotarev A. (SPbSU), Danilov L. (Bioinformatics Institute)

## Problem & Results

### Introduction

SARS-CoV2 virus caused an epidemic with more than a million people deaths and the number grows. Predicting the evolution of a new virus can help fight disease more successfully. The studies on virus receptor binding domain (RBD) revealed the increase in affinity to angiotensin converting enzyme 2 (ACE2) compared to SARS-CoV that caused 2002-2004 SARS outbreak. Thus we attempted to simulate further evolution of RBD in this direction and suggest the class of drugs that could possibly inhibit mutant RBD (RBD-mut).
We analysed RBD-ACE2 interface in the PyMOL visualization system and identified key amino acid residues of their interaction. Using the FoldX based pipeline we went through all possible missense mutations in these codons and selected combinations of them that: a) preserve the stability of RBD; b) increase the stability of the RBD/ACE2 complex. The resulting RBD-mut variants showed increased hydrophobicity of interface.

FDA-approved molecules were docked into two selected RBD-muts using AutoDock software, and the resulting RBD-mut/ligand complexes were ranked by interaction energy. Possible intermolecular interactions in the first 100 RBD-mut/ligand complexes for each RBD-mut were manually analyzed in PyMOL. The molecules that could successfully bind the interface region by the formation of polar contacts and hydrophobic interactions were suggested as potentially effective competitive inhibitors of each RBD-mut. For RBD-mut 1 molecules with heterocyclic core and polar groups on the sides were chosen. For RBD-mut 2 molecules with larger aromatic parts were suggested instead of heterocycles. Thus, we identified the classes of substances that could act as the most effective inhibitors for each RBD-mut.

### Objective

Prediction of mutations in SARS-CoV-2 S-protein RBD-domain, that increase its affinity for ACE2, and the search for potentially effective competitive inhibitors of the resulting proteins

### Plan

1) Introduction to the PDB database and the basics of visualization of three-dimensional structures in PyMOL, search for mutation variants that increase affinity for ACE2.
2) Selection of the most stable mutant variants of RBD (RBD-mut) with increased affinity using FoldX.
3) Docking of FDA-approved drug molecules using AutoDock and search for potential competitive inhibitors of selected RBD-mut, detailed analysis of the interaction of selected ligands.

### Results

All results you can find in the __results__ folder. Using FoldX based pipeline we performed wide mutations screening on perspective positions. Full report is presented in __mutscreen_report.csv__. Afterwards, each team chose its own mutated RBD variant to keep going on study with: __6lzg_Repair_msk.pdb__, __6lzg_Repair_spb.pdb__. Finally, the second AutoDock based pipeline was applied. As a result, we got sorted by affinity list of ligands. You can see an example of it in __docking_report.csv__.

![](images/ligands.gif)


## Methods

Visualisation - PyMOL [[1]](#1), mutations screening - FoldX [[2]](#2), molecular docking - AutoDock Vina [[3]](#3).

## Requirements

...

## User Guide

...

## Refereneces
<a id="1">[1]</a> 
PyMOL
The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.

<a id="1">[2]</a> 
Schymkowitz J, Borg J, Stricher F, Nys R, Rousseau F, Serrano L. The FoldX web server: an online force field. Nucleic Acids Res. 2005;33(Web Server issue):W382-W388. doi:10.1093/nar/gki387

<a id="1">[3]</a> 
O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010) 455-461
