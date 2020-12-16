# SARS Cov-2 mutations modelling and molecular docking

Materials for semester project in Bioinformatics Institute (fall 2020)

__Students__:   Ivanova E.,  Varchenko K. (SPb team), Akimenkova M., Shemyakina A., Vavulov O. (Moscow team)

__Supervisors__: Zolotarev A. (SPbSU), Danilov L. (Bioinformatics Institute)

## Abstract

SARS-CoV2 virus caused an epidemic with more than a million people deaths and the number grows. Predicting the evolution of a new virus can help fight disease more successfully. The studies on virus receptor binding domain (RBD) revealed the increase in affinity to angiotensin converting enzyme 2 (ACE2) compared to SARS-CoV that caused 2002-2004 SARS outbreak. Thus we attempted to simulate further evolution of RBD in this direction and suggest the class of drugs that could possibly inhibit mutant RBD (RBD-mut).
We analysed RBD-ACE2 interface in the PyMOL visualization system and identified key amino acid residues of their interaction. Using the FoldX based pipeline we went through all possible missense mutations in these codons and selected combinations of them that: a) preserve the stability of RBD; b) increase the stability of the RBD/ACE2 complex. The resulting RBD-mut variants showed increased hydrophobicity of interface.

FDA-approved molecules were docked into two selected RBD-muts using AutoDock software, and the resulting RBD-mut/ligand complexes were ranked by interaction energy. Possible intermolecular interactions in the first 100 RBD-mut/ligand complexes for each RBD-mut were manually analyzed in PyMOL. The molecules that could successfully bind the interface region by the formation of polar contacts and hydrophobic interactions were suggested as potentially effective competitive inhibitors of each RBD-mut. For RBD-mut 1 molecules with heterocyclic core and polar groups on the sides were chosen. For RBD-mut 2 molecules with larger aromatic parts were suggested instead of heterocycles. Thus, we identified the classes of substances that could act as the most effective inhibitors for each RBD-mut.

## Objective

Prediction of mutations in SARS-CoV-2 S-protein RBD-domain, that increase its affinity for ACE2, and the search for potentially effective competitive inhibitors of the resulting proteins

## Plan

1) Introduction to the PDB database and the basics of visualization of three-dimensional structures in PyMOL, search for mutation variants that increase affinity for ACE2.
2) Selection of the most stable mutant variants of RBD (RBD-mut) with increased affinity using FoldX.
3) Docking of FDA-approved drug molecules using AutoDock and search for potential competitive inhibitors of selected RBD-mut, detailed analysis of the interaction of selected ligands.

## Results

All results you can find in the __results__ folder. Using FoldX based pipeline we performed wide mutations screening on perspective positions. Full report is presented in __mutscreen_report.csv__. Afterwards, each team chose its own mutated RBD variant to keep going on study with: __6lzg_Repair_msk.pdb__, __6lzg_Repair_spb.pdb__. Finally, the second AutoDock based pipeline was applied. As a result, we got sorted by affinity list of ligands. You can see an example of it in __docking_report.csv__. The result of docking for the selected ligands is presented below.

![](images/ligands.gif)


## Methods

Visualization - PyMOL [[1]](#1), mutations screening - FoldX [[2]](#2), molecular docking - AutoDock Vina [[3]](#3).

## Requirements

Compatibility is guaranteed for followed python packages versions:

    numpy==1.19.4
    pandas==1.1.4
    python-dateutil==2.8.1
    pytz==2020.4
    six==1.15.0

## User Guide

### foldx_pipeline

A FoldX based pipeline for finding mutations that increase the affinity of the SARS-Cov-2 S-protein for ACE2

#### Description

1. Correction of the structure of the complex with FoldX command RepairRBD

2. Running FoldX command PositionScan on the specified positions. Pre-selection __only SNP__ based on the RBD nucleotide sequence data. Compilation of mutations sets according to the specified parameters

3. Modification of the the RBD structure with the obtained sets of mutations, exclusion of S-protein destabilizing combinations

4. Modification of the RBD+ACE2 complex structure with the obtained sets of mutations 

5. Running FoldX command AnalyseComplex on the resulting mutant complexes, as well as on the original corrected structure of the complex. Saving results

#### First launch

__cd path_to_project/foldx_pipeline__

Setting up the environment

__pip install -r requirements.txt__

Make a correction of the complex structure using the RepairPDB command, step 1:

__python Repair.py -f config.cfg__

To make the screening of mutations at the specified positions (PositionScan, BuildModel, AnalyseComplex), steps 2 - 5:

__python Screening.py -f config.cfg__

The results of the pipeline are contained in the __reports__ folder

Search parameters are set via a config file (\*.cfg) of the following form

        complex=6lzg.pdb
        rbd=RBD.pdb
        positions=NB487h,SB477h,VB503p,TB500p,FB456p,YB489p
        posRange=2:6
        mutTake=all
        nthreads=12
        repair=done
        posscan=todo
        bmRBD=todo
        bmComplex=todo
        analyseComplex=todo
        rs=19

Parameters description:

__complex__ - file name for the RBD and ACE2 complex structure (before Repair)

__rbd__ - file name for the RBD structure (after Repair)

__positions__ - list of positions in the RBD for screening. It is served in the format of the positions parameter for the PositionScan command in FoldX ([format description](http://foldxsuite.crg.eu/command/PositionScan)). The d key (24 amino acids) is not supported yet

__posrange__ {int:int} - range of target mutations size by the mutated positions number. For example, posrange=2:3 means that we will only check for mutation sets in which the remainder substitution occurred at two or three positions

__muttake__ {float, int, "input", "all"} - since the potential number of combinations of suitable mutations can be very large, we can limit the sample for final analysis using _muttake_. If a decimal fraction from 0 to 1 is supplied, the corresponding fraction of the total number of mutations is randomly selected. If a positive integer is supplied, a fixed number of mutations is randomly selected. With the "input" parameter, the user sets the parameter manually (float or int) during the pipeline operation. When it is set to "all", all combinations of mutations are analyzed

__nthreads__ {int} - number of parallel pipeline execution threads

__repair__ {"todo", "done"} - flag for performing structure correction. The value "todo" blocks screening

__posscan__, __bmrbd__, __bmcomplex__, __analysecomplex__ {"todo", "done", "skip"} - status flags of the corresponding pipeline stages for separate execution of stages. If you want to perform a certain stage separately, it is necessary that the previous stages have the "done" flag set, and the subsequent ones - "skip" (skip the stage)

__rs__ {int} - parameter for reproducibility of the intermediate results distribution in the project catalog, does not affect the result of the pipeline

### virtual_screening

An AutoDock Vina based pipeline for finding the most affine conformation for the ligands from the FDA-approved database for the determined protein binding site

#### Description



#### First launch




## Refereneces
<a id="1">[1]</a> 
PyMOL
The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.

<a id="1">[2]</a> 
Schymkowitz J, Borg J, Stricher F, Nys R, Rousseau F, Serrano L. The FoldX web server: an online force field. Nucleic Acids Res. 2005;33(Web Server issue):W382-W388. doi:10.1093/nar/gki387

<a id="1">[3]</a> 
O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010) 455-461
