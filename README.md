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

Visualization - PyMOL [[1]](#1), AutoDockTools [[2]](#2), mutations screening - FoldX [[3]](#3), molecular docking - AutoDock Vina [[4]](#4).

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

1. Config files and task shell script preparation according to given options and config template

2. Parallel virtual screening with AutoDock Vina

3. Get the most affine ligands. Saving results

#### First launch

__cd path_to_project/virtual_screening__

Configs and script preparation, step 1:

__python prepare_configs.py -db {folder with ligands library} -threads {parallel threads number}__

That script provides config files for every ligand docking task according with __config_template.txt__ file of following content

        receptor  = ../output/RBD_mut_protein.pdbqt
        ligand    = ../FOLDER/LIGAND.pdbqt
        center_x = -32.799
        center_y = 24.877
        center_z = 5.741
        size_x = 10
        size_y = 12
        size_z = 8
        out       = ../output/LIGAND.docked.pdbqt
        log       = ../output/LIGAND.docked.log
        num_modes = 10
        cpu = 4
        seed = 19 

Config template description:

__receptor__ - relative path to the receptor pdbqt-structure from the folder containing vina binary file

__ligand__ - similar to _receptor_ relative path to ligand directory; FOLDER and LIGAND are placeholders and replaced automatically

__center\_{x, y, z}, size\_{x, y, z}__ - docking box center coordinates and width; can be found from AutoDockTools visualization software

__out, log__ - template for the output structure and log files

__num_modes__ - number of single ligand conformations to analyse

__cpu__ - CPU number given to a single vina task

__seed__ - reproducibility factor

To take up step 2 just launch new generated shell script: __start\_screening\_\*.sh__

The final results extraction may be performed by __get_top.py__ script:

__python get_top.py -n {number of best structures to save} -from {docked ligands directory} -to {output folder}__

To renew the project you can use __wash\_project.py__ script.


## Refereneces
<a id="1">[1]</a> 
PyMOL
The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.

<a id="2">[2]</a> 
Morris, G. M., Huey, R., Lindstrom, W., Sanner, M. F., Belew, R. K., Goodsell, D. S. and Olson, A. J. (2009) Autodock4 and AutoDockTools4: automated docking with selective receptor flexiblity. J. Computational Chemistry 2009, 16: 2785-91.

<a id="3">[3]</a> 
Schymkowitz J, Borg J, Stricher F, Nys R, Rousseau F, Serrano L. The FoldX web server: an online force field. Nucleic Acids Res. 2005;33(Web Server issue):W382-W388. doi:10.1093/nar/gki387

<a id="4">[4]</a> 
O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, Journal of Computational Chemistry 31 (2010) 455-461
