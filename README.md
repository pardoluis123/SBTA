Systems Biology Trajectory Analysis Package (SBTA)
SBTA is a Systems Biology approach to analyzing Molecular Dynamics Trajectories. The program employs a network-theoretic framework currently structured using mainly the MDTraj trajectory analysis package.

Workflow Overview
The workflow is modular and divided into two major sections:

1. Data Generation
-Input and process molecular dynamics trajectory data into H_PRINT (Trajectory.py).
-Generate and output adjacency matrices. (Trajectory_Processor).

2. Data Processing and Visualization
-Create fingerprints of pairwise comparisons of residues in a trajectory using heatmapped adjacency matrices (fingerprint_maker.py).
-Cluster frames of individual Trajectories(Data_manip.py).
-Filter Adjacency matrices for residues of interest(Data_manip.py).

Convenient file management operations(File_Management.py), such as:
    -Creating symbolic links to files without full read/write permissions.
    -Stripping and saving new trajectories.

Data management file operations(Data_manip.py) include:
    -Filtering specific residues from adjacency matrices.
    -Stripping zero values for sparse matrices.
    -Performing matrix-wise operations.
    -Kmeans Clustering.
    -Importing CPPTRAJ Data (Hbond as of now).
    

There are two example runfiles for simplified usage of the program:

Template_trajectory.py: A "fill in the blank" Python file where you point to specific trajectories, and the program runs on its own using a streamlined operations subclass.
Template_FingerPrint.py: Another "fill in the blank" Python file that allows you to point to previously created adjacency arrays and create "whole" fingerprints or "difference" fingerprints (difference between two arrays).


Important Notes
The program requires trajectories with topological information, which must be loaded accordingly.
It supports trajectory formats with inherent topological information, such as PDB files.


Development Details

SBTA relies on several Python packages:
os
numpy
scipy
MDTraj
matplotlib
sys

Additionally, SGMD includes a convenience file (Convenience.py) that stores useful items such as:
Lists of residues.
Residue-name pairings.
Paths to specific trajectories (stored as Python lists and converted to numpy arrays).
