
************************************GENERAL USER DETAILS************************************************

Spectral Graph Theoretic Analysis Package for Molecular Dynamics Simulations


SGMD (H_PRINT) is a network theoretic program built to integrate into various molecular dynamics simulations workflows and is currently structured using trajectory
analysis package MDTraj.


At this point in development the workflow is modular and is broken into two major sections as follows:

Data Generation-

1. Input and Process molecular dynamics trajectory data into H_PRINT (Trajectory.py)
2. Generate Adjacency Matrices using (Trajectory_Processor.py)

Data Processing-

1. One could create either "whole" fingerprints of all pairwise comparisons of residues in the system or filter out specific residues of interest to be compared (fingerprint_maker.py)
2. There are other convenient "file management" section that could preform useful operations such as creating symbolic links to files without full (rw) permissions or stripping and saving new trajectories
3. There is also a data management file that proves to be useful in the way of manipulating adjacency matrices such as filtering out specific residues, stripping zero values for sparse matrices,
   and preforming operations on matrix-wise operations.

There exist two "example" runfiles for the simplified use of the program at the moment
1.Example_trajectory.py - is a "fill in the blank" style python file where you point to your specific trajectories of interest and the program runs on its own via a subclass for streamlining operations
2.FingerPrintExample.py - is a "fill in the blank" style python file where you point to your (previously created) adjacency arrays and create "whole" and if desired "difference" (difference of the two arrays) fingerprints.
************************************GENERAL USER DETAILS************************************************

************************************Important Notes************************************************
1. Its important to note the program requires trajectories with topological information and it must be loaded as such
However it does account for trajectory formats that include inherent topological information such as pdb files



************************************Important Notes************************************************

************************************DEV DETAILS************************************************

SGMD Relies on several python packages and they are listed here

1. OS
2. NUMPY
3. SCIPY
4. MDTRAJ
5. Matplot
6. sys
7. os

SGMD also contains a convenient file for holding conveninent items such as list of residues, reisude-name pairings, or even paths to certain trajectories
(all in the style of python lists that are eventually converted to numpy arrays)

This file is called
-Convenience.py
************************************DEV DETAILS************************************************

