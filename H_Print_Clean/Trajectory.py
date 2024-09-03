import mdtraj as md
import numpy as np


class Trajectory:
    
    def __init__(self,trajectory=None,replicates=None,topology=None,residues=None, slice_option=0):
        """
        A custom "trajectory" class which wraps a typical mdtraj trajectory object and stores it
        and provides useful features for future processing such as a "template" array
        or a dictionary matching atom indexes to residue indexes

        Parameters
        ----------
        trajectory : str
            the path to a molecular dynamics trajectory
        
        replicates : list
            A list of paths leading to molecular dynamics trajectories 
            (assumed to be replicates of the same system)

        topology : str
            Path to an mdtraj topology that 
        
        residues : list
            A list containing  "residues of interest"
        
        slice_option : int
            Determines what operation should be done on the residues of interest
            1=analyze only thoose residues 
            0=analyze everything except theese residues

        """
        self.replicates=np.array(replicates) if replicates is not None else None
        self.residues = np.array(residues) if residues is not None else None
        self.slice_option = slice_option if slice_option is not None else 0
        self.topology=topology if topology is not None else None

        if trajectory is not None:
            self.trajectory=self.load_trajectory(trajectory,topology)
            #updates to filtered trajectory if residues are provided
            if self.residues is not None:
                self.trajectory=self.filter(self.trajectory,residues,slice_option)

            #initializes array elements for processing
            self.array_of_residues,self.Trajectory_array,self.atom_to_residue=self.initialize_array_elements()

        if replicates is not None:
            base_trajectory=replicates[0] #dummy trajectory
            self.trajectory=self.load_trajectory(base_trajectory,topology)

            #updates to filtered trajectory if residues are provided
            if self.residues is not None:
                self.trajectory=self.filter(self.trajectory,residues,slice_option)

            #initializes array elements for processing
            self.array_of_residues,self.Trajectory_array,self.atom_to_residue=self.initialize_array_elements()

    def load_trajectory(self,trajectory=None,topology=None):
        """
        Returns an initialized mdtraj trajectory
        -if the user gives no topology and the trajectory does not contain inherent topological
        information the program will not run

        Parameters
        ----------
        trajectory : str
            the path to a molecular dynamics trajectory

        topology : str
            Path to an molecular dynamics topology 

        """
        trajectory=trajectory if trajectory is not None else self.trajectory
        topology=topology if topology is not None else self.topology

        if topology:
            traj = md.load(trajectory,top=topology)
        else:
            try:
                traj = md.load(trajectory)
                topology=traj.topology
                print(f"{topology}, properly loaded")
            except:
                print("Filetype does not contain topology information, please load topology alongside or load different filetype")
                os._exit(0)
        return traj

    def filter(self,trajectory,residues,slice_option): 

        """
        Returns a filtered mdtraj trajectory with only residues of interest

        Parameters
        ----------
        trajectory : mdtraj.trajectory
            an initialized mdtraj trajectory
        
        residues : list
            A list containing  "residues of interest"

        slice_option : int
            Determines what operation should be done on the residues of interest
            1=analyze only thoose residues 
            0=analyze everything except theese residues

        """

        residues=residues if residues is not None else self.residues

        residues={residue-1 for residue in residues}
        
        if (slice_option==0):
            residues_mask = {residue.index for residue in trajectory.topology.residues if residue.resSeq not in residues}
        
        if (slice_option==1):
            residues_mask = {residue.index for residue in trajectory.topology.residues if residue.resSeq in residues}

        residues_query = " or ".join([f"residue == {resnum}" for resnum in residues_mask])
        residue_indices_selection = trajectory.top.select(residues_query) 
        filtered_trajectory = trajectory.atom_slice(residue_indices_selection) 
        
        return filtered_trajectory

    def initialize_array_elements(self,topology=None):
        """
        -a template array with residue indexes on the axis and zeros elsewhere
        -an empty list to become the "Trajectory Array"
        -atom_to_residue dictionary for a trajectory
        
        Parameters
        ----------
        topology for indexing of atoms and residues

        """
        topology=topology if topology is not None else self.trajectory.topology

        Trajectory_array=[]

        atom_to_residue = {atom.index: atom.residue.resSeq for atom in topology.atoms}
        indexes = [residue.resSeq for residue in topology.residues]
        
        array_of_residues = np.zeros((len(indexes)+1, len(indexes)+1))

        array_of_residues[0,1:]=indexes
        array_of_residues[1:,0]=indexes

        return array_of_residues,Trajectory_array,atom_to_residue

if __name__ == "__main__":
    from Convenience import test_trajectory,test_topology,test_residues,test_slice_option
    import os
    
    #Recommended Usage
    test=Trajectory(trajectory=test_trajectory,topology=test_topology,residues=test_residues,slice_option=test_slice_option)

    #Test if array was initialized correctly
    print(f"axis one \n{test.array_of_residues[0,:]}\n axis two \n{test.array_of_residues[:,0]}\n")
    