import numpy as np
import sys
import os
import mdtraj as md
import pandas as pd
from time import time
from File_Management import FileManagement
from matplotlib import pyplot
import scipy.sparse as sp
from Trajectory import *
from Convenience import test_trajectory,test_topology,test_residues,test_slice_option,restrained_residue_list
import os

class Processor:
    
    def __init__(self, Trajectory=None,name=None,datatype=None):
        
        """
        Returns a filtered array 

        Parameters
        ----------
        Trajectory:SOCS.Trajectory
            Propietary object of interest 
        name:str
            Preferred base-name for all output files
        datatype:int
            Dictates wether output arrays will be saved as sparse matrices or txt files
            0=sparse 1=txt 
        """

        self.Trajectory=Trajectory if Trajectory is not None else None 
        self.name=name if name is not None else "Temp_Name"
        self.datatype=datatype if datatype is not None else 0

        if name is not None:
            parent_directory=os.getcwd()
            self.data_directory=f"{parent_directory}/H_Print_{self.name}"
        else:
            self.data_directory=self.name
        
    def Create_Hydrogen_Matrix(self,Baker_hubard,array_template,atom_to_residue=None):
        """
        Returns an intiialized hydrogen bond matrix 

        Parameters
        ----------
        Baker_hubard:mdtraj.baker_hubbard 
            Takes as input the resulting matrix of the baker_hubbard function in mdtraj 
        array_template:numpy.ndarray
            A template array with residues of interest intiialized on each axis for adjacency matrix creation
        atom_to_residue:dict
            A dictionary of atom to residue pairings so the Baker_Hubbard atom indexes can be mapped to their
            respective residues
        """
        atom_to_residue=atom_to_residue if atom_to_residue is not None else self.Trajectory.atom_to_residue
        array_template=array_template if array_template is not None else self.Trajectory.array_of_residues
        
        index=np.copy(array_template[0,1:])
        for d_i, h_i, a_i in Baker_hubard:
            
            donor_residue_seq = atom_to_residue[d_i]
            acceptor_residue_seq = atom_to_residue[a_i]
            
            donor_indices = np.where(index == donor_residue_seq)[0] 
            acceptor_indices = np.where(index == acceptor_residue_seq)[0]
            print(index)
            print(f"our donor_residue_seq:{donor_residue_seq} and our acceptor_residue_seq{acceptor_residue_seq}")
            print(f"our donor location: {donor_indices}, our acceptor location: {acceptor_indices}")
            
            if donor_residue_seq != acceptor_residue_seq:
                array_template[donor_indices,acceptor_indices]+=1
                array_template[acceptor_indices,donor_indices]+=1
        
        return array_template

    def Process_H_Frame(self,frame,template=None):
        
        """
        Returns an initiated Hydrogen bond adjacency matrix with values for a given frame

        Parameters
        ----------
        frame:mdtraj.trajectory.frame
            Takes as input a frame of a trajectory to process
        array_template:numpy.ndarray
            A template array with residues of interest intiialized on each axis for adjacency matrix creation
        """

        template=template if template is not None else self.Trajectory.array_of_residues
        
        array_template=np.copy(template)
        Baker_Hubard=md.baker_hubbard(frame)
        
        frame_array=self.Create_Hydrogen_Matrix(Baker_Hubard,array_template)
        return frame_array

    def Trajectory_Array(self,trajectory=None):
        """
        Returns An array containing a hydrogen bonding adjacency matrix for each frame in the supplied trajectory

        Parameters
        ----------
        trajectory:mdtraj.trajectory
            Takes as input the trajectory the user wishes to process
        """
        trajectory=trajectory if trajectory is not None else self.Trajectory.trajectory
        
        Trajectory_array=np.zeros(shape=(len(trajectory),),dtype=object)
        n=0

        #iterate through all frames in the trajectory and add adjacency matrix for every frame
        for i in trajectory:
            frame_array = self.Process_H_Frame(i, self.Trajectory.array_of_residues)
            Trajectory_array[n] = frame_array
            n += 1

        #quick check of all frames
        j=0
        n=3200
        if len(trajectory)>3201:
            sum_concat=np.zeros(shape=(30,),dtype=object)

            for i in range(0,29):
                sum_concat[i]=np.sum(Trajectory_array[j:n])
                j+=3200
                n+=3200

        else:
            sum_concat=None
            
        self.sum_concat=sum_concat

        #return trajectory_array
        return Trajectory_array

    def replicate_array(self,replicates=None,name=None,datatype=None):
        """
        Returns An array containing a hydrogen bonding adjacency matrix for each frame 
        of each trajectory supplied

        Parameters
        ----------
        replicates:list
            Takes as input a list of paths leading to the different replicates being tested
        name:str
            Takes as input the base_name the user would like to use for output data
        datatype:int
            Dictates wether output arrays will be saved as sparse matrices or txt files
            0=sparse 1=txt 
            
        """

        replicates = replicates if replicates is not None else self.Trajectory.replicates
        name = name if name is not None else self.name
        datatype=datatype if datatype is not None else self.datatype

        n=1

        #Set some Numpy arrays to use as our final fill in the blanks
        Final_Array=np.zeros(shape=(len(replicates),),dtype=object) #Final_Array  
        Sum_Array=np.zeros(shape=(len(replicates),),dtype=object) #Sum_Array

        System_Start=time()

        #iterate through replicates
        for i in replicates:
            Replicate_Start=time() #checking time
            print(f"Generating Data For Replicate {n}:\n{i}\n")

            #initiate Trajectory object
            current_trajectory=Trajectory(trajectory=i,topology=self.Trajectory.topology,residues=self.Trajectory.residues,slice_option=self.Trajectory.slice_option)
            #initiate Processor object
            current_processor=Processor(current_trajectory,name,datatype)
            #Pull Trajectory Array (NP Array of Frames) from Processor 
            current_trajectory_array=current_processor.Trajectory_Array()
            #Add array to a third array of trajectory arrays so this is now a 4dimensional array
            Final_Array[n-1]=current_trajectory_array
            #Count sum for T-Testing
            Sum_Array[n-1]=np.sum(current_processor.Average_Array()[1:,1:])
            Replicate_End=time()#just another time check
            n+=1
            print(f"the total time taken for processing{i} was\n: {(Replicate_Start-Replicate_End)/60} Minutes")
        
        

        System_End=time()
        print(f"the total time taken for processing was: {(System_Start-System_End)/60} Minutes")
        Final_Array=np.array(Final_Array)
        Final_Sum_Array=np.array(Sum_Array)
        return Final_Array,Final_Sum_Array

    def Average_Array(self,array=None,template=None):
        """
        Returns an "averaged array" where the value of multiple arrays at the same
        point (i,j) will be averaged across all arrays and a "average" array is returned

        Parameters
        ----------
        array:numpy.ndarray
            Takes as input a 3d numpy array containing the 2d arrays of interest
        array_template:numpy.ndarray
            A template array with residues of interest intiialized on each axis for adjacency matrix creation
       
        """
        array=array if array is not None else self.Trajectory_Array()
        template = template if template is not None else self.Trajectory.array_of_residues
       
        #Initiate 
        data_to_process = np.copy(array)
        data_to_process=np.mean(data_to_process,axis=0)
        #data_to_process=np.round(data_to_process,decimals=2)

        data_to_process[0,:]=self.Trajectory.array_of_residues[0,:] 
        data_to_process[:,0]=self.Trajectory.array_of_residues[:,0] 

        return data_to_process
    
    def Save_Average_replicates(self,name=None,array=None):
        """
        Saves the "average array" of a given 4d numpy array containing
        3D numpy arrays

        Parameters
        ----------
        name:str
            base name for output files
        array:numpy.ndarray
            A 4D numpy array (symbolizing the replicates)
            containing 3d arrays (each trajectory)
            containing 2d arrays (each frame)
       
        """
        array,summed_averages = array if array is not None else self.replicate_array()
        name = name if name is not None else self.name
        
        array=self.Average_Array(array=array)
        compressed_array=sp.csr_matrix(array)
        sp.save_npz(f"{self.data_directory}/{name}_Replicate_Average.npz", compressed_array)
        np.savetxt(f"Counts_{name}.txt",summed_averages,fmt='%.10f')

        return

    def Save_Average_Trajectory(self,name=None,array=None):
        """
        Saves the "average array" of a given 3d numpy array containing
        2D numpy arrays

        Parameters
        ----------
        name:str
            base name for output files
        array:numpy.ndarray
            A 3D numpy array (symbolizing a trajectory)
            containing 2d arrays (symbolizing frames)
       
        """
        array = array if array is not None else self.Trajectory_Array()
        name = name if name is not None else self.name
        
        averaged_array=self.Average_Array(array=array)
        compressed_array=sp.csr_matrix(averaged_array)

        sp.save_npz(f"{self.data_directory}/{name}_Trajectory_Average.npz", compressed_array)
        if self.sum_concat is not None:
            np.savetxt(f"Counts_{name}.txt",fmt='%.10f')
        
        return
    
    def Save_Full_Trajectory(self,name=None,array=None):
        """
        Saves a complete Trajectory array (3D) as sparse matrix

        Parameters
        ----------
        name:str
            base name for output files
        array:numpy.ndarray
            A 3D numpy array (symbolizing a trajectory)
            containing 2d arrays (symbolizing frames)
        """
        array = array if array is not None else self.Trajectory_Array()
        name = name if name is not None else self.name
        compressed_array=sp.csr_matrix(array)
        sp.save_npz(f"{self.data_directory}/{name}_Replicate_Average.npz", compressed_array)

        return
    
    def Save_Full_Replicate_array(self,name=None,array=None):
        """
        Saves a complete Replicate Array (3D) as sparse matrix

        Parameters
        ----------
        name:str
            base name for output files
        array:numpy.ndarray
            A 4D numpy array (symbolizing the replicates)
            containing 3d arrays (each trajectory)
            containing 2d arrays (each frame)
        """
        
        array,summation = array if array is not None else self.replicate_array()
        name = name if name is not None else self.name

        compressed_array=sp.csr_matrix(array)
        sp.save_npz(f"{self.data_directory}/{name}_Replicate_Average.npz", compressed_array)
        
        return

        
if __name__ == "__main__":

    #Example Use for a multiple replicates of a trajectory we have the following options
    # -note functions can be used without the initial trajectory object 
    # -however initializing with Trajectory object streamlines process

    #parameter examples
    name="CCUGCU_G34"
    directory_input="/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/BKUP_5JUP_N2_GCU/"
    CCUGCU_topology="/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop"
    filter_parameter="1in50"

    nombre_dos="CCUCGU_G34"
    directory_input_two="/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/BKUP_5JUP_N2_CGU/"
    CCUCGU_topology="/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/TLEAP/5JUP_N2_CGU_nowat.prmtop"

    #running only thirty replicates
    Files_one=FileManagement(directory_input=directory_input,name=name,filter_parameter=filter_parameter)
    upto30_one=Files_one.toomany_reps()

    Files_two=FileManagement(directory_input=directory_input_two,name=nombre_dos,filter_parameter=filter_parameter)
    upto30_two=Files_two.toomany_reps()

    #create trajectory object with replicate list in place of singular trajectory
    #test=Trajectory(replicates=upto30,topology=CCUGCU_topology,residues=test_residues,slice_option=test_slice_option)
    CCUGCU=Trajectory(trajectory="/zfshomes/lperez/fingerprint/H_Print/symlink/CCUGCU_G34/mdcrd_1in50_20.mdcrd",topology=CCUGCU_topology,residues=restrained_residue_list,slice_option=test_slice_option)
    CCUCGU=Trajectory(trajectory="/zfshomes/lperez/fingerprint/H_Print/symlink/CCUCGU_G34/mdcrd_1in50_20.mdcrd",topology=CCUCGU_topology,residues=restrained_residue_list,slice_option=test_slice_option)

    #intiialize subsequent processor object with Trajectory in place
    process_CCUGCU=Processor(Trajectory=CCUGCU,name="CCUGCU",datatype=0)
    process_CCUCGU=Processor(Trajectory=CCUCGU,name="CCUCGU",datatype=0)

    #save Average Replicates
    process_CCUGCU.Save_Average_Trajectory()
    process_CCUCGU.Save_Average_Trajectory()

    process_CCUGCU.Save_Full_Replicate_array(name="Full_CCUGCU_array")
    process_CCUCGU.Save_Full_Replicate_array(name="Full_CCUCGU_array")



