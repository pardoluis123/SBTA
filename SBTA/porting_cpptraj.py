import numpy as np
import mdtraj as md
import os
import scipy.sparse as sp
from time import time


class Cpptraj_Porter:

    def __init__(self,file=None,topology=None) -> None:
        self.file=file if file is not None else None
        self.topology=topology if topology is not None else None
        if topology is not None:
            top=md.load_topology(topology)
            self.indexes = np.array([residue.resSeq+1 for residue in top.residues])

    def intialize_hbond_matrix(self,indexes=None)->np.array:
        """
        Returns a more succint array without atom identifiers from a Cpptraj input file and

        Parameters
        ----------
        file : str
            The path to the file from which data will be read
        
        """

        indexes = indexes if indexes is not None else self.indexes

        hbond_cpptrajarray=np.zeros(shape=(len(indexes)+1,len(indexes)+1),dtype=float)
        
        hbond_cpptrajarray[0,1:]=indexes
        hbond_cpptrajarray[1:,0]=indexes

        return hbond_cpptrajarray

    def populate_hbond_matrices(self,directory_input=None,matrix=None)->np.array:
        """
        Populates an initiated hydrogen bonding matrix with data from the cpptraj hbond
        outfile.

        Parameters
        ----------
        file : str
            The path to the file from which data will be read
        
        """
        directory_input = directory_input if directory_input is not None else self.file
        matrix=matrix if matrix is not None else self.intialize_hbond_matrix()


        frame_matrix=np.copy(matrix)
        frac_matrix=np.copy(matrix)
        test_matrix=np.copy(matrix)


        with open(directory_input,"r") as datafile:
            lines = datafile.readlines()
            n=0
            empty_list_unique=[]
            empty_list_res=[]
            empty_dict={}

            for line in lines[1:]:
                #format['#Acceptor', 'DonorH', 'Donor', 'Frames', 'Frac', 'AvgDist', 'AvgAng']
                line=line.split()

                acceptor=line[0]
                donor=line[2]
                frames=float(line[3])
                frac=float(line[4])
                
                if (acceptor,donor) not in dict:
                    empty_dict[(acceptor,donor)]=0
                else:
                    empty_dict[(acceptor,donor)]+=frames

                empty_list_res.append((acceptor,donor))
                empty_list_res.append((donor,acceptor))
                print(acceptor,donor)

                if donor != acceptor:
                    test_matrix[int(donor),int(acceptor)]+=frames
                    test_matrix[int(acceptor),int(donor)]+=frames

                acceptor=float(line[0].split('@')[0].split('_')[1])
                donor=float(line[2].split('@')[0].split('_')[1])
                empty_list_unique.append((acceptor,donor))
                empty_list_unique.append((donor,acceptor))
                
               

                print(acceptor,donor)

                
                frames=float(line[3])
                frac=float(line[4])

                if donor != acceptor:
                    frame_matrix[int(donor),int(acceptor)]+=frames
                    frame_matrix[int(acceptor),int(donor)]+=frames
                    frac_matrix[int(donor),int(acceptor)]+=frac
                    frac_matrix[int(acceptor),int(donor)]+=frac

                #move forward in loop
                n+=1
        print(f"length of residue list: {len(empty_list_res)}, length of unique list{len(empty_list_unique)}")
        print(f"first row of res: {empty_list_res[0,:]}, first row of list{empty_list_unique[0,:]}")
        print(f"last row of res: {empty_list_res[0,:]}, last row of list{empty_list_unique[0,:]}")
        os._exit(0)
        frame_and_frac=np.array([frame_matrix,frac_matrix])

        print(f"finalized arrays, frame_matrix \n{frame_matrix}\n\n frac_matrix\n{empty_dict}")
        
        return frame_and_frac

    def differences_per_atoms(self,directory_input=None)->dict:
        """
        Filler function just meant to take a look at how different the entries are when
        simplified to the residue level vrs the atomic level

        Parameters
        ----------
        file : str
            The path to the file from which data will be read
        
        """
        directory_input = directory_input if directory_input is not None else self.file

        with open(directory_input,"r") as datafile:
            lines = datafile.readlines()
            n=0
            unique_dict={}
            res_dict={}

            for line in lines[1:]:
                #format['#Acceptor', 'DonorH', 'Donor', 'Frames', 'Frac', 'AvgDist', 'AvgAng']
                line=line.split()

                acceptor=line[0]
                donor=line[2]
                frames=float(line[3])
                frac=float(line[4])
                
                if acceptor not in unique_dict:
                    unique_dict[acceptor]=0
                else:
                    unique_dict[acceptor]+=frames

                acceptor=float(line[0].split('@')[0].split('_')[1])
                donor=float(line[2].split('@')[0].split('_')[1])
                frames=float(line[3])
                frac=float(line[4])

                if acceptor not in res_dict:
                    res_dict[acceptor]=0
                else:
                    res_dict[acceptor]+=frames


                #move forward in loop
                n+=1
        print(f"length of residue list: {len(unique_dict)}, length of unique list{len(res_dict)}")
        
        os._exit(0)
        frame_and_frac=np.array([frame_matrix,frac_matrix])

    def compress_converted_array(self,name=None,array=None):
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
        array = array if array is not None else self.populate_hbond_matrices()
        name = name if name is not None else "filler"

        compressed_array=sp.csr_matrix(array)
        sp.save_npz(f"{name}_Replicate_Average.npz", compressed_array)

        return

    def save_converted_array(self,name=None,matrix=None):
        """
        Saves a complete Trajectory array (3D) as txt file

        Parameters
        ----------
        name:str
            base name for output files
        array:numpy.ndarray
            A 3D numpy array (symbolizing a trajectory)
            containing 2d arrays (symbolizing frames)
        """
        name=name if name is not None else "filler_name"
        matrix=matrix if matrix is not None else self.populate_hbond_matrices()[0]

        np.savetxt(f"{name}_testmatrix.txt",matrix)

        return

if __name__=="__main__":
    from Convenience import CCUCGU_CONCAT_HBOND,CCU_Topologies,CCUGCU_CONCAT_HBOND
    from Trajectory_Processor import Processor
    import os
    
    CCUCGU_TOP=CCU_Topologies[1]
    CCUGCU_TOP=CCU_Topologies[0]

    test_cpptraj_hbond=Cpptraj_Porter(file=CCUGCU_CONCAT_HBOND,topology=CCUGCU_TOP)
    test_processor=Processor()
    
    #frame_and_frac=test_cpptraj_hbond.populate_hbond_matrices()
    
    test_cpptraj_hbond.differences_per_atoms()


   
