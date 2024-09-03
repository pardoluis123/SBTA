import numpy as np
import scipy.sparse as sp
import mdtraj as md
from Convenience import n2_residue_numbers, test_topology,restrained_residue_list,test_residues
import os

class Data_Manipulator:
    
    def __init__(self,array_one=None,array_two=None,res_filtered=None,residues_to_filter=None,topology=None,threshold=None):
        """
        A Class for manipulating arrays as one sees fit akin to "data management" functions in 
        data software.

        Parameters
        ----------
        array_one : arraylike
            first array of interest (two slots are provided incase one wants to preform operations on the two arrays)
        
        array_Two : arraylike
            first array of interest (two slots are provided incase one wants to preform operations on the two arrays)
    
        res_filtered : str
            List of residues filtered from trajectory previously (soon to be edited out but for now important to function)
        
        residues_to_filter : str
            List of residues IN *CURRENT ARRAY* that one wishes to delete
        
        threshold : float
            An optional value for setting threshholds of significant "average" interactions
            between residues in array. 


        
        """
        self.array_one=self.process_input(array_one) if array_one is not None else None
        self.array_two=self.process_input(array_two) if array_two is not None else array_one
        self.topology=topology if topology is not None else None
        self.residues_to_filter=residues_to_filter if residues_to_filter is not None else None
        self.res_filtered=res_filtered if res_filtered is not None else None
        self.threshold = threshold if threshold is not None else None

        if topology is not None:
            self.res_filtered=[residue.resSeq for residue in topology.residues if residue.resSeq+1 not in self.res_filtered]

            if residues_to_filter is n2_residue_numbers: #checking if its a dictionry with is
                self.namedict=n2_residue_numbers
                self.residues_to_filter = list(n2_residue_numbers.keys())
                self.residues_to_filter = [residue-1 for residue in self.residues_to_filter]
            else:
                self.residues_to_filter=[residue.resSeq for residue in topology.residues if residue.resSeq+1 in residues_to_filter]
        else:
            self.residues=None
        
    def create_difference_array(self,array_one=None,array_two=None):
        """
        Returns a difference array

        Parameters
        ----------
        Array One:numpy.ndarray
            first array of interest
        Array Two:numpy.ndarray
            second array of interest         
        """
        

        array_one=array_one if array_one is not None else self.array_one
        array_two=array_two if array_two is not None else self.array_two
        
        difference_array=np.copy(array_two)
        difference_array[1:,1:]=array_two[1:,1:]-array_one[1:,1:]


        return difference_array
    
    def filter_all_over_diff(self,array=None,threshhold=None):
        """
        Returns a filtered array

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest
        threshhold:float
            float value at which the array will be filtered for "significant" values     
        """

        array = array if array is not None else self.process_input(array=self.array_one)

        if threshhold is not None:
            threshhold=self.threshold
            upperboundary=self.threshold
            lowerboundary=self.threshold*-1

        else:
            midpoint=np.median(array[1:,1:])
            max=np.max(array[1:,1:])
            min=np.min(array[1:,1:])
            upperboundary=midpoint+(.5*max)
            lowerboundary=midpoint-(-.5*min)


        print(f"midpoint,{midpoint},max{max},min{min},upperboundary{upperboundary},lowerboundary{lowerboundary}")


        row_mask = np.any((array[:, 1:] > upperboundary) | (array[:, 1:] < lowerboundary), axis=1)
        col_mask = np.any((array[1:, :] > upperboundary) | (array[1:, :] < lowerboundary), axis=0)        


        filtered_array=array[row_mask][:,col_mask]
        return filtered_array
            
    @staticmethod
    def process_input(array):
        """
        Returns a Numpy array from processed sparse matrix

        Parameters
        ----------
        Array:scipy.csr_matrix
            Array of interest
        """
        
        #array
        array=array if array is not None else None

        #filter an array
        processed_array=np.array(array.toarray())#239 long
        

        return processed_array

    def filter_array(self,array=None,residues_to_filter=None,res_filtered=None):
        """
        Returns a filtered array 

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest
        residues_to_filter:list
            A list of residues that one wants to filter out of the matrix being inputed
        res_filtered:list
            list of residues filtered from the original trajectory in step 1
        """

        #Residues to filter
        residues_to_filter=residues_to_filter if residues_to_filter is not None else self.residues_to_filter
        array=array if array is not None else self.process_input(array=self.array_one)
        res_filtered= res_filtered if res_filtered is not None else self.res_filtered
        
        #Just for correcting numbering
        if residues_to_filter=="ALL":
            filtered_array=array
            return filtered_array
        

        residues_to_filter = [0]+residues_to_filter 
        
        # Create a mask that marks the rows and columns to keep
        row_mask = np.isin(array[:, 0], residues_to_filter)

        col_mask = np.isin(array[0, :], residues_to_filter)
        
        filtered_rows=array[row_mask,:]
        filtered_array=filtered_rows[:,col_mask]

        residues_to_filter.pop()

        return filtered_array

     #filter an array
     
if __name__=="__main__":
    #load in a a sparse matrix
    
    array_CCUGCU=sp.load_npz("testtxt/CCUGCU_Replicate_Average.npz")
    array_CCUCGU=sp.load_npz("testtxt/CCUCGU_Replicate_Average.npz")
    
    #Class Initiated
    Test=Data_Manipulator(array_one=array_CCUCGU,residues_to_filter=n2_residue_numbers,res_filtered=restrained_residue_list,topology=md.load_topology("/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/TLEAP/5JUP_N2_CGU_nowat.prmtop"))
    degrees=Test.create_difference_array()

    
    


       
        
