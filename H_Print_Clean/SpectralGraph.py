import numpy as np
import scipy.sparse as sp
import mdtraj as md
from Convenience import n2_residue_numbers, test_topology,restrained_residue_list,test_residues
import os
from Data_Manip import Data_Manipulator 
process_input=Data_Manipulator.process_input

class Spectral_Analysis:

    def __init__(self,array=None):
        self.array=process_input(array=array) if array is not None else None  
        self.degreematrix=self.create_degree_matrix()
        self.neighbors=self.calculate_neighbors()
        return

    #Create degree matrix
    def create_degree_matrix(self,array=None):
        
        array=array if array is not None else self.array
        degrees=np.zeros(shape=(len(array[0,:]),len(array[0,:])),dtype=object)
        degrees[0,:]=array[0,:]
        degrees[:,0]=array[:,0]

        for i in range(1,len(array[0,:])):
            current_row=array[i,1:]
            count=np.count_nonzero(current_row!=0)
            degrees[i,i]=count

        return degrees

    # Calculare clustering coefficients
    def calculate_neighbors(self,array=None):
        array=array if array is not None else self.array

        neighbors={}
        
        for i in range(1,array.shape[0]):
            current_row=array[i,1:]
            count=np.where(current_row>0)[0] + 1#adjust for skipped row
            neighbors[i]=count

        return neighbors
    
    # Calculate clustering coefficients
    def calculate_clustering_coefficients(self,neighbors=None,array=None,degree=None):

        array=array if array is not None else self.array
        neighbors=neighbors if neighbors is not None else self.calculate_neighbors()
        degree=degree if degree is not None else self.degreematrix

        clustering_coefficients={}
        average_coefficient=0

        for index,current_neighbors in neighbors.items():
            if len(current_neighbors)>1:
                array_to_check = array[np.ix_(current_neighbors, current_neighbors)]
                number_of_edges=np.count_nonzero(array_to_check)
                number_of_edges=number_of_edges/2
                top_fraction=2*number_of_edges
                bottom_fraction=degree[index,index]*(degree[index,index]-1)
                corellation_coefficient=top_fraction/bottom_fraction
                clustering_coefficients[index]=corellation_coefficient
                average_coefficient+=corellation_coefficient
            else:
                clustering_coefficients[index]=0
        
        average_coefficient=average_coefficient/len(array[0,:])

        return clustering_coefficients,average_coefficient



if __name__ == "__main__":
    
    #load in a a sparse matrix
    array_CCUGCU=sp.load_npz("testtxt/CCUGCU_Replicate_Average.npz")
    array_CCUCGU=sp.load_npz("testtxt/CCUCGU_Replicate_Average.npz")
    
    
    #Class Initiated
    test=Spectral_Analysis(array=array_CCUGCU)
    print(f"our array is being loaded and initialized:\n{test.array[0,:]}")
    test.create_degree_matrix()
    test.calculate_clustering_coefficients()
    