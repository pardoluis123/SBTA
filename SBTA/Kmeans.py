import numpy as np
import scipy.sparse as sp
import scipy.cluster.vq as vq
import scipy.cluster as cl
from scipy.stats import zscore as zs
import mdtraj as md
from Data_Manip import Data_Manipulator
from Convenience import n2_residue_numbers, test_topology,restrained_residue_list,test_residues
import os

class Kmeans_Cluster():

    def __init__(self,array=None):
        self.process_input=Data_Manipulator.process_input
        self.array=self.process_input(array) if array is not None else None

    def format_traj_for_clust(self,array=None)->np.ndarray:
        """
        returns a formatted array of choice for kmeans clustering
        *preferrably this is a trajectory array of frame arrays*

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest (should be a 3dimensional array holding 2dimensional arrays i.e. xyz->xy)     
        """
        
        array=array if array is not None else self.array

        for idx in range(len(array)): 
            array[idx] = array[idx].flatten() 

        return array
    
    def format_replicate_for_clust(self,array=None)->np.ndarray:
        """
        Returns a formatted array of replicates for kmeans clustering
        *Input should be a 4d array holding an array of replicates, each holding 2d frames)

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest (should be a 3dimensional array holding 2dimensional arrays i.e. xyz->xy)     
        """

        array=array if array is not None else self.array

        n=0
        for i in array:#iterate through replicates
            self.format_traj_for_clust(array=i) #reformat the frames and return
        
        return array

    def format_array_for_whitening(self,array=None)->np.ndarray:
        """
        Returns a frame array which has removed all columns containing only zeroes so whitening
        is not affected as theese variances are not meaningful for clustering and may affect data
        *Input should be a 4d array holding an array of replicates, each holding 2d frames*
        *whitened=normalized*

        Parameters
        ----------
        Array:numpy.ndarray
            Frame array of interest (2d adjacency matriz)  
        """
        array=array if array is not None else None
        

        # Remove zero-variance columns 
        non_zero_mask = array != 0

        array = array[non_zero_mask]  # Filter out zero-variance columns
        return array

    def whiten_frames(self,array=None)->np.ndarray:
        """
        Returns trajectory array in which all frames have been iterated through and "whitened"
        as in all zero values have been removed and the frames have been normalized

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest a one dimensional array holding flattened adjacency matrices 
        """
        
        array=array if array is not None else self.array
        #for frame in replicate
        for idx in range(len(array)):
            current_frame = self.format_array_for_whitening(array=array[idx])
            whitened_frame=vq.whiten(current_frame)
            array[idx]=whitened_frame

        return array
    
    def zscore_frames(self,array=None)->np.ndarray:
        """
        Returns trajectory array in which all frames have been iterated through and zscore normalized
        as in all zero values have been removed and the frames have been normalized

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest a one dimensional array holding flattened adjacency matrices 
        """
        array=array if array is not None else self.array
        
        #for frame in replicate
        for idx in range(len(array)):
            zframe=zs(array[idx])
            array[idx]=zframe

        return array

    def whiten_replicates_cluster(self,array=None)->np.ndarray:
        """
        Returns array replicate array after iterating through each replicated and "whitening" the frames
        so the values hold similar weight for kmeans clustering

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest (one dimensial array of frames holding flattened adjacency matrix frames)
        """

        array=array if array is not None else self.format_replicate_for_clust()
        
        #for replicate
        for idx in range(len(array)): 
            current_replicate=array[idx]
            array[idx]=self.whiten_frames(array=current_replicate)
            
        return array
    
    def zscore_replicates_cluster(self,array=None)->np.ndarray:
        """
        Returns array replicate array after iterating through each replicated and zscore normalizing frames
        so the values hold similar weight for kmeans clustering

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest (one dimensial array of frames holding flattened adjacency matrix frames)
        """
        array=array if array is not None else self.format_replicate_for_clust()
        
        #for replicate
        for idx in range(len(array)): 
            current_replicate=array[idx]
            array[idx]=self.zscore_frames(array=current_replicate)

        return array
    
    def create_empty_clustering_table(self,array=None,n=None)->np.ndarray:
        """
        Returns an empty array to hold data derived from kmeans procedure

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest (one dimensial array of frames holding flattened adjacency matrix frames)
        """

        array= array if array is not None else self.zscore_replicates_cluster()
        n = n if n is not None else 2

        #create final new array
        labels=[num for num in range(0,n)]
        final_array=np.empty((n,2), dtype=object)
        final_array[:,0]=labels
        return final_array     

    def iterate_through_labels(self,array=None,n=None)->np.ndarray:
        
        pass

    def preform_clust_w(self,array=None,n=None)->np.ndarray:
        """
        Preforms Kmeans Clustering on an array of choice(whitening normalization)
        *preferrably this is an array of arrays either a trajectory array of frame arrays or replicate array of replicates)

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest
        threshhold:float
            float value at which the array will be filtered for "significant" values     
        """

        array=array if array is not None else self.whiten_replicates_cluster()
        cluster_output=[]
        for idx in range(1,len(array)):
            array[idx] = np.vstack(array[idx])  
            k_array=vq.kmeans(obs=array[idx],k_or_guess=8)
            cluster_output.append(k_array)
        
        return k_array

    def preform_clust_z(self,array=None,n=None)->np.ndarray:
        """
        Returns array of kmeans data after preforming Kmeans Clustering on an array of choice(zscore normalization)
        *preferrably this is an array of arrays either a trajectory array of frame arrays or replicate array of replicates)

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest
        threshhold:float
            float value at which the array will be filtered for "significant" values     
        """

        array=array if array is not None else self.zscore_replicates_cluster()
        n=n if n is not None else 2

        final_data_holder=[]
        
        #iterate through and preform Kmeans on each array
        for idx in range(1,len(array)):
            #for each array preform kmeans using scipy's kmeans function
            empty_table=self.create_empty_clustering_table(n=n)
            array[idx] = np.vstack(array[idx])  
            codebook, distortion=vq.kmeans(obs=array[idx],k_or_guess=2) 

            #iterate through and update each new empty table
            j=0
            for i in codebook:  
                empty_table[j,1]=i
                j+=1
            
            final_data_holder.append(empty_table)

        final_data_holder=np.array(final_data_holder,dtype=object)

        return final_data_holder

    def save_cluster_array(self,ntype=None,array=None)->None:
        """
        Returns None but saves Kmeans Data Array as numpy datafile

        Parameters
        ----------
        Array:numpy.ndarray
            Array of interest
        ntype:the type of kmeans desired
            float value at which the array will be filtered for "significant" values     
        """

        ntype = ntype if ntype is not None else 0

        #based on normalization type generate data
        if ntype == 0:
            array=array if array is not None else self.preform_clust_z()
        elif ntype ==1:
            array=array if array is not None else self.preform_clust_w()
        
        np.save(f"{out_GCU}_reparray_first_pass_cluster",array)


        '''#save generated data
        n=1
        for i in array:
            np.save(f"{out_GCU}_rep_{n}_first_pass_cluster",i[0])
            print(i[1])
            n+=1'''
        
        return 

if __name__=="__main__":
    array_CCUGCU="/zfshomes/lperez/fingerprint/H_Print/H_Print_CCUCGU_G34/CCUCGU_G34_Replicate_Array.npy"
    out_GCU="/zfshomes/lperez/fingerprint/H_Print/cluster_output_CCUGCU/"
    
    #array_CCUCGU="/zfshomes/lperez/fingerprint/H_Print/testtxt/test_frame_Replicate_Average.npz"    
    
    #Class Initiated
    Test=Kmeans_Cluster(array=array_CCUGCU)
    Test.save_cluster_array()
    
    #path to previously saved array
    realpath="/zfshomes/lperez/fingerprint/H_Print/cluster_output_CCUGCU/_reparray_first_pass_cluster.npy"
    