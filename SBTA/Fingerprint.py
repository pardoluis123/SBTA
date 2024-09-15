import matplotlib as mpl
from Convenience import n2_residue_numbers, test_topology, restrained_residue_list, test_residues
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np
import mdtraj as md

class Fingerprint_maker:
    def __init__(self,array,residues_of_interest=None,top=None):
    
        self.array = array
        self.namedict = None #Weir lab convenience case
        self.top = md.load_topology(top) if top is not None else None

        if residues_of_interest is n2_residue_numbers:
            self.residues=residues_of_interest
        else:
            self.residues = {residue.resSeq : f"{residue.name}:{residue.resSeq+1}" for residue in self.top.residues}
    
    def Base_Fingerprint(self,name,array=None):
        residues_to_filter= residues_to_filter if residues_to_filter is not None else self.residues
        array=array if array is not None else self.array

        tick_positions,ticklabs=self.filtered_ticklab(array,residues=self.residues,topology=self.top)
        scalemin,scalemax=self.scaleminmax(array)
    
        vcenter=(scalemin+scalemax)/2

        if scalemax >= 0:
            # Define the custom colormap
            colors = [(0, 'white'),  # Blue for negative values
                    (0.5, 'white'),  # White for zero
                    (1, 'orange')]  # Orange for positive values
            cmap = mcolors.LinearSegmentedColormap.from_list('orange_blue', colors)
            norm = mcolors.TwoSlopeNorm(vmin=scalemin, vcenter=vcenter, vmax=scalemax)
        if scalemax <= 0:
            # Define the custom colormap
            colors = [(0, 'orange'),  # Blue for negative values
                    (0.5, 'white'),  # White for zero
                    (1, 'white')]  # Orange for positive values
            cmap = mcolors.LinearSegmentedColormap.from_list('orange_blue', colors)
            norm = mcolors.TwoSlopeNorm(vmin=scalemin, vcenter=vcenter, vmax=scalemax)
        

        plot_array=array[1:,1:]

        plt.imshow(plot_array, cmap=cmap, interpolation='nearest', aspect='auto', norm=norm)        
        cbar = plt.colorbar()
        cbar.set_label('Scale')

        plt.title('Average Hydrogen Bonding')

        plt.xticks(np.arange(0, len(ticklabs), 1), labels=ticklabs, rotation=90, ha='right', fontsize=13)
        plt.yticks(np.arange(0, len(ticklabs), 1), labels=ticklabs)

        plt.gca().xaxis.set_ticks_position('top')
        print(f"before saving our name is {name}")

        plt.savefig(f'{name}_Average_Heatmap.png', dpi=300, bbox_inches='tight')
        plt.clf()

    def Filtered_Fingerprint(self,name,array=None,residues_to_filter=None):
        residues_to_filter= residues_to_filter if residues_to_filter is not None else self.residues
        array=array if array is not None else self.array

        tick_positions,ticklabs=self.filtered_ticklab(array,residues=self.residues,topology=self.top)
        scalemin,scalemax=self.scaleminmax(array)
    
        vcenter=(scalemin+scalemax)/2

        if scalemax <= 0:
            print(scalemax)
            # Define the custom colormap
            colors = [(0, 'white'),  # Blue for negative values
                    (0.5, 'white'),  # White for zero
                    (1, 'orange')]  # Orange for positive values
            cmap = mcolors.LinearSegmentedColormap.from_list('orange_blue', colors)
            norm = mcolors.TwoSlopeNorm(vmin=scalemin, vcenter=vcenter, vmax=scalemax)
        

        plot_array=array[1:,1:]

        plt.imshow(plot_array, cmap=cmap, interpolation='nearest', aspect='auto', norm=norm)        
        cbar = plt.colorbar()
        cbar.set_label('Scale')

        plt.title('Average Hydrogen Bonding')

        plt.xticks(np.arange(0, len(ticklabs), 1), labels=ticklabs, rotation=90, ha='right', fontsize=13)
        plt.yticks(np.arange(0, len(ticklabs), 1), labels=ticklabs)

        plt.gca().xaxis.set_ticks_position('top')
        print(f"before saving our name is {name}")

        plt.savefig(f'{name}_Average_Heatmap.png', dpi=300, bbox_inches='tight')
        plt.clf()
    
    def difference_Fingerprint(self, name, array=None, residues_to_filter=None):
        array = array if array is not None else self.array
        residues_to_filter = residues_to_filter if residues_to_filter is not None else self.residues
        
        #pulling together some tick positions and labels; and the colorbar min and max
        tick_positions, ticklabs = self.filtered_ticklab(array, residues=self.residues, topology=self.top)
        scalemin,scalemax = self.scaleminmax(array) 
        plot_array = array[1:, 1:]

        # Normalize the colormap to center at the average of scalemin and scalemax
        vcenter = (scalemin + scalemax) / 2
        print(f"scalemin:{scalemin}\nscalemax:{scalemax}\ncenter={0}")

        # Define the custom colormap
        colors = [(0, 'blue'),  # Blue for negative values
        (.5, 'white'), # White for zero
        (1, 'orange')] # Orange for positive values
        cmap = mcolors.LinearSegmentedColormap.from_list('blue_orange', colors)
        
        #wether we center the color at 0 or at the middle depends on if we have negative values
        if scalemin<=0.0 and scalemax>0.0:
            norm = mcolors.TwoSlopeNorm(vmin=scalemin, vcenter=0, vmax=scalemax)
        else:
            norm = mcolors.TwoSlopeNorm(vmin=scalemin,vcenter=vcenter, vmax=scalemax)

        #plot
        plt.imshow(plot_array, cmap=cmap, interpolation='nearest', aspect='auto',norm=norm)

        #Righthand Side Colorscale labelling
        cbar = plt.colorbar()
        cbar.set_label('Scale') 
        custom_ticks = self.difference_colorbar_ticks() #custom ticks
        cbar.set_ticks(custom_ticks)

        #title and ticks
        plt.title('Average Hydrogen Bonding')
        plt.xticks(np.arange(0, len(ticklabs), 1), labels=ticklabs, rotation=90, ha='right', fontsize=13)
        plt.yticks(np.arange(0, len(ticklabs), 1), labels=ticklabs)
        plt.gca().xaxis.set_ticks_position('top')

        #Output
        print(f"before saving our name is {name}")
        plt.savefig(f'{name}_Average_Heatmap.png', dpi=300, bbox_inches='tight')
        plt.clf() 
        return

    def filtered_ticklab(self,array,residues=None,topology=None):
        residues = residues if residues is not None else self.residues
        topology = topology if topology is not None else self.top
        if residues is n2_residue_numbers:
            tick_values = array[0,1:] # Values from every 10th column of the first row
            ticklabs =[f'{self.residues[i+1]}' if i != 0 else ' ' for i in tick_values]
        else:
            tick_values = array[0,1:] 
            ticklabs = [f'{self.residues[i-1]}' if i != 0 else ' ' for i in tick_values] 

        return tick_values,ticklabs
    
    def difference_colorbar_ticks(self,array=None):
        array=array if array is not None else self.array
        amin,amax=self.scaleminmax(array)
        scale=[]
        if abs(amin)>abs(amax):
            
            while amin<amax:
                scale.append(amin)
                amin=amin+amax
            scale.append(amax)
            return scale
        else:
            while amax>amin:
                scale.insert(0,amax)
                amax=amax-amin
            scale.append(amax)
            return scale

    def scaleminmax(self,array):
        array=array[1:,1:]
        scalemin=np.min(array)
        scalemax=np.max(array)
        return scalemin, scalemax
    
    def multi_print_maker(self,array=None,out_directory=None,name=None):
        """
        Creates multiple fingerprints at once from an array of arrays

        Parameters
        ----------
        Array:numpy.ndarray
            Array containing other adjacency matrices which will be used to make comparisons
        Array Two:out_directory
            The directory where the images should all be stored     
        """
        
        array=array if array is not None else self.array
        name=name if name is not None else "filler_name"

        n=1
        for i in array:
            array=np.loadtxt(i)
            cpptraj_fmaker=Fingerprint_maker(array=cpptraj_array,top=self.top,residues_of_interest=n2_residue_numbers)
            cpptraj_fmaker.Base_Fingerprint(name="n2_cpptrajdiff")

    def average_aggregator(self,array,residues_of_interest):
        array=self.filter_array(array,residues_of_interest)
        array=array[1:,1:] #1dimension
        array=array.flatten()
        array_values=np.unique(array) #no duplicates so 1 side of array
        total_average_hydrogen=np.divide((np.sum(array_values)),len(residues_of_interest))
        print(f"CGU filtered array\n {array}\n \nvalues{array_values}\n total average {total_average_hydrogen}")

if __name__ == "__main__":
    from Convenience import CAR_plus1_Codon, n2_residue_numbers,p_site,restrained_residue_list,test_residues
    import os
    import scipy.sparse as sp 
    from Data_Manip import Data_Manipulator
  
    cpptraj_array=np.loadtxt("/zfshomes/lperez/fingerprint/H_Print/mdtraj_difference.txt")    
    cpptraj_fmaker=Fingerprint_maker(array=cpptraj_array,top="/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop",residues_of_interest=n2_residue_numbers)
    cpptraj_fmaker.Base_Fingerprint(name="n2_cpptrajdiff")
    

