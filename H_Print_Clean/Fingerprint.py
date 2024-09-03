import matplotlib as mpl
from Convenience import n2_residue_numbers, test_topology, restrained_residue_list, test_residues
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np
import mdtraj as md

class Fingerprint_maker:
    def __init__(self,array,residues_of_interest=None,top=None):
    
        self.array_one = array
        self.namedict = None #Weir lab convenience case
        self.top = md.load_topology(top) if top is not None else None

        if residues_of_interest is n2_residue_numbers:
            self.residues=residues_of_interest
        else:
            self.residues = {residue.resSeq : f"{residue.name}:{residue.resSeq+1}" for residue in self.top.residues}
         
    def Base_Fingerprint(self,name,array=None):
        array = array if array is not None else self.array_one
        tick_positions,ticklabs=self.filtered_ticklab(array,residues=self.residues,topology=self.top)

        plot_array = array[1:, 1:]
        
        # Calculate absolute values
        abs_plot_array = np.abs(plot_array)
        
        scalemin, scalemax = self.scaleminmax(abs_plot_array)

        # Set the figure size (width, height in inches)
        plt.figure(figsize=(12, 8))

        print(f"{scalemin} {scalemax}")

        plt.imshow(abs_plot_array, cmap='Oranges', interpolation='nearest', aspect='auto', vmin=scalemin, vmax=scalemax)
        
        cbar = plt.colorbar()
        cbar.set_label('Scale') 

        plt.title('Average Hydrogen Bonding', loc='center')
        
        # Extract every 10th value from the first row for tick labels
        tick_values = array[0, ::10]  # Values from every 10th column of the first row
        tick_positions = np.arange(0, len(tick_values) * 10, 10)  # Corresponding positions for these values
        ticklabs = [f'{int(v):.0f}' for v in tick_values]  # Convert the values to string labels

        plt.xticks(np.arange(0, len(tick_positions), 1), labels=ticklabs, rotation=90, ha='right', fontsize=13)
        plt.yticks(np.arange(0, len(tick_positions), 1), labels=ticklabs)

        plt.gca().xaxis.set_ticks_position('top')
        print(f"before saving our name is {name}")

        plt.savefig(f'{name}Average_Heatmap.png', dpi=300, bbox_inches='tight')
        plt.clf()
        return

    def Filtered_Fingerprint(self,name,array=None,residues_to_filter=None):
        residues_to_filter= residues_to_filter if residues_to_filter is not None else self.residues
        array=array if array is not None else self.array_one

        tick_positions,ticklabs=self.filtered_ticklab(array,residues=self.residues,topology=self.top)
        scalemin,scalemax=self.scaleminmax(array)
    
        vcenter=(scalemin+scalemax)/2
        if scalemax <= 0:
            print(scalemax)
            # Define the custom colormap
            colors = [(0, 'orange'),  # Blue for negative values
                    (0.5, 'white'),  # White for zero
                    (1, 'blue')]  # Orange for positive values
            cmap = mcolors.LinearSegmentedColormap.from_list('orange_blue', colors)
            norm = mcolors.TwoSlopeNorm(vmin=scalemin, vcenter=vcenter, vmax=scalemax)
        else:
            # Define the custom colormap
            colors = [(0, 'blue'),  # Blue for negative values
                    (0.5, 'white'),  # White for zero
                    (1, 'orange')]  # Orange for positive values
            cmap = mcolors.LinearSegmentedColormap.from_list('blue_orange', colors)
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
        array = array if array is not None else self.array_one
        residues_to_filter = residues_to_filter if residues_to_filter is not None else self.residues
        
        tick_positions, ticklabs = self.filtered_ticklab(array, residues=self.residues, topology=self.top)
        scalemin,scalemax = self.scaleminmax(array)
        print(scalemin, scalemax)
        
        plot_array = array[1:, 1:]


        # Normalize the colormap to center at the average of scalemin and scalemax
        vcenter = (scalemin + scalemax) / 2
        print(vcenter)

        if scalemin<=0 and scalemax>=0:
            # Define the custom colormap
            colors = [(0, 'blue'),  # Blue for negative values
            (.5, 'white'), # White for zero
            (1, 'orange')] # Orange for positive values
            cmap = mcolors.LinearSegmentedColormap.from_list('blue_orange', colors)
            norm = mcolors.TwoSlopeNorm(vmin=scalemin, vcenter=vcenter, vmax=scalemax)
        
        else:
            print(scalemax)
                                # Define the custom colormap
            colors = [(0, 'white'),  # Blue for negative values
            (.5, 'white'), # White for zero
            (1, 'orange')] # Orange for positive values
            cmap = mcolors.LinearSegmentedColormap.from_list('orange_blue', colors)
            norm = mcolors.TwoSlopeNorm(vmin=scalemin,vcenter=vcenter, vmax=scalemax)

        plt.imshow(plot_array, cmap=cmap, interpolation='nearest', aspect='auto')
        
        cbar = plt.colorbar()
        cbar.set_label('Scale') 

        plt.title('Average Hydrogen Bonding')

        plt.xticks(np.arange(0, len(ticklabs), 1), labels=ticklabs, rotation=90, ha='right', fontsize=13)
        plt.yticks(np.arange(0, len(ticklabs), 1), labels=ticklabs)

        plt.gca().xaxis.set_ticks_position('top')
        print(f"before saving our name is {name}")

        plt.savefig(f'{name}_Average_Heatmap.png', dpi=300, bbox_inches='tight')
        plt.clf() 
        return

    def filtered_difference_Fingerprint(self, name, array=None, residues_to_filter=None):
        residues_to_filter= residues_to_filter if residues_to_filter is not None else self.residues
        array=array if array is not None else self.array_one

        tick_positions,ticklabs=self.filtered_ticklab(array,residues=self.residues,topology=self.top)
        scalemin,scalemax=self.scaleminmax(array)
    
        vcenter=(scalemin+scalemax)/2
        if scalemax <= 0:
            print(scalemax)
            # Define the custom colormap
            colors = [(0, 'orange'),  # Blue for negative values
                    (0.5, 'white'),  # White for zero
                    (1, 'blue')]  # Orange for positive values
            cmap = mcolors.LinearSegmentedColormap.from_list('orange_blue', colors)
            norm = mcolors.TwoSlopeNorm(vmin=scalemin, vcenter=vcenter, vmax=scalemax)
        else:
            # Define the custom colormap
            colors = [(0, 'blue'),  # Blue for negative values
                    (0.5, 'white'),  # White for zero
                    (1, 'orange')]  # Orange for positive values
            cmap = mcolors.LinearSegmentedColormap.from_list('blue_orange', colors)
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

    def filtered_ticklab(self,array,residues=None,topology=None):
        residues = residues if residues is not None else self.residues
        topology = topology if topology is not None else self.top
        if residues is n2_residue_numbers:
            tick_values = array[0,1:] # Values from every 10th column of the first row
            ticklabs =[f'{self.residues[i+1]}' if i != 0 else ' ' for i in tick_values]
        else:
            tick_values = array[0,1:] 
            ticklabs = [f'{self.residues[i]}' if i != 0 else ' ' for i in tick_values] 

        return tick_values,ticklabs
    
    def scaleminmax(self,array):
        array=array[1:,1:]
        scalemin=np.min(array)
        scalemax=np.max(array)
        return scalemin, scalemax
        
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
  
    #load previously created data
    array_AAUGCU=sp.load_npz("/home66/mraval/tRNAmod/Hbond_fingerprint/AAUGCU_nomod_Trajectory_Average.npz")
    array_AAUCGU=sp.load_npz("/home66/mraval/tRNAmod/Hbond_fingerprint/AAUCGU_nomod_Trajectory_Average.npz")

    #Initiate Data_Manipulator
    array_AAUGCU=Data_Manipulator(array_one=array_AAUGCU,residues_to_filter=n2_residue_numbers,res_filtered=restrained_residue_list,topology=md.load_topology("/home66/mraval/tRNAmod/34GUU36/AAUGCU_nomod/TLEAP/5JUP_N2_nowat.prmtop"))
    array_AAUCGU=Data_Manipulator(array_one=array_AAUCGU,residues_to_filter=n2_residue_numbers,res_filtered=restrained_residue_list,topology=md.load_topology("/home66/mraval/tRNAmod/34GUU36/AAUCGU_nomod/TLEAP/5JUP_N2_nowat.prmtop"))
    #Initiate Data_Manipulator

    #Create Difference Array between CCU and GCU all arrays
    full_array_one=array_AAUGCU.process_input()
    full_array_two=array_AAUCGU.process_input()
    difference_AAU= np.copy(full_array_one)
    difference_AAU[1:, 1:] =  full_array_two[1:, 1:]-full_array_one[1:, 1:]
    print(len(difference_AAU[0,:]))

    #Create Difference Array between CCU and GCU all arrays
    diff_AAUCGU=Data_Manipulator(array_one=difference_AAU,residues_to_filter=n2_residue_numbers,res_filtered=restrained_residue_list,topology=md.load_topology("/home66/mraval/tRNAmod/34GUU36/AAUCGU_nomod/TLEAP/5JUP_N2_nowat.prmtop"))

    final_filtered=diff_AAUCGU.filter_all_over_diff(array=difference_AAU)
    print(len(final_filtered[0,:]))
    

    Nozero_AAUCGU=Fingerprint_maker(array=final_filtered,residues_of_interest=None,top=md.load_topology("/home66/mraval/tRNAmod/34GUU36/AAUCGU_nomod/TLEAP/5JUP_N2_nowat.prmtop"))

    Nozero_AAUCGU.filtered_difference_Fingerprint("AAU_Fingerprints_filtfing")
    Nozero_AAUCGU.Base_Fingerprint("AAU_Fingerprints_baseging")
    
    os._exit(0)
    Full_CCUGCU.Base_Fingerprint("CCUGCU_Unrestrained")
    Full_CCU_Diff.Base_Fingerprint("CCU_Diff")

    Filtered_CCUGCU.Filtered_Fingerprint("Filtered_CCUGCU_Unrestrained")
    Filtered_CCU_Diff.Filtered_Fingerprint("Filtered_CCU_Diff")

    
    
    os._exit(0)

    Test_three.difference_Fingerprint("CCUGCU_P_Site_diff")
    full_array_one=Test_uno.process_input()
    full_array_two=Test_two.process_input()


    Filtered_array_one=Test_uno.filter_array()
    Filtered_array_two=Test_two.filter_array()
        

    Test_one.Base_Fingerprint("CCUCGU_TestLabs")
    os._exit(0)
    Test_two.Base_Fingerprint("CCUGCU_P_Site_")


    Test_one.Filtered_Fingerprint("CCUCGU_P_Site_")
    Test_two.Filtered_Fingerprint("CCUGCU_P_Site_")



    os._exit(0)

    #load in a a sparse matrix
    #array_uno=sp.load_npz('/zfshomes/lperez/fingerprint/H_Print/CCUCGU_G34_test_sparse_matrix.npz')

    #Class Initiated


    #Initiate fingerprint maker class


    #from visualization import display_arrays_as_video

    # a nice weir lab sanity check is built in simply set
    # residues of interest to n2_residue_numbers in order
    # to assure the data is valid
    
    
    H_Print_CCUCGU_G34=Fingerprint_maker(
        array_one=np.loadtxt('/zfshomes/lperez/fingerprint/H_Print/H_Print_CCUCGU_G34/CCUCGU_G34_average.txt'),
        top=md.load_topology('/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/TLEAP/5JUP_N2_CGU_nowat.prmtop'),
        residues_of_interest=n2_residue_numbers
        )
    
    H_Print_CCUGCU_G34=Fingerprint_maker(
            array_one=np.loadtxt('/zfshomes/lperez/fingerprint/H_Print/H_Print_CCUGCU_G34/CCUGCU_G34_average.txt'),
            top=md.load_topology('/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop'),
            residues_of_interest=n2_residue_numbers
            )
    
    #full sized fingerprints with correct indexes added shortly
    H_Print_CCUCGU_G34.Raw_Fingerprint("H_Print_CCUCGU_G34_raw")
    H_Print_CCUGCU_G34.Raw_Fingerprint("H_Print_CCUCGU_G34_raw")
    #full sized fingerprints with correct indexes added shortly

    

    filtered_CCUCGU_G34=H_Print_CCUCGU_G34.filter_array(residues_to_filter=n2_residue_numbers)
    filtered_CCUGCU_G34=H_Print_CCUGCU_G34.filter_array(residues_to_filter=n2_residue_numbers)


    test_difference_array = H_Print_CCUCGU_G34.difference_array(first_matrix=filtered_CCUCGU_G34,second_matrix=filtered_CCUGCU_G34)


    H_Print_CCUCGU_G34.difference_Fingerprint(name="test_difference",array=test_difference_array)


    #example using weir lab sanity check
    H_Print_CCUGCU_G34.residues_of_interest=n2_residue_numbers
    filtered_array=H_Print_CCUGCU_G34.filter_array()
    np.savetxt('N2_Residue_CCUGCU_G34.txt', filtered_array, fmt='%.8f')

    H_Print_CCUCGU_G34.residues_of_interest=n2_residue_numbers
    filtered_array_CCUCGU_=H_Print_CCUCGU_G34.filter_array()
    np.savetxt('N2_Residue_CCUCGU_G34.txt', filtered_array, fmt='%.8f')

    #example using weir lab sanity check
    H_Print_CCUGCU_G34.Filtered_Fingerprint("Results/CCUGCU_A-Site_Residues")
    H_Print_CCUCGU_G34.Filtered_Fingerprint("Results/CCUCGU_A-Site_Residues")
    #example using regular numbering on Residues of CAR and the +1 Codon
    H_Print_CCUGCU_G34.Filtered_Fingerprint("Results/CCUGCU_Car-+1_Residues",list(CAR_plus1_Codon.keys()))  
    H_Print_CCUCGU_G34.Filtered_Fingerprint("Results/CCUCGU_Car-+1_Residues",list(CAR_plus1_Codon.keys()))
    
  
    #example using regular numbering on Residues of CAR and the +1 Codon
    topology=md.load_topology("/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop")
    file=np.loadtxt("/zfshomes/lperez/fingerprint/H_Print/CCUGCU_G34_KB_Average.txt")
    KB_Print_CCUCGU_G34=Fingerprint_maker(array_one=file,top=topology,residues_of_interest=n2_residue_numbers)
    KB_Print_CCUCGU_G34.Raw_Fingerprint("KBTEST_H_Print_CCUCGU_G34")


    if test == True:

        #producing some stats to compare with what we know is true GCU
        H_Print_CCUGCU_G34.filter_op = 0
        H_Print_CCUGCU_G34.average_aggregator(H_Print_CCUGCU_G34.array_one,list(CAR_plus1_Codon.keys()))

        #producing some stats to compare with what we know is true
        H_Print_CCUCGU_G34.filter_op = 0
        H_Print_CCUCGU_G34.average_aggregator(H_Print_CCUCGU_G34.array_one,list(CAR_plus1_Codon.keys()))

        
        #Call the function with just one heatmap
       
        #Use with H_Print
        H_Print_CCUGCU_G34.Filtered_Fingerprint("")
        single_heatmap=H_Print_CCUGCU_G34.filter_array(residues_to_filter=CAR_plus1_Codon)
        single_heatmap=single_heatmap[1:,1:]
        #display_arrays_as_video([single_heatmap],res_indicies=list(range(10)),seconds_per_frame=1,outfile_prefix='testing',)



        



        

