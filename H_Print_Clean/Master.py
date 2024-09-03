import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib as mpl
import time
import os
import tempfile
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from visualization import create_axis_labels

class Trajectory:
    #A class for importing and manipulating trajectories as neeeded

    #takes 4 arguments
        # - the trajectory file
        # - the topology file
        # - A list of residues of interest
        # - The option to either remove "0" or highlight only "1" said residues, if none are specified the analysis is done with all residues
            #-note the only mandatory parameters are trajectory and topology

    #Is initiated with four attributes
        # - The Trajectory
        # - An array with X and Y axis set to the residues of interest
        # = A dictionary containing all atom-residue parings,
        # = An empty list for processing

    #initiating instance of object
    def __init__(self,trajectory,topology,residues=None, slice_option=None):
        #Attributes
        self.residues = residues if residues is not None else []
        self.slice_option = slice_option if slice_option is not None else "default"
        self.trajectory=self.filter(md.load(trajectory,top=topology),self.residues,self.slice_option)
        self.array_of_residues,self.atom_to_residue,self.Trajectory_array=self.initiate_array()

    #Function For Filtering
    def filter(self,trajectory,residues,slice_option): 
    #fixing indexing because they are all one indexed in mdtraj

        #0 indexing input files from Amber
        original_residues=[residue.index for residue in trajectory.topology.residues]
        residues_query = {resnum-1 for resnum in residues} 

        if (slice_option==1):
            residues_query = " or ".join([f"residue == {resnum}" for resnum in residues]) 
            atom_indices_selection = trajectory.top.select(residues_query) 
            filtered_trajectory = trajectory.atom_slice(atom_indices_selection) 

        if (slice_option==0):
            mask = [atom.index for atom in trajectory.top.atoms if atom.residue.index not in residues]
            residues_query = " or ".join([f"residue == {resnum}" for resnum in residues]) 
            atom_indices_selection = trajectory.top.select(residues_query) 
            filtered_trajectory = trajectory.atom_slice(atom_indices_selection) 

        else:
            filtered_trajectory=trajectory
        
        self.residues=[residue.index for residue in filtered_trajectory.topology.residues]
        #print(f"our original residues: {original_residues}\n Our new set: {self.residues}")
        return filtered_trajectory

    #Initiate Array
    def initiate_array(self,):
        Trajectory_array=[]
        atom_to_residue = {atom.index: atom.residue.resSeq for atom in self.trajectory.topology.atoms}
        indexes = [residue.resSeq for residue in self.trajectory.topology.residues]
        residue_to_index = [len(indexes)]

        
        array_of_residues = np.zeros((len(indexes), len(indexes)))
        array_of_residues[0,:]=indexes
        array_of_residues[:,0]=indexes

        return array_of_residues,atom_to_residue,Trajectory_array

class H_Print:
    #A class for creating Hydrogen Bond Prints from the data 

    #Takes 2 Parameters 
        # Trajectory object must be initiated with "Trajectory"
        # The number of frames the trajectory contains
    #Is initiated with four attributes
        # - The Trajectory
        # - An array with X and Y axis set to the residues of interest
        # = A dictionary containing all atom-residue parings,
        # = An empty list for processing

    def __init__(self, Trajectory):
        self.Trajectory=Trajectory
        self.Trajectory_array,self.final_average_Hbond=self.frame_iterator()
        
    #Update hbond counts
    def update_hbond_counts(self,hydrogen_matrix,frame_array):

        indexes=frame_array[0,:]
        print(indexes)

        for d_i, h_i, a_i in hydrogen_matrix:

            donor_residue_seq = self.Trajectory.atom_to_residue[d_i]
            acceptor_residue_seq = self.Trajectory.atom_to_residue[a_i]
            
            indices_1 = np.where(indexes == donor_residue_seq)[0]
            indices_2 = np.where(indexes == acceptor_residue_seq)[0]
            print(f"our donor {donor_residue_seq} and {acceptor_residue_seq}\n followed by indices {indices_1} and {indices_2} ")

            frame_array[indices_1,indices_2]+=1
            frame_array[indices_2,indices_1]+=1

                 
        frame_array=frame_array[1:,1:]
        return frame_array

    #iterate through frames of trajectory
    def Trajectory_Processor(self,Trajectory_array,array_of_residues):

        #Initialize Arrays
        average_Hbond= np.copy(Trajectory_array)
        
        average_Hbond= np.mean(average_Hbond,axis=0)
        average_Hbond= np.round(average_Hbond, decimals=2)

        final_average_Hbond=np.copy(array_of_residues)
        final_average_Hbond[1:,1:]=average_Hbond[1:,1:]

        return final_average_Hbond

    #iterate through frames of trajectory
    def frame_iterator(self):

        n=0
        for i in self.Trajectory.trajectory:
            print(f"Processing Frame {n+1} of the current trajectory\n")
            frame_array=np.copy(self.Trajectory.array_of_residues)
            current_trajectory_frame=i
            hydrogen_matrix=md.baker_hubbard(current_trajectory_frame)
            print(f"our current matrix \n{hydrogen_matrix}")
            frame_array[1:,1:]=self.update_hbond_counts(hydrogen_matrix,frame_array) 
            self.Trajectory.Trajectory_array.append(frame_array)

            n+=1

        self.Trajectory.Trajectory_array = np.array(self.Trajectory.Trajectory_array)

        final_average_Hbond=self.Trajectory_Processor(self.Trajectory.Trajectory_array,self.Trajectory.array_of_residues)

        return self.Trajectory.Trajectory_array,final_average_Hbond

    #Export_Hbonds
    def Export_Hbonds(self,name):
        i=1
        for array in self.Trajectory_array:
            array = pd.DataFrame(array)
            array.to_csv(name+str(i)+'.csv', index=False, header=False)
            i+=1

        final_average_Hbond = pd.DataFrame(self.final_average_Hbond)
        final_average_Hbond.to_csv(name+'final_average_Hbond.csv', index=False, header=False)

    def Fingerprint(self,name):
        #test_plot=self.final_average_Hbond[self.final_average_Hbond>0]
        plt.imshow(self.final_average_Hbond, cmap='Oranges', interpolation='nearest', aspect='auto',vmin=0, vmax=4)
        
        # Add colorbar
        cbar = plt.colorbar()
        cbar.set_label('Scale')

        axis=len(self.final_average_Hbond[0,:])


        # Add title
        plt.title('Average Hydrogen Bonding')
        #plt.xticks(ticks=range(axis), labels=[f'Label {i}' for i in range(axis)])

        tick_locations, tick_labels=create_axis_labels(self.Trajectory.residues)

        plt.xticks(tick_locations, tick_labels)
        plt.yticks(tick_locations, tick_labels)

        plt.gca().xaxis.set_ticks_position('top')
        print(f"before saving our name is {name}")
        print(f"also testing before saving our name is {name}/")
        plt.savefig(f'{name}Average_Heatmap.png', dpi=300, bbox_inches='tight')
        plt.clf() 
        return

    def Replicate_average(self):
        pass

    def Histogram(self,Trajectory_array,Final_percent_frames,final_average_Hbond):
        final_average_Hbond=final_average_Hbond[1:,1:]
        final_average_Hbond = final_average_Hbond.flatten()
        bin_edges = np.histogram_bin_edges(final_average_Hbond, bins=5)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        plt.xticks(bin_centers, x_labels, rotation=45, ha='right')
        plt.title('Average Hydrogen Bonding')
        plt.savefig('histogram.png', dpi=300, bbox_inches='tight')
        plt.clf() 
        return

class FileManagement:
    #A class for Creating a "quick" Object with trajectories and lists of files with features like "strip_only"

    #Takes 2 Parameters 
        # The directory Trajectories are found in 
        # The parmtop file for the version of the trajectory being used
    #Is initiated with four attributes
        # - The path to the directory
        # - The path to the topology
        # = A list of trajectories
        # = Currently a list of only stripped trajectories (very modular could contain more)


    def __init__(self,directory_input,top_input):
        self.directory_input=directory_input
        self.topology=top_input
        self.filelist=self.Trajectory_management()
        self.strip_only=self.Strip_only()

    def Trajectory_management(self):
        filelist=os.listdir(self.directory_input)
        filelist= [self.directory_input+file for file in filelist]
        return filelist
    
    def Strip_only(self):
        strips = [file for file in self.filelist if "1in50" in file]
        return strips
        
class Run_Program:
    #A class for Creating a "quick" Object with trajectories and lists of files with features like "strip_only"

    #Takes 2 Parameters 
        # The directory Trajectories are found in 
        # The parmtop file for the version of the trajectory being used
    #Is initiated with four attributes
        # - The path to the directory
        # - The path to the topology
        # = A list of trajectories
        # = Currently a list of only stripped trajectories (very modular could contain more)
    
    def __init__(self,directory_input,top_input, residues, slice_option, name):
        self.directory_input=directory_input
        self.Files=FileManagement(directory_input,top_input)
        self.residues=residues
        self.slice_option=slice_option
        self.name=name

    def symbolic_link(self,i):

        symlink_dir=f"{os.getcwd()}/symlink/{self.name}/"
        os.makedirs(symlink_dir, exist_ok=True)
        symlink_path = os.path.join(symlink_dir, f"{os.path.basename(i)}.mdcrd")

        try:
            os.symlink(i, symlink_path)
            print(f"Created symbolic link: {i} -> {symlink_path}")
            return symlink_path

        except FileExistsError:
            return symlink_path

    def process_trajectory(self,i,exportdir):

        print(f"now accessing {i} ...\n")
        try:
            if ".mdcrd" not in i:
                symlink_path=self.symbolic_link(i)
                current_trajectory = Trajectory(symlink_path, self.Files.topology, self.residues, self.slice_option)
                print(f"Loaded trajectory for {i} ...\n")
            else:
                current_trajectory = Trajectory(i, self.Files.topology, self.residues, self.slice_option)
                print(f"Loaded trajectory for {i} ...\n")

            current_fingerprint=H_Print(current_trajectory)
            print(f"Fingerprint Object for {i} created ...\n")

            current_fingerprint.Export_Hbonds(f"{exportdir}/")
            print(f"files for {i} exported  ...\n")

            current_fingerprint.Fingerprint(f"{exportdir}/")
            print(f"Average hbonding heatmap for {i} exported  ...\n")

        finally:
            pass

    def create_directories(self):
            n=1
            current_directory=os.getcwd() 
            core_directory=f"{current_directory}/H_Print_{self.name}"
            try:
                os.listdir(core_directory)
            except FileNotFoundError:
                print(f"{core_directory} not found so creating a new directory\n")
                os.mkdir(core_directory)

            for i in self.Files.strip_only: 
                new_directory=f"{current_directory}/H_Print_{self.name}/Replicate_{n}"
                try:
                    dir=os.listdir(new_directory)
                    print(f"the following files were found inside {new_directory} \n {dir}")
                except FileNotFoundError:
                    print(f"{new_directory} not found so creating a new directory\n")
                    os.mkdir(new_directory)
                n+=1

    def iterate_replicates(self): 
        n=1
        self.create_directories()
        print(self.Files.strip_only)
        for i in self.Files.strip_only:
            current_directory=os.getcwd()
            destination=f"{current_directory}/H_Print_{self.name}/Replicate_{n}/"
            print(f"Data files for replicate {n} will be created at {destination} \n")
            self.process_trajectory(i,destination)
            n+=1

class Post_Processor:
    def __init__(self,directory,name):
        self.directory=directory
        self.name=name

    def iterate(self):

        working_directory=os.listdir(self.directory)
        master="/zfshomes/lperez/fingerprint/hydrogen/H_Print_CCUCGU_G34/Replicate_1/final_average_Hbond.csv"
        master = np.genfromtxt(master, delimiter=',', skip_header=1)
        print(master)
        array=[]
        for i in working_directory:
            if ".csv" not in i:
                current_replicate=f"{self.directory}/{i}/final_average_Hbond.csv"
                print(current_replicate)
                data = np.genfromtxt(current_replicate, delimiter=',', skip_header=1)
                data=data[1:,1:]
                array.append(data)
        array=np.array(array)
        master[1:,1:]=self.aggregate(array)

        final_avg = pd.DataFrame(master)
        final_avg.to_csv(f"{self.directory}/{self.name}.csv", index=False, header=False)
        
    def aggregate(self,array):

        average_Hbond= np.mean(array,axis=0)
        average_Hbond= np.round(average_Hbond, decimals=2)
        return average_Hbond

    def Fingerprint(self,name,array):
        #test_plot=self.final_average_Hbond[self.final_average_Hbond>0]
        plt.imshow(array, cmap='Oranges', interpolation='nearest', aspect='auto',vmin=0, vmax=4)
        
        # Add colorbar
        cbar = plt.colorbar()
        cbar.set_label('Scale')


        # Add title
        plt.title('Average Hydrogen Bonding')
        #plt.xticks(ticks=range(axis), labels=[f'Label {i}' for i in range(axis)])

        tick_locations, tick_labels=create_axis_labels(array)

        plt.xticks(tick_locations, tick_labels)
        plt.yticks(tick_locations, tick_labels)

        plt.gca().xaxis.set_ticks_position('top')
        print(f"before saving our name is {name}")
        print(f"also testing before saving our name is {name}/")
        plt.savefig(f'{name}Average_Heatmap.png', dpi=300, bbox_inches='tight')
        plt.clf() 
        return
        
        

if __name__ == "__main__":


    restrained_residue_list = [
    1, 6, 7,
    8, 9, 10, 11, 12, 17, 18, 19, 20, 21, 22, 32, 33,
    34, 35, 36,
    37, 38, 39, 40, 41, 42, 43, 44, 45, 52, 53,
    54, 55, 56, 57, 58, 59, 60, 61, 62, 64, 74, 77, 78,
    79, 80, 81, 82,
    83, 84, 85, 86, 87, 105, 106, 107, 108, 109, 110, 111,
    112, 113, 114, 115,
    116, 117, 118,
    119, 120, 121, 142, 143, 144, 145, 146, 147, 148, 149,
    150, 151, 152, 153, 154,
    155, 156, 157, 158, 159, 160, 161, 162, 163, 176, 177, 178, 179, 180, 181,
    182, 183, 196, 199, 200, 201, 202,
    203, 204, 205, 206,
    207, 208,
    209, 210, 211, 212, 213, 214, 215, 224, 225, 226, 227,
    228, 229, 247, 248, 249, 250, 251, 252, 253,
    254, 255, 270, 271,
    272, 273, 274, 275, 276, 277,
    278, 279, 280, 281, 282, 283, 284,
    285, 286, 287, 290, 291, 292, 293, 294, 295, 296,
    297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312,
    313, 314, 315, 316,
    317, 318, 319, 331,
    332, 333, 334, 335, 336,
    337, 338, 339, 340, 341,
    342, 343, 344, 345, 363, 364, 365,
    366, 367, 368, 369, 390, 391, 392, 393,
    394, 395, 396, 397,
    398, 399, 400, 416, 417, 418,
    419, 431,
    432, 433, 434, 435,
    436, 437, 438, 439, 440, 441, 442,
    443, 453, 454, 455, 456, 457, 458,
    459, 460, 469, 470, 471,
    472, 473, 474, 475,
    476, 477, 478, 479, 490, 491,
    492, 493, 494]

    CCUGCU_G34=Run_Program('/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/BKUP_5JUP_N2_GCU/',
                             '/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop',
                             restrained_residue_list,
                             0,
                             "CCUGCU_G34")

    CCUCGU="/zfshomes/lperez/fingerprint/hydrogen/H_Print_CCUCGU_G34/"
    CCUCGU_structure=Post_Processor(CCUCGU,"CCUCGU")
    
    CCUCGU_structure.Fingerprint()

    CCUGCU="/zfshomes/lperez/fingerprint/hydrogen/H_Print_CCUGCU_G34/"
    CCUGCU_structure=Post_Processor(CCUGCU,"CCUGCU")
    CCUGCU_structure.Fingerprint()





if __name__ != "__main__":

    restrained_residue_list = [
    1, 6, 7,
    8, 9, 10, 11, 12, 17, 18, 19, 20, 21, 22, 32, 33,
    34, 35, 36,
    37, 38, 39, 40, 41, 42, 43, 44, 45, 52, 53,
    54, 55, 56, 57, 58, 59, 60, 61, 62, 64, 74, 77, 78,
    79, 80, 81, 82,
    83, 84, 85, 86, 87, 105, 106, 107, 108, 109, 110, 111,
    112, 113, 114, 115,
    116, 117, 118,
    119, 120, 121, 142, 143, 144, 145, 146, 147, 148, 149,
    150, 151, 152, 153, 154,
    155, 156, 157, 158, 159, 160, 161, 162, 163, 176, 177, 178, 179, 180, 181,
    182, 183, 196, 199, 200, 201, 202,
    203, 204, 205, 206,
    207, 208,
    209, 210, 211, 212, 213, 214, 215, 224, 225, 226, 227,
    228, 229, 247, 248, 249, 250, 251, 252, 253,
    254, 255, 270, 271,
    272, 273, 274, 275, 276, 277,
    278, 279, 280, 281, 282, 283, 284,
    285, 286, 287, 290, 291, 292, 293, 294, 295, 296,
    297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312,
    313, 314, 315, 316,
    317, 318, 319, 331,
    332, 333, 334, 335, 336,
    337, 338, 339, 340, 341,
    342, 343, 344, 345, 363, 364, 365,
    366, 367, 368, 369, 390, 391, 392, 393,
    394, 395, 396, 397,
    398, 399, 400, 416, 417, 418,
    419, 431,
    432, 433, 434, 435,
    436, 437, 438, 439, 440, 441, 442,
    443, 453, 454, 455, 456, 457, 458,
    459, 460, 469, 470, 471,
    472, 473, 474, 475,
    476, 477, 478, 479, 490, 491,
    492, 493, 494]


    CCUGCU_G34=Run_Program('/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/BKUP_5JUP_N2_GCU/',
                             '/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop',
                             restrained_residue_list,
                             0,
                             "CCUGCU_G34")
    
    hydrogen/H_Print_CCUCGU_G34/CCUCGU.csv
    

    CCUCGU_G34=Run_Program('/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/BKUP_5JUP_N2_CGU/',
                            '/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/TLEAP/5JUP_N2_CGU_nowat.prmtop',
                            restrained_residue_list,
                            0,
                            "CCUCGU_G34")
    

    start=time.time()
    CCUGCU_G34.iterate_replicates()
    parallelized_time = time.time() - start 
    print(f"our time for loading and processing all frames of all trajectories in CCUGCU_G34 \n Time = {parallelized_time}")


    start=time.time()
    CCUCGU_G34.iterate_replicates()
    parallelized_time = time.time() - start 


    print(f"our time for loading and processing all frames of all trajectories in CCUCGU_G34 \n Time = {parallelized_time}")

    #Histogram(Trajectory_array,Final_percent_frames,final_average_Hbond)    

    #Filestore=FileManagement("/zfshomes/lperez/SUBSYSTEMS/G_34/CC/CCUGCU_G34/NEUTRAL/Trajectories/")
    

    

    
