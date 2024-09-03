
#import dependencies
import numpy as np
import mdtraj as md
import pandas as pd
from time import time
from matplotlib import pyplot
import matplotlib.pyplot as plt
import numpy as np
import mdtraj as md

#Program Specific classes and files
from File_Management import FileManagement
from Trajectory import Trajectory
from Trajectory_Processor import Processor

#Example Use for a multiple replicates of a trajectory we have the following options
# -note functions can be used without the initial trajectory object 
# -however initializing with Trajectory object streamlines process

#Convenience file holds many useful lists 
from Convenience import CCU_System_directories,CCU_System_names,CCU_Topologies,sh_job_template,test_residues,test_slice_option

#define variables (for ease of readability mostly but also speeding things up w/ numpy)
directories=np.array(CCU_System_directories)
names=np.array(CCU_System_names)
CCU_Topologies=np.array(CCU_Topologies)

#Step 1 File Management - Initializing this class also includes an inherent operation based on provided attributes
# self.filelist=np.array(self.list_files_in_dir()) takes whatever directory you input and produces a list of all the files
# self.filtered_filelist=np.array(self.Filter_files()) filteres files based on input 


#parameter examples
name="CCUGCU_G34"
directory_input="/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/BKUP_5JUP_N2_GCU/"
CCUGCU_topology="/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop"
filter_parameter="1in50"

#running only thirty replicates
Files=FileManagement(directory_input=directory_input,name=name,filter_parameter=filter_parameter)
upto30=Files.toomany_reps()
#create trajectory object with replicate list in place of singular trajectory
test=Trajectory(replicates=upto30,topology=CCUGCU_topology,residues=test_residues,slice_option=test_slice_option)

#intiialize subsequent processor object with Trajectory in place
process=Processor(Trajectory=test,name=name)

#run_iterate_replicates
_30_Replicate_array=process.create_replicate_array()





