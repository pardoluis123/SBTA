from Trajectory import Trajectory
from Trajectory_Processor import Processor
from File_Management import FileManagement
from Convenience import restrained_residue_list
import os
import numpy as np

class data_creator():
    def __init__(self,outname=None,Trajectory_or_Replicate_Path=None,Topology=None,filter_parameter=None,Residues_of_interest=None,slice_option=None):
        self.outname=outname if outname is not None else None
        self.path=Trajectory_or_Replicate_Path if Trajectory_or_Replicate_Path is not None else None
        self.Topology=Topology if Topology is not None else None
        self.filtparm=filter_parameter if filter_parameter is not None else None
        self.res=Residues_of_interest if Residues_of_interest is not None else None
        self.slice_option=slice_option if slice_option is not None else None

    def File_management_routine(self,path=None, name=None, filter_parameter=None):
        path = path if path is not None else self.path
        name = name if name is not None else self.outname
        filter_parameter = filter_parameter if filter_parameter is not None else self.filtparm

        #if we are dealing with a directory of replicates
        if os.path.isdir(path):
            Files_Manager=FileManagement(directory_input=path,name=name,filter_parameter=filter_parameter)
            Files_Manager.create_symbolic_directory()
            if len(Files_Manager.filelist)>30:
                filelist=Files_Manager.toomany_reps()
            else:
                filelist=Files_Manager.Filter_files()
        if not os.path.isdir(path):
            filelist=path

        return filelist
    
    def Trajectory_object_routine(self,traj=None,top=None,res=None,slice_option=None):
        traj=traj if traj is not None else self.File_management_routine()
        top = top if top is not None else self.Topology
        res = res if res is not None else self.res
        slice_option  = slice_option if slice_option is not None else self.slice_option

        if isinstance(traj, np.ndarray):
            traj_object=Trajectory(replicates=traj,topology=top,residues=res,slice_option=slice_option)
        if not isinstance(traj, np.ndarray):
            traj_object=Trajectory(trajectory=traj,topology=top,residues=res,slice_option=slice_option)
        return traj_object

    def Processor_object_routine(self,traj_object=None,outname=None,datatype=0):
        traj_object = traj_object if traj_object is not None else self.Trajectory_object_routine()
        outname=outname if outname is not None else self.outname

        processor_object=Processor(Trajectory=traj_object,name=outname,datatype=0)
        
        return processor_object

    def Final_Concatenated_routine(self,proc_object=None):
        proc_object=proc_object if proc_object is not None else self.Processor_object_routine()

        if proc_object.Trajectory.replicates is not None:
            proc_object.Save_Average_replicates()

        else:
            proc_object.Save_Average_Trajectory()

        return


if __name__ == "__main__":
    #parameter examples
    outname="CCUGCU_G34"
    Trajectory_or_Replicate_Path="/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/BKUP_5JUP_N2_GCU/"
    #Trajectory_or_Replicate_Path="/zfshomes/lperez/fingerprint/H_Print/symlink/CCUGCU_G34/mdcrd_1in50_9.mdcrd" #path to trajectory file
    Topology = "/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop" #Topology file for trajectory
    filter_parameter="1in50" #filter parameter here is a string of values you might need for filtering
    Residues_of_interest=restrained_residue_list #list of residues to filter out
    slice_option=0 #slice option is 0 = delete(residue list above) or 1 = keep(keep only residue list above)
 
    test_run = data_creator(outname,Trajectory_or_Replicate_Path,Topology,filter_parameter,Residues_of_interest=restrained_residue_list,slice_option=0)
    test_run.Final_Concatenated_routine()
       





