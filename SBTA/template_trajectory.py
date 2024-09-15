from Run_Data_Creation import data_creator
from Convenience import restrained_residue_list


outname="CCUCGU_G34"
Trajectory_or_Replicate_Path="/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/BKUP_5JUP_N2_CGU/"
Topology="/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/TLEAP/5JUP_N2_CGU_nowat.prmtop"
slice_option=0
filter_parameter="mdcrd_1in50" #filter parameter here is a string of values you might need for filtering


creationist_two = data_creator(outname,Trajectory_or_Replicate_Path,Topology,filter_parameter,Residues_of_interest=restrained_residue_list,slice_option=slice_option)
creationist_two.Full_array_routine()

#parameter examples
outname="CCUGCU_G34"
Trajectory_or_Replicate_Path="/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/BKUP_5JUP_N2_GCU/"
Topology="/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop"
slice_option=0
filter_parameter="mdcrd_1in50" #filter parameter here is a string of values you might need for filtering


creationist = data_creator(outname,Trajectory_or_Replicate_Path,Topology,filter_parameter,Residues_of_interest=restrained_residue_list,slice_option=slice_option)
creationist.Full_array_routine()


