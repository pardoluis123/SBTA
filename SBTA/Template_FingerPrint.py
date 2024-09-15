from Run_Fingerprints import Run_Fingerprint
from Convenience import CCU_Topologies

one="/zfshomes/lperez/fingerprint/H_Print/test_frame_GCU_Replicate_Average.npz"
two="/zfshomes/lperez/fingerprint/H_Print/test_frame_Replicate_Average.npz"
top="/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop"
nombre="CCU_CG_G34_CPPTRAJ"
thresh=0

test_run = Run_Fingerprint(file_one=one,file_two=two,topology=top,name=nombre,threshold=None)
test_run.create_difference_fingerprint()
