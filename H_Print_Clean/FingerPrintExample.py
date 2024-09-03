from Run_Fingerprints import Run_Fingerprint

one="testtxt/CCUGCU_Replicate_Average.npz"
two="testtxt/CCUCGU_Replicate_Average.npz"
top="/home66/mraval/tRNAmod/34GUU36/AAUGCU_nomod/TLEAP/5JUP_N2_nowat.prmtop"
nombre="CCU_G34_nomod"
thresh=.6

test_run = Run_Fingerprint(file_one=one,file_two=two,topology=top,name=nombre,threshold=None)
test_run.create_difference_fingerprint()
