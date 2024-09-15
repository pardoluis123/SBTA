from Data_Manip import Data_Manipulator
from Fingerprint import Fingerprint_maker
import scipy.sparse as sp
import mdtraj as md
from Convenience import restrained_residue_list,n2_residue_numbers

class Run_Fingerprint():

    def __init__(self,file_one=None,file_two=None,topology=None,name=None,threshold=None):
        self.file_one=self.load_file_routine(file_one) if file_one is not None else None
        self.file_two=self.load_file_routine(file_two) if file_two is not None else None
        self.top= md.load_topology(topology) if topology is not None else None
        self.name=name if name is not None else "filler name"
        self.threshold=threshold if threshold is not None else None
    
    def load_file_routine(self,file=None):

        file = file if file is not None else None
        file=sp.load_npz(file)


        return file

    def Initiate_Data_Managementobj_routine(self,file_one=None,file_two=None):
        file_one= file_one if file_one is not None else self.file_one
        file_two = file_two if file_two is not None else self.file_two

        Data_Management = Data_Manipulator(array_one=file_one,array_two=file_two,res_filtered=n2_residue_numbers,residues_to_filter=restrained_residue_list,topology=self.top,threshold=self.threshold)
        return Data_Management

    def create_difference_array(self,file_one=None,file_two=None,name=None):
        name=name if name is not None else None

        current_obj = self.Initiate_Data_Managementobj_routine()
        difference_array = current_obj.create_difference_array()
        return difference_array
    
    def create_difference_fingerprint(self,file_one=None,file_two=None):

        difference_array = self.create_difference_array()
        print(f"\nThis is our {difference_array}\n")
        diff=Data_Manipulator(array_one=self.file_one,residues_to_filter=n2_residue_numbers,res_filtered=restrained_residue_list,topology=self.top,threshold=self.threshold)
        
        final_filtered=diff.filter_all_over_diff(array=difference_array)
        print(final_filtered)
        Nozero_AAUCGU=Fingerprint_maker(array=final_filtered,residues_of_interest=None,top=self.top)

        Nozero_AAUCGU.Base_Fingerprint(self.name)
        Nozero_AAUCGU.filtered_difference_Fingerprint(self.name)


        return


if __name__ == "__main__":

    one="/zfshomes/lperez/fingerprint/H_Print/CCUGCU_I34_Replicate_Average.npz"
    two="/zfshomes/lperez/fingerprint/H_Print/CCUCGU_I34_Replicate_Average.npz"
    top="/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/TLEAP/5JUP_N2_CGU_nowat.prmtop"
    nombre="testing_new"
    threshold=0

    #test_run = Run_Fingerprint(file_one=two,file_two=one,topology=top,name=nombre,threshold)
    #test_run.create_difference_fingerprint()







    