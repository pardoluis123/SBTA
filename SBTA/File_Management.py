import os
import numpy as np
from Convenience import unrestrained_residue_list,sh_job_template,imports_for_program

class FileManagement:
    
    def __init__(self,directory_input,name,topologies_input=None,filter_parameter=None,working_directory=None,sh_job_template=None,residues=None):
        """
        A class for various filemanagement tasks including symbolic links, and filtering 

        Parameters
        ----------
        directory_input : str
            The directory from which data will be coming
        
        name : str
            basename for subdirectory creation
    
        res_filtered : str
            List of residues filtered from trajectory previously (soon to be edited out but for now important to function)
        
        residues_to_filter : str
            List of residues IN *CURRENT ARRAY* that one wishes to delete
        
        threshold : float
            An optional value for setting threshholds of significant "average" interactions
            between residues in array. 
        
        """
        self.directory_input=directory_input
        self.topologies_input=topologies_input if topologies_input is not None else None
        self.working_directory=working_directory if working_directory is not None else os.getcwd()
        self.name=name
        self.filter_parameter=filter_parameter if filter_parameter is not None else None
        self.filelist=np.array(self.list_files_in_dir())
        self.filtered_filelist=np.array(self.Filter_files()) if filter_parameter is not None else self.filelist
        self.sh_job_template=sh_job_template if not None else None
        self.sh_job_template=sh_job_template if sh_job_template is not None else None
        self.residues=residues if residues is not None else None
    
    def list_files_in_dir(self,directory=None):
        """
        returns a list of all files contained in given directory

        Parameters
        ----------
        directory_input : str
            The directory from which data will be coming
        
        """
       
        directory=directory if directory is not None else self.directory_input

        try: #quick clause in here for cases where we wanna use just the directory given if its pointing to a file
            filelist=os.listdir(self.directory_input)
            filelist= [directory+file for file in filelist]
        except:
            filelist=[]

        return filelist
    
    def Filter_files(self,filelist=None,parameter=None):
        """
        returns a filtered list of files based on parameter

        Parameters
        ----------
        filelist : list
            The list of files (each a string path) which is being filtered
        parameter : str
            An identifier by which to identify the files of interest in the list
        """
 
        filelist=filelist if filelist is not None else self.filelist
        parameter=parameter if parameter is not None else self.filter_parameter

        filtered =[filtered_file for filtered_file in filelist if parameter in filtered_file]
        np.array(filtered)
        return filtered
    
    def toomany_reps(self,filelist=None,parameter=None,working_directory=None):
        """
        returns a filtered list of paths when there are more than 30 replicates per version of system

        Parameters
        ----------
        filelist : list
            The list of files (each a string path) which is being filtered
        parameter : str
            An identifier by which to identify the files of interest in the list
        working_directory:str
            A path pointing to the current working directory or directory where files would like to be created 
        """

        parameter=parameter if parameter is not None else 30
        filter_parameter = "1in50"
        working_directory=working_directory if working_directory is not None else self.working_directory
        
        name=self.name

        symlink_dir=f"{working_directory}/symlink/{name}/"
        filelist = self.list_files_in_dir(symlink_dir)
        filtered_filelist = self.Filter_files(filelist,filter_parameter)
        lessthan30_filelist = np.empty(30, dtype=object)

        i=0

        for file in filtered_filelist:
            if file[-2:].isdigit() and int(file[-2:]) <= parameter:
                lessthan30_filelist[i]=file
                i+=1

            if file[-1:].isdigit() and not file[-2:].isdigit():
                lessthan30_filelist[i]=file
                i+=1

        lessthan30_filelist=[f"{file}.mdcrd" for file in lessthan30_filelist]
        return np.array(lessthan30_filelist)

    def create_symbolic_link(self,i,symlink_path=None):
        """
        creates a symbolic link pointing at a file (useful for mdcrd's with no .mdcrd extension)

        Parameters
        ----------
        i : path
            A string path which the symbolic link will point to
        symlink_path:str
            The desired path for the new symlink pointer 
        """
        
        symlink_path=symlink_path if symlink_path is not None else os.getcwd()
        i=i if i is not None else None
    
        try:
            os.symlink(i, symlink_path)
            print(f"Created symbolic link: {i} -> {symlink_path}")
            return 

        except FileExistsError:
            print(f"symbolic link exists: {i} -> {symlink_path}")
            return 

    def create_symbolic_directory(self,Trajectories=None,working_directory=None,name=None):
        """
        Creates a directory within which to store all the new symbolic links

        Parameters
        ----------
        Trajectories : list
            A list holding paths to trajectories that need to have symbolic links created for them
        name:str
            a string name for the trajectory of analysis        
        """
 
        Trajectories=Trajectories if Trajectories is not None else self.filtered_filelist
        working_directory=working_directory if working_directory is not None else os.getcwd()
        name=name if name is not None else self.name

        for i in Trajectories:
                symlink_dir=f"{working_directory}/symlink/{name}/"
                symlink_path = os.path.join(symlink_dir, f"{os.path.basename(i)}.mdcrd")
                os.makedirs(symlink_dir, exist_ok=True)
                self.create_symbolic_link(i,symlink_path)

    def create_data_directories(self,name=None):
        """
        Creates a directory within which to store all resulting data

        Parameters
        ----------

        name:str
            a string name for the trajectory of analysis        
        """
        
            
        name=name if name is not None else self.name

        current_directory=os.getcwd() 
        core_directory=f"{current_directory}/H_Print_{self.name}"
        self.evaluate_directory(core_directory)

    def create_replicate_subdirectories(self,parent_directory=None,filelist=None):
        """
        Creates subdirectories form which to store all resulting data

        Parameters
        ----------

        Parent_Directory:str
            A main directory in which to create the subdirectories
        
        filelist:list
            Currently being very inefficiently used as a gauge of how many directories to make
               
        """
        
        filelist=filelist if filelist is not None else filelist
        parent_directory=parent_directory if parent_directory is not None else os.getcwd()

        n=1
        for i in filelist: 
            new_directory=f"{parent_directory}/H_Print_{self.name}/Replicate_{n}"
            try:
                dir=os.listdir(new_directory)
            except FileNotFoundError:
                print(f"{new_directory} not found so creating a new directory\n")
                os.mkdir(new_directory)
            n+=1
        pass

    def evaluate_directory(self,directory):
        """
        evaluates if directory exists (why? not sure, a little redundant will re-check)
        Parameters
        ----------

        Directory:str
            A full path to the directory you are evaluating
        
        filelist:list
            Currently being very inefficiently used as a gauge of how many directories to make
               
        """

        try:
            os.listdir(directory)
        except FileNotFoundError:
            print(f"{directory} not found so creating a new directory\n")
            os.mkdir(directory)

    def cpptraj_hbond(self,directory_input=None,name=None):
        directory_input = directory_input if directory_input is not None else self.directory_input
        name = name if name is not None else self.name
        
        with open(directory_input,"r") as datafile:
            for line in datafile:
                #format['#Acceptor', 'DonorH', 'Donor', 'Frames', 'Frac', 'AvgDist', 'AvgAng']
                line=line.split()
                print(line)       
        pass

    #For running multiple systems "psuedo parallel"
    def multiple_systems(self,directory_input=None,name=None,topologies=None,filter_parameter=None,working_directory=None):
        
        #A function for creating batch jobs of different systems 
        #Each follows a template but the idea is to ~eventually~ integrate a 
        #Simple assesment of the cluster nodes at a given timepoint
        #and generate different templates for which nodes would be most likely to be open etc

        #Input
            #Type->arraylike:directory_input
                #a list of directories
            #Type->arraylike:names
                #a list of names
            #Type->arraylike:topologies
                # a list of topologies
            


        directory_input = directory_input if directory_input is not None else self.directory_input
        topologies = topologies if topologies is not None else self.topologies_input
        names= names if names is not None else self.name
        working_directory = working_directory if working_directory is not None else os.getcwd()
        filter_parameter = filter_parameter if filter_parameter is not None else None

        systems=zip(directory_input=None,names=None,topologies=None,filter_parameter=None,working_directory=None)

        for directory_input,name,topology,filter_parameter,working_directory in enumerate(systems):            
        
            self.create_sh(self,replicates=directory_input,name=name,topologies=topology,filter_parameter=filter_parameter,working_directory=working_directory,sh_job_template=sh_job_template)

        pass

    #Create Python File raw charachter data for dynamically generating charachter data 
    def create_python(self,directory_input=None,name=None,topologies=None,filter_parameter=None,working_directory=None,residues=None):
        directory_input = directory_input if directory_input is not None else self.directory_input
        topologies = topologies if topologies is not None else self.topologies_input
        name = name if name is not None else self.name
        working_directory = working_directory if working_directory is not None else os.getcwd()
        filter_parameter = filter_parameter if filter_parameter is not None else None
        residues=residues if residues is not None else self.residues

        #Initiate Trajectory Object (can be found in convenience and are a 10 frame trajectory)

        adding_imports=f"{imports_for_program}\n\n"
        
        file_management=adding_imports+f"Files=FileManagement(directory_input='{directory_input}',name='{name}',filter_parameter='{filter_parameter}',working_directory='{working_directory}')\n\n"
        
        #pull filtered filenames
        if name == "CCUGCU_G34_unrestrained_1in50" or name == "CCUCGU_G34_unrestrained_1in50":
            list_of_files= file_management + f"\n\nupto30 = Files.toomany_reps()"
        else:
            list_of_files= file_management + f"\n\nupto30=Files.filtered_filelist"

        Trajectory_setup=list_of_files+f"\n\n{name}_Trajectory = Trajectory(\nreplicates=upto30,\ntopology='{topologies}',\nresidues={list(residues)})"

        Processor_Object=Trajectory_setup+f"\n\n{name}_Processor=Processor(Trajectory={name}_Trajectory,name='{name}')\n"

        Processor_iterate=Processor_Object+f"\n\nReplicate_array={name}_Processor.create_replicate_array()"

        average_array=Processor_iterate+f"\n\naverage_replicate_array={name}_Processor.create_average_replicate_array()"

        export_arrays=average_array+f"\n\n{name}_Processor.Export_Individual(data_directory='{working_directory}',name='{name}_average_array',array=average_replicate_array,filetype=0) \n\n "

        export_arrays=average_array+f"\n\n{name}_Processor.Export_Individual(data_directory='{working_directory}',name='{name}_Replicate_array',array=Replicate_array,filetype=0) \n\n "

        return export_arrays

    #Create batch file that dynamically creates its own H_Print python file and runs it 
    def create_sh(self,systems=None,name=None,topologies=None,filter_parameter=None,working_directory=None,sh_job_template=None):
        
        systems=systems if systems is not None else None
        name=name if name is not None else self.name
        topologies=topologies if topologies is not None else None
        filter_parameter=filter_parameter if filter_parameter is not None else None
        working_directory= working_directory if working_directory is not None else self.working_directory
        sh_job_template=sh_job_template if sh_job_template is not None  else self.sh_job_template
        
        sh_file_name=f"{working_directory}/{name}_1in50_unrestrained.sh"
        python_file_name=f"run_{name}_1in50_unrestrained.py"

        python_stuff = self.create_python(systems,name,topologies,filter_parameter,working_directory)
        with open(sh_file_name,'w') as outfile:

            outfile.write(sh_job_template)

            outfile.write(f"cd {working_directory}")
            outfile.write(f"\n\ncat<< EOF >{python_file_name}\n")

            outfile.write(f"\n{python_stuff}\n")
            outfile.write("\nEOF")
            outfile.write(f"\n\nchmod 777 {python_file_name} ")
            outfile.write(f"\n\npython {python_file_name}  >  {name} ")
            #does things here
            outfile.close()


        pass



if __name__=="__main__":

    NNU_GCU_G34='/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/BKUP_5JUP_N2_GCU/'

    Files=FileManagement(directory_input='/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/BKUP_5JUP_N2_GCU/',
                         name="CCUGCU_G34_unrestrained_1in50",
                         topologies_input= "/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop",
                         filter_parameter="1in50",
                         sh_job_template=sh_job_template,
                         residues=unrestrained_residue_list)
    

    Files.create_sh()
    
    testing_cpptraj=FileManagement(directory_input="/zfshomes/lperez/SUBSYSTEMS/I_34/CCCCGU/NEUTRAL/NEUTRAL_1/DATA/nhb_AVE_94_all_428_all.dat",
                   name="Testrun")
    
    testing_cpptraj.cpptraj_hbond()
#
#    #Example Uses
#
#    
#    #Initiate Filemanagement Object
#    Files=FileManagement(directory_input='/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/BKUP_5JUP_N2_GCU/',name="CCUGCU_G34",filter_parameter="1in50")
#
#
#
#    Files.create_symbolic_directory()
    
    
