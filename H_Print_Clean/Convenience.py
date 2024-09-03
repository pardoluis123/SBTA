#For Uasy Use Inside Of Fingerprint Visualization
restrained_residue_list = (1, 6, 7, 8, 9, 10, 11, 12, 17, 18, 19, 20, 21, 22, 32, 33, 34, 35, 36,37, 38, 39, 40, 41, 42, 43, 44, 45, 52, 53,54, 55, 56, 57, 58, 59, 60, 61, 62, 64, 74, 77, 78,79, 80, 81, 82,83, 84, 85, 86, 87, 105, 106, 107, 108, 109, 110, 111,112, 113, 114, 115,116, 117, 118,119, 120, 121, 142, 143, 144, 145, 146, 147, 148, 149,150, 151, 152, 153, 154,155, 156, 157, 158, 159, 160, 161, 162, 163, 176, 177, 178, 179, 180, 181,182, 183, 196, 199, 200, 201, 202,203, 204, 205, 206,207, 208,209, 210, 211, 212, 213, 214, 215, 224, 225, 226, 227,228, 229, 247, 248, 249, 250, 251, 252, 253,254, 255, 270, 271,272, 273, 274, 275, 276, 277,278, 279, 280, 281, 282, 283, 284,285, 286, 287, 290, 291, 292, 293, 294, 295, 296,297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312,313, 314, 315, 316,317, 318, 319, 331,332, 333, 334, 335, 336,337, 338, 339, 340, 341,342, 343, 344, 345, 363, 364, 365,366, 367, 368, 369, 390, 391, 392, 393,394, 395, 396, 397,398, 399, 400, 416, 417, 418,419, 431,432, 433, 434, 435,436, 437, 438, 439, 440, 441, 442,443, 453, 454, 455, 456, 457, 458,459, 460, 469, 470, 471,472, 473, 474, 475,476, 477, 478, 479, 490, 491,492, 493, 494)
unrestrained_residue_list=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 
 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 
 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 
 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 
122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 
141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 
160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 
179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 
198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 
217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 
236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 
256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 
276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 
296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 
316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 
336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 
356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 
376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 
396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 
416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 
436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 
456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 
476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493)
n2_residue_numbers = {
    94: 'C of Car',
    127: 'A of Car', 
    240: 'R of Car',
    408: 'tRNA Wobble',
    409: 'tRNA nt 2',
    410: 'tRNA nt 1',
    423: 'mRNA nt 1',
    424: 'mRNA nt 2',
    425: 'mRNA Wobble',
    426: '+1 Codon n1',
    427: '+1 Codon n2',
    428: '+1 Codon n3',
}
CAR_plus1_Codon = {
                    94:'E. coli 16S',
                    127:'E. coli 16S', 
                    240:'S. cerevisiae rps3',
                    408:'Trna Wobble',#bc of schematic
                    426:'+1codon n1',
                    427:'+1codon n2',
                    428:'+1codon n3',
                    }
p_site = [411,412,413,410,409,408,94,127,240,420,421,422,423,424,425,426,427,428]

#Lists of directories with their names and topologies ordered as such
CCU_System_directories=["/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/BKUP_5JUP_N2_GCU/",
                    "/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/BKUP_5JUP_N2_CGU/",
                    "/zfshomes/lperez/SUBSYSTEMS/I_34/CCUCGU_I34/ANALYSIS/BACKUP/",
                    "/zfshomes/lperez/SUBSYSTEMS/I_34/CCUGCU_I34/ANALYSIS/BACKUP/"]

CCU_System_names=["CCU_GCU_G34",
              "CCU_CGU_G34",
              "CCUCGU_I34",#thirty
              "CCUGCU_I34"]#thirty

CCU_Topologies=[
    "/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop",
    "/home66/kscopino/AMBER22/CODONS/CCUCGU_G34/TLEAP/5JUP_N2_CGU_nowat.prmtop",
    "/zfshomes/lperez/SUBSYSTEMS/I_34/CCUCGU_I34/TLEAP/5JUP_N2_CCUCGU_I34_nowat.prmtop",
    "/zfshomes/lperez/SUBSYSTEMS/I_34/CCUGCU_I34/TLEAP/5JUP_N2_CCUGCU_I34_nowat.prmtop"]

#For Testing Residue Creation
test_trajectory="testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd"
test_topology="testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop"
test_residues=restrained_residue_list
test_slice_option=0
#For Testing Residue Creation

#Mitsu's AAU files 
AAU_GCU_NoMod="/home66/mraval/tRNAmod/34GUU36/AAUGCU_nomod/5JUP_N2_2fpns.mdcrd"
AAU_CGU_NoMod="/home66/mraval/tRNAmod/34GUU36/AAUCGU_nomod/5JUP_N2_2fpns.mdcrd"
AAU_GCU_Top="/home66/mraval/tRNAmod/34GUU36/AAUGCU_nomod/TLEAP/5JUP_N2_nowat.prmtop"
AAU_CGU_Top="/home66/mraval/tRNAmod/34GUU36/AAUCGU_nomod/TLEAP/5JUP_N2_nowat.prmtop"

#For writing outfiles
sh_job_template='''
#SBATCH -N 1
#SBATCH --partition=mwgpu
#SBATCH --mem=30G
#SBATCH -B 1:1:1
#SBATCH --cpus-per-gpu=1

source /zfshomes/lperez/miniconda3/etc/profile.d/conda.sh
conda activate


'''
#For importing into the program
imports_for_program='''
#import dependencies
import os
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

#Convenience file holds many useful lists 
from Convenience import restrained_residue_list


#We begin with a multi-system run of the program
    # -note functions can be used without the initial trajectory object 
    # -however initializing with Trajectory object streamlines process
'''

