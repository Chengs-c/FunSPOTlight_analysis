import STRIDE
#from STRIDE import deconvolution
#deconvolve()
from STRIDE.Deconvolution import DeconvolveParser, Deconvolve
import os
path="D:/research/undergraduate/fda_ST/STRIDE_data_set/data2/STRIDE"
os.chdir(path)
if __name__=='__main__':
    Deconvolve(sc_count_file="sc_count.txt",
               sc_anno_file="sc_anno.txt",
               st_count_file="st_count.txt",
               sc_scale_factor=1000,
               st_scale_factor=1000,
               normalize=False,
               out_dir="result",
               out_prefix="result")
