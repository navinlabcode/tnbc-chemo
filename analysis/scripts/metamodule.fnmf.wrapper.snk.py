#!/volumes/USR1/yyan/apps/anaconda3/bin/python

import glob
import os
import numpy as np

proj = './metamodule_fnmf'
dir_in = proj + '/expr_matrix/each_sample'
f_in_list = glob.glob(dir_in + '/ARTC*.gc.NonNegCenterMat.rds')
dir_out = proj + '/nmf_each_sample'
# K=20
K_options = [i+2 for i in np.arange(0, 30, 1)] # 2,4,6,...,100
samples_list = [os.path.basename(x) for x in f_in_list]
samples_list = [x.replace('.gc.NonNegCenterMat.rds', '') for x in samples_list]
print(samples_list, flush=True)
print(len(samples_list), flush=True)

# ø dir_out
#     ø nmf_each_sample
#         ø ARTCxxx
#             ø Kxxx
#                 - reconstruction_err.txt
                
RCMD = '/usr/bin/Rscript fastnmf.R '

rule all:
    input:
        expand(dir_out + '/{sample}' + '/K' + '{K}' + '/reconstruction_err.txt', sample=samples_list, K=K_options)

rule nmf:
    input:
        mtx = dir_in + '/{sample}.gc.NonNegCenterMat.rds',
    params:
        k = '{K}'
    output:
        txt = dir_out + '/{sample}' + '/K' + '{K}' + '/reconstruction_err.txt'
    run:
        cmd = "{} {} {} {}".format(
            RCMD,
            input.mtx,
            str(params.k), 
            os.path.dirname(output.txt))
        print(cmd, flush=True)
        shell(cmd)
        
        
        
    
