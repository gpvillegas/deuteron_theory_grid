#
# calculate cross sections using Misaks code
#
# prepare calculation for kinematic grid
#
#
from LT import datafile as df
import numpy as np
import shutil as sh
import os
import sys
import glob as G

from pathlib import Path


#%% helper functions
def check_dir(this_dir):
    if os.path.isdir(this_dir):
        print(this_dir, " exists, will use it ")
    else:
        # create output directory (if necessary)
        try:
            os.mkdir(this_dir)
        except Exception as err:
            print(f"problem creating directory {this_dir}: ", err)
            sys.exit()
        
#----------------------------------------------------------------------
def make_link(src,dst):
    # create symbolic links
    print(f'make link from {src} to {dst}')
    try:
        os.symlink(src,dst)
        print("created symbolic link for : ", src)
    except:
        msg = sys.exc_info()[1]
        if msg.errno == 17 :
            print(dst, " exists, will replace it ")
            os.remove(dst)
            # recusive call
            make_link(src,dst)
        else:
            print("make_link problem : ", msg, " for ", dst , "->", src)
            sys.exit()

#----------------------------------------------------------------------
mev2gev = 1.e-3
dtr = np.pi/180.
#----------------------------------------------------------------------


# number of processes, number of xterms

nproc = 50

#
BASE_DIR = '/home/gvill/deuteron/deuteron_theory'

EDPNGO_DIR = '/edepngo4_v1'

LOC_DIR = '/grid' # get the current directory

# setup file names and directories


calc_dir = '/calc_grid_all'
kin_dir = BASE_DIR + LOC_DIR + '/kin_grid'


root_dir = BASE_DIR + LOC_DIR + calc_dir

grid_dir = '/grid/'
result_dir = root_dir + grid_dir

k_files = G.glob(f'{kin_dir}/*.data')

check_dir(root_dir)

#%% split kinematics in groups
n_kin = len(k_files)

# print exp_files

# split kinematics input files over processes

kin_groups = np.array_split(np.arange(n_kin), nproc)
# print file_groups

if len(kin_groups) == 0:
    sys.exit("no kinematics!" )

#%% setup file groups


#
file_group = []
# create input files
for i,kg in enumerate(kin_groups):
    file_group.append(np.array(k_files)[kg])

#%% file groups are 

script_list = []


check_dir(result_dir)
# create shell scripts to run all
# loop over groups and setup directories
for i, fg in enumerate(file_group):
    dir_name = root_dir + '/d{:d}'.format(i)
    check_dir(dir_name)
    command = f"""
python ../../python/calc_crosecs.py \\
    --base_dir {BASE_DIR} \\
        --edpngo_dir {EDPNGO_DIR} \\
            --local_dir {LOC_DIR}{calc_dir}/d{i:d} \\
                --result_dir /..{grid_dir} -S 
"""
    o = open(dir_name + '/setup.sh','w')
    o.write(command)
    o.close()
    os.chmod(dir_name + '/setup.sh',0o755)

    # copy the input files
    kin_dir_loc = dir_name + '/kin_files'
    check_dir(kin_dir_loc)
    for f_l in fg:
        sh.copy(f_l, kin_dir_loc + '/.')
        
    # write command to run one calculation
    run_command = f"""
python ../../python/calc_crosecs.py \\
    --base_dir {BASE_DIR} \\
        --edpngo_dir {EDPNGO_DIR} \\
            --local_dir {LOC_DIR}{calc_dir}/d{i:d} \\
                --result_dir /..{grid_dir} $1 
"""
    o = open(dir_name + '/run_one.sh','w')
    o.write(run_command)
    o.close()
    os.chmod(dir_name + '/run_one.sh',0o755)
   
    # write the script to run the calculations for all files
    sh_command = f"""
for i in {kin_dir_loc}/*.data; do ./run_one.sh $i ; done
"""
    o = open(dir_name + '/run_all_kin.sh','w')
    o.write(sh_command)
    o.close()
    os.chmod(dir_name + '/run_all_kin.sh',0o755)
    #
    # prepare final scripts
    #
    script_name = f'/run_all_d{i}.sh'
    o = open(root_dir + script_name,'w')
    o.write(f'cd {dir_name}\n')
    o.write('./setup.sh\n')
    o.write('./run_all_kin.sh')
    o.close()
    os.chmod(root_dir + script_name,0o755)
    script_list.append(script_name)

#%% create script to run all as xterms
run_all_name = root_dir +'/run_all.sh' 
so = open(run_all_name ,'w')
for l in script_list:
    so.write('xterm -hold -e ./{} &\n'.format(l))
so.close()
os.chmod(run_all_name,0o755)
