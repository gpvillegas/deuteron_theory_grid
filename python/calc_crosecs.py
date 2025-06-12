"""
 calculate cross sections using Misaks code

 compile and link code: ./run_sigma_ed4.mak

 to try out: 1) make a subdirectory in edepngo4 called calc and cd into it
             2) make sure your working directory in spyder is now calc
             3) run the following command:
                    %run ../python/calc_crosecs.py -S ../python/kinematics_list.data
                    
                    this should setup the necessary links and run the code
                    
             4) check the other command line options by doing
                
                    %run ../python/calc_crosecs.py -h
                    

             5) loopk at the kinbematics file in the python subdirectory: kinematics_list.data 
                    
    
"""
#
#
import LT.box as B
import numpy as np
import run_command as rc
import os
import sys
import time

# argument parser
import argparse as AG


#%% Default values
# These are the detault values


# all path are relative to this one
BASE_DIR = '/home/gvill/deuteron/deuteron_theory'

# location of edpngo relative to BASE_DIR: localtion of the program
EDPNGO_DIR = '/edepngo4_v1'

# local directory relative to BASE_DIR: where he calculation is carried out
LOC_DIR = '/grid'

# result directory relative to LOC_DIR
RES_DIR = '/results'

# program name
PROGRAM = '/run_sigma_ed4'


#%% helper functions
def check_dir(out_dir):
    # create output directory (if necessary)
    try:
        os.mkdir(out_dir)
    except:
        msg = sys.exc_info()[1]
        if msg.errno == 17 :
            print(out_dir, " exists, will use it ")
        else:
            print("problem : ", msg)
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


def setup(dir_name):
    # setup the links for the directory dir_name
    
    src1 = BASE_DIR + EDPNGO_DIR + '/B_ampl_pn.dat'; dst1 = '/B_ampl_pn.dat'
    src2 = BASE_DIR + EDPNGO_DIR + '/B_ampl_pp.dat'; dst2 = '/B_ampl_pp.dat'
    src3 = BASE_DIR + EDPNGO_DIR + '/B_ampl_pp_nocoul.dat'; dst3 = '/B_ampl_pp_nocoul.dat'    
    
    # setup the symbolic links to the data files
    make_link(src1, dir_name + dst1)
    make_link(src2, dir_name + dst2)
    make_link(src3, dir_name + dst3)
    
    return
    
#%%----------------------------------------------------------------------
# conversion factors
mev2gev = 1.e-3  
dtr = np.pi/180.

comment1 = 'nkin, iw, ics'
comment2 = 'ei, pr, thr0, phir0, ktf, noff, npv, npro, icon, isd, ipp, iv1, iv2, q2 , ID'



#%% setup command line arguments

parser = AG.ArgumentParser(prog = '/calc_crosecs', 
                           description = 'Calculate cross sections using sigma_ed4',
                           formatter_class=AG.ArgumentDefaultsHelpFormatter)
# nargs = '?' use 1 argument and if missing assign defaul value
#parser.add_argument("file", nargs='?', help="my help text for this", default = 'data.dat')
# nargs = '*' use all positional argumensts e.g. wildcards
parser.add_argument("kin_file", nargs = '?', help="kinematics file name", default = 'kinematics_list.data')

parser.add_argument('-o', '--output_header', help="header for output file", default = 'csec_calc_')
parser.add_argument('-S', '--run_setup', help="run setup to make sure the links are correct", action='store_true')
parser.add_argument('-D', '--base_dir', help="set the base directory", default = BASE_DIR )
parser.add_argument('-E', '--edpngo_dir', help="edpngo director relative to base", default = EDPNGO_DIR )
parser.add_argument('-L', '--local_dir', help="set the local (working) directory relative to base", default = LOC_DIR )
parser.add_argument('-R', '--result_dir', help="set the result directory relative to local", default = RES_DIR )
parser.add_argument('-P', '--program_name', help="program name", default = PROGRAM )

args = parser.parse_args()


result_dir = args.result_dir

# store the new values as defaults
RES_DIR = args.result_dir
BASE_DIR = args.base_dir
LOC_DIR = args.local_dir
EDPNGO_DIR = args.edpngo_dir

kin_file = args.kin_file

output_header = args.output_header


#%%----------------------------------------------------------------------
# setup the directory if selected

# make the links to the data files needed


if args.run_setup:
    setup(BASE_DIR + LOC_DIR)
    # command for running the code

    command = BASE_DIR + EDPNGO_DIR + PROGRAM



    # setup directoes for log files and error logs
    log_dir = BASE_DIR + LOC_DIR + '/log'
    err_dir = BASE_DIR + LOC_DIR + '/err'
    make_link(command, BASE_DIR + LOC_DIR + PROGRAM)

    # directory for the results
    result_dir = BASE_DIR + LOC_DIR + RES_DIR

    # make sure the directories exist
    check_dir(log_dir)
    check_dir(err_dir)
    check_dir(result_dir)
    
    print('Setup complete!')
    sys.exit(0)

#%%----------------------------------------------------------------------
# command for running the code

command = BASE_DIR + EDPNGO_DIR  + PROGRAM


# setup directoes for log files and error logs
log_dir = BASE_DIR + LOC_DIR + '/log'
err_dir = BASE_DIR + LOC_DIR + '/err'

# directory for the results
result_dir = BASE_DIR + LOC_DIR + RES_DIR

# make sure the directories exist
check_dir(log_dir)
check_dir(err_dir)
check_dir(result_dir)

 
# open kinematic file (given as argument)
kin_data = B.get_file(kin_file)

k_name = kin_data.par['kinematics'] 

# get the list of files for different theta values available

# kinematics file
nkin = len(kin_data)
if nkin == 0:
    print(70*'-')
    print('nothing to calculate for ', kin_data.filename)
    print(70*'-')
    sys.exit(0)
print(f'will calculate {nkin} kinematic settings ')
"""
# loop over all kinematic variables and calculate cross sections
"""


# create input file for run_sigma_ep

# get calculation control parameters from parameter section of datafile

ktf = kin_data.par['ktf']
noff= kin_data.par['noff']
npv= kin_data.par['npv']
npro= kin_data.par['npro']
icon = kin_data.par['icon']
iw = kin_data.par['iw']
isd = kin_data.par['isd']
ipp = kin_data.par['ipp']
iv1 = kin_data.par['iv1']
iv2 = kin_data.par['iv2']
ics = kin_data.par['ics']

# use user comments if present
try:
    comment1 = kin_data.par['comment1']
except:
    pass

try:
    comment2 = kin_data.par['comment2']
except:
    pass


# creat output file name
# result_file = result_dir + "csec_calc_%s_%d_%d_%d_%d.data"%(kin,iw, noff, npv, icon)
 
result_file = result_dir + f"{output_header}{k_name}_{icon}_{iw}_{noff}_{npv}_{ics}.data"


# prepare kinematics file for calculation
o =open('sigma_ed.kin','w')
o.write(comment1 + '\n')
o.write( "%d  %d %d\n"%(nkin, iw, ics) )
o.write(comment2 + '\n')
# now loop over the good kinematics
for k in kin_data:
    # prepare input file  for 
    # read(15,*) ei, pr, thr0, phir0, ktf, noff, npv, npro, icon, isd, ipp, iv1, iv2, q2 , ID
    # phir is phi of the recoiling system
    
    ol = f"{k['E_i']:g} {k['p_miss']:g}  {k['theta_r']:g} {k['phi_r']:g} {ktf:d} {noff:d} {npv:d} {npro:d} {icon:d} {isd:d} {ipp:d} {iv1:d} {iv2:d} {k['Q2']:g} {k['ID']:s}\n"
    o.write(ol)
o.close()
# now run the code
outfile = log_dir +  '/calc_output_'+k_name+'.dat'
errfile = err_dir + '/calc_err_'+k_name+'.dat'

# run the code
print("running : ", command)
t_start = time.time()
rc.run_command(command, outfile, errfile)
t_end = time.time()
print(f"finished, writing results, time used : {t_end - t_start:.2f} s")

# check if there were errors (error file is not empty)
if os.stat(errfile).st_size == 0:
    print('No errors logged')
else:
    print(70*'-')
    print(f'Warning: errors logged in {errfile}')
    print(70*'-')
    sys.exit(-1)

# prepare output file
# write header information
o = open(result_file,'w')
l = "# %d  %s\n"%(ktf,   "# ktf the knock-out nucleon is proton(1), neutron(-1)")
o.write(l)
l = "# %d  %s\n"%(noff,  "# noff Off-Shell effects included (1) not included (0)")
o.write(l)
l = "# %d  %s\n"%(npv ,  "# npv only pole term in FSI (0), pole+PV (1) ")
o.write(l)
l = "# %d  %s\n"%(npro,  "# npro it calculates d\sigma/dE'_e d\Omega_e d\Omega_pf")
o.write(l)
l = "# %d  %s\n"%(icon,  "# icon = 0 initialization, 1 - PWIA, 2 -forward rescattering, 3 - charge exchange rescattering ")
o.write(l)
l = "#         12 - PWIA+FSI, 13 - PWIA+CHEX, 123 - PWIA+FSI+CHEX\n"
o.write(l)
l = "# %d  %s\n"%(iv2, " #i v2  - version of parameterization of CHEX amplitude - ")
o.write(l)
l = "#         10 - Gibbs-Loiseau (has only real part), 21 - SAID, 1021 - SAID+GL\n"
o.write(l)
l = "#         10021 - SAID+GL (with Im part at p>p0 modelled using SAID's value\n"
o.write(l)
l = "#         -21 SAID   - For details look at fpn_chex.f\n"
o.write(l)
l = "# %d  %s\n"%(iv1 ,   "# iv1  - version of FSI; 1- diagonal FSI using fpn.f, 2- diagonal FSI using f_NN_sd.f ")
o.write(l)
l = "#         which uses SAID parameterization for plab<3.68, 3 - nondiagonal approximation \n"
o.write(l)
l = "#         using f_NN_sd.f. At p>1.4 it uses diffractive parametrization like fpn.f\n"
o.write(l)
l = "# %d  %s\n"%(isd  , "# isd  - isd = 1 spectator mechanism, isd=2 direct mechanism isd =0 both   ")
o.write(l)
l = "# %d  %s\n"%(iw  , "# iw   - 1 - Paris 2 - V18 deuteron wave funcion, 3- cd Bonn, 4 - AV18sb")
o.write(l)
l = "# %d  %s\n"%(ipp ,  "# ipp  - 1 - quasi \"pp\" wave functio - only S-parital wave, 0, both S and D waves")
o.write(l)
l = "# %d  %s\n"%(ics ,  "# ics  - 1 - SLAC, 2 Kelly, 3 - Bodek, BB Arrington for-factor paremeterizations")
o.write(l)
l = "\n#    crs  - cross section in nb/GeV/str^2\n\n"
o.write(l)
l = f"#\ icon = {icon}; iw = {iw}\n"
o.write(l)
# read the result file of the calculation
sigma_ed4 = B.get_file('sigma_ed4.data')
# now write all to the result_file
sigma_ed4.write_all(o)
o.close()
# all done
