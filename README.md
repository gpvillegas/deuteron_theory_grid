Description:
------------
This code uses Misak Sargsians edepngo4 cross section code to create grids in Q2, pr, and thr. 
This was designed for deuteron break up SF and xsec calculations. 

Code setup:
-----------
The cross section code lives in the directory 'edepngo4_v1/', this is fortran code that needs to be compiled before using.
To do this run these commands in the main directory:

	cd edepngo4_v1/
 	make

How to calculate a grid:
------------------------
Important: the shell scripts created are bash scripts

1. The following scripts have to be run in order if you want to to create a new grid:

- **prep_grid.py** : prepares input files for the calculation
Here you can select the specific parameters for your desired grid, like

     	icon -> (line 83) type of model: 1 = PWIA, 2 = forward rescattering, 3 = charge exchange rescattering, 12 = PWIA+FSI, 13 =  PWIA+CHEX, 123 = PWIA+FSI+CHEX
     	iw -> (line 103) wave function model: 1 = Paris, 2 = V18, 3 = CD-Bonn, 4 = AV18sb
  
	Then choose your kinematic range in pr, q2 and thr. Important, the # of intervals will determine the max size of your grid, most likely the resulting grid will be smaller than this number due to forbidden kinmatics or other calculation failures that elimiante some points.

- **prep_calc_all_crosec_KG.py** : prepares new directory with scripts and linbks to run edpngo4_v1
Here you can set up the number of processes (number of xterms) you want to split up the job into. Here are some stats for some grids I have made. Generally, PWIA runs much faster than FSI (duh)
     for a PWIA grid with ~ 150000 points (3721 files) in 50 terminals it lasted ~ 1h
	for the same grid with FSI (3721 files in 50 terminals) it lasted ~ 6d (yep, run at your own risk)

In this script setup the following BEFORE you run it:
There are 3 directory locations you have to set up, right now these are set up to my local directories: (lines 58-62)
      
      BASE_DIR = '/home/gvill/deuteron/deuteron_theory'

      EDPNGO_DIR = '/edepngo4_v1'

      LOC_DIR = '/grid' # get the current directory

After running this script, skip to point 2.

- calc_crosecs.py : steering script to run edpngo4_v1 for a set of kinenatics in a file
Contains the names the resulting grid will have (line 245)

  		result_file = result_dir + f"{output_header}{k_name}_{icon}_{iw}_{noff}_{npv}_{ics}.data"

This name is for the files in the kin_grid/ directory.
Refer to any of the grid files headers to know the meaning of the control parameters icon, iw, noff, npv, ics.

- make_grid_file.py: combines the output of the full calculation to create an ascii input file used in the interpolation code
Here you have to point to the grid files directory (line 17)

     	data_dir = './calc_grid_all/grid/'

And also modify the file pattern parameters to match the grid you are making (line 25)

    file_patt = f'csec_calc_ThQ*_{icon}_{iw}_{noff}_{npv}_{ics}.data'
where the meaning of the pattern is in calc_crossecs.py above.

Also modify the output grid name pattern to your liking, right now is (line 66):

    grid_name = f'strfun_grid_{icon}_{iw}_{N_q2}_{N_pr}_{N_thr}.data'
    
you can get N_q2, N_pr, N_thr from the initial prep_grid.py you ran.

2. Once the prep_calc_all_crosec_KG.py has been run, cd to the calculation directory e.g. ./calc_grid_all and 
run the script ./run_all.sh  this will perform the entire calculation.
After the calculation is complete then use make_grid_file.py to make the data file used
in interpolation.                 

_Optional:_

- **kill_xterms.py** : script to kill all x-terms created once calculation is complete (this needs to be adjusted to the OS used)
- **clean_grid.py**  : deletes all files in ./kin_grid and ./calc_grid_all, I recommend doing this to avoid confusion with the kin_grid files.
