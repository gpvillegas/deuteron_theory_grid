#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 16:29:46 2025

# make grid calculation input files


# Grid is in pr, q2, thr

@author: boeglinw
"""

import numpy as np
import LT.box as B

import os
import sys

# general setup
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
        

#%%
def set_kin(kin):
    kin_info = f"""
# ----------------------------------------------------------------------
# control parameters for edpngo4_v1 calculation

# name of kinematics (no blanks)

# ----------------------------------------------------------------------
#\kinematics = {kin}
# ----------------------------------------------------------------------
"""
    return kin_info
    

#%% set comment

def set_comment(com1, com2):
    comment = f"""
# optional comments inserted into kinematics file for code (not important)

# \ comment1 = {com1}

# \ comment2 = {com2}

# ----------------------------------------------------------------------
# calculation control
# ---------------------------------------------------------------------- 
"""
    return comment
    


#%%

parameters = """

#   ktf: the knock-out nucleon is proton(1), neutron(-1)
#\   ktf  = 1  

#  noff: Off-Shell effects included (1) not included (0)
#\   noff = 1  

#  npv: only pole term in FSI (0), pole+PV (1)
#\   npv  = 1

#  npro: it calculates d\sigma/dE'_e d\Omega_e d\Omega_pf
#\   npro = 1  

#  icon: icon = 0 initialization, 1 - PWIA, 2 -forward rescattering, 3 - charge exchange rescattering
#         12 - PWIA+FSI, 13 - PWIA+CHEX, 123 - PWIA+FSI+CHEX
#\   icon = 12


#  iv2:  - version of parameterization of CHEX amplitude -
#         10 - Gibbs-Loiseau (has only real part), 21 - SAID, 1021 - SAID+GL
#         10021 - SAID+GL (with Im part at p>p0 modelled using SAID's value
#         -21 SAID   - For details look at fpn_chex.f          
#\   iv2 = 10  


#  iv1:  - version of FSI; 1- diagonal FSI using fpn.f, 2- diagonal FSI using f_NN_sd.f
#         which uses SAID parameterization for plab<3.68, 3 - nondiagonal approximation
#         using f_NN_sd.f. At p>1.4 it uses diffractive parametrization like fpn.f
#\  iv1= 3    

#  isd:   - isd = 1 spectator mechanism, isd=2 direct mechanism isd =0 both
#\  isd = 1   

#  iw:    - 1 - Paris 2 - V18 deuteron wave funcion, 3- cd Bonn, 4 - AV18sb
#\  iw = 3    

#  ipp:   - 1 - quasi "pp" wave functio - only S-parital wave, 0, both S and D waves
#\  ipp = 0   

#  ics:  - 1 - SLAC, 2 Kelly, 3 - Bodek, BB Arrington for-factor paremeterizations
#\   ics = 1

#         crs  - cross section in nb/GeV/str^2

#----------------------------------------------------------------------
# kinematic settings
#----------------------------------------------------------------------

# theta_r: angle between recoil and q (degrees)
# phi_r: angle of the recoil plane and the electron scattering plane (degrees) 
# p_miss: missing momentum (GeV/c)
# Q2:  - 4-momemtum transfer squared (GeV/c)**2
# E_I: incident energy (GeV)
# ID: kinematics identification string (no blanks) 

"""

header = '#! theta_r[f,0]/ phi_r[f,1]/ p_miss[f,2]/ Q2[f,3]/ E_i[f,4]/  ID[s,5]/' 

#%% incident energy: only needed for code
Ei = 10.542

phr = 0.

# location of kinematics files
kin_dir = './kin_grid/'

check_dir(kin_dir)

# Q2 range
q2_min = 2.0
q2_max = 8.0

N_q2 = 60 # number of intervals

q2_step = (q2_max - q2_min)/N_q2

# pr range
pr_min = 0.
pr_max = 1.5

N_pr = 50

pr_step = (pr_max - pr_min)/N_pr

# thr range
thr_min = 0.
thr_max = 180.

N_thr = 60

thr_step = (thr_max - thr_min)/N_thr

for i_thr in range(N_thr + 1):
    thr = thr_min + i_thr * thr_step 
    for i_q2 in range(N_q2 + 1):
        q2 = q2_min + i_q2 * q2_step
        f_name = kin_dir+ f'/{thr:.2f}_{q2:.3f}.data'
        o = open(f_name, 'w')
        o.write(set_kin(f'ThQ_{thr:.2f}_{q2:.3f}'))
        o.write(set_comment(f'{q2_min} < Q2 < {q2_max}',f'{pr_min} < pr < {pr_max}'))
        o.write(parameters)
        o.write(f'\n{header}\n')
        for i_pr in range(N_pr + 1):
            pr = pr_min + i_pr * pr_step
            o.write(f'{thr} {phr} {pr} {q2} {Ei}  k_{i_thr}_{i_q2}_{i_pr}\n')

o.close()
