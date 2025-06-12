#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 12:08:53 2025

combine the output of the grid calculation to make the data file for interpolation

@author: boeglinw
"""

import LT.box as B
import numpy as np
import os
import sys
import glob as G

data_dir = './calc_grid_all/grid/'

icon = 12
iw = 3
noff = 1
npv = 1
ics = 1

file_patt = f'csec_calc_ThQ*_{icon}_{iw}_{noff}_{npv}_{ics}.data' # select files with the correct calc. type

files = G.glob(data_dir + file_patt)


data = []
for f in files:
    d = B.get_file(f)
    th = d['thr0'][0]
    q2 = d['q2'][0]
    data.append([th, q2, d])
#%%
data_a = np.array(data, dtype = object)
th_a = data_a[:,0]
q2_a = data_a[:,1]


i_s_th = np.argsort(th_a)


#%%
u_v, u_i, u_c = np.unique(data_a[:,0][i_s_th], return_index = True, return_counts = True)  # get unique values

th_slices = [slice(*zz) for zz in zip(u_i, u_i + u_c)]  # find the slices with constant theta (data_a[:,0]) values

#%% make a new list of arrays so that the th values are sorted and grouped in groups
i_s_all = []
for sl in th_slices:
    i_ss = i_s_th[sl] 
    i_s_q2 = np.argsort(q2_a[i_ss])  # sorte values of constant theta according to q2
    i_s_th_q2 = i_ss[i_s_q2]
    i_s_all += list(i_s_th_q2)

# i_s_all array of indices into data_a with sorted theta and sorted q2 _values
i_s_all = np.array(i_s_all) 
#%% loop over sorted data and write output file

N_q2 = 60
N_pr = 50
N_thr = 60

grid_name = f'strfun_grid_{icon}_{iw}_{N_q2}_{N_pr}_{N_thr}.data'
o = open(grid_name,'w')
for th, q2, d in data_a[i_s_all]:
    n_dat = d['pr'].shape[0]
    o.write(f"{n_dat}\n")
    o.write(f"{q2} {th}\n")
    o.write(f" 1 icon={d.par['icon']}/iw={d.par['iw']}\n")
    for i,pr in enumerate(d['pr']):
        o.write(f"{pr:.12g} {d['w_l'][i]:.12g} {d['w_t'][i]:.12g} {d['w_lt'][i]:.12g} {d['w_tt'][i]:.12g}\n")
    # next data file

    
     
