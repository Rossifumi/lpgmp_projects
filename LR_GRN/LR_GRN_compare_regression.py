#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 14:29:16 2016

@author: Jiajia
"""

import numpy as np
import pandas as pd
import random
import cal_rsq_interp

# random pick n genes for testing
df_bna = pd.read_csv("/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Bna_exp_fpkm.csv",header=None,index_col=0)

df_bna.columns = ["0h","0.5h","1h","2h","3h","4h","5h","6h","8h","10h","12h","14h","16h","20h"
                  ,"24h","28h","32h","36h","40h","44h","48h","56h","64h","72h"]
            
df_bna_filter = df_bna[(df_bna.sum(axis=1)>24) & (df_bna.sum(axis=1)<48000)]
          
n = 100    
m = len(df_bna_filter)

sample_loc = np.sort(random.sample(range(m),n))

df_test = df_bna_filter.iloc[sample_loc]

# time points                       
measured_time = np.array([0,0.5,1,2,3,4,5,6,8,10,12,14,16,20,24,28,32,36,40,44,48,56,64,72])
computed_time = np.arange(0,72.5,0.5)

# calculate r square
rsq_quadratic_all = []
rsq_cubic_all = []
rsq_spline_k2_all = []
rsq_spline_k3_all = [] 

prop = 0.8
iter_num = 10

for index, row in df_test.iterrows():
    measures = np.array(row)
    [quadratic_avg, cubic_avg, spline_k2_avg, spline_k3_avg] = cal_rsq_interp.cal_rsq_interp(measured_time, measures, computed_time, prop, iter_num)
    rsq_quadratic_all.append(quadratic_avg)
    rsq_cubic_all.append(cubic_avg)
    rsq_spline_k2_all.append(spline_k2_avg)
    rsq_spline_k3_all.append(spline_k3_avg)
    
#print average output of sampling
print('Average root square for quadratic regression with %d genes, bootstrap rate of %2.f and %d times iterations: %d' % (n, prop, iter_num, np.average(rsq_quadratic_all)))
print('Average root square for cubic regression with %d genes, bootstrap rate of %2.f and %d times iterations: %d' % (n, prop, iter_num, np.average(rsq_cubic_all)))
print('Average root square for B-Spline k2 regression with %d genes, bootstrap rate of %2.f and %d times iterations: %d' % (n, prop, iter_num, np.average(rsq_spline_k2_all)))
print('Average root square for B-Spline k3 regression with %d genes, bootstrap rate of %2.f and %d times iterations: %d' % (n, prop, iter_num, np.average(rsq_spline_k3_all)))   



