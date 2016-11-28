#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 15:58:04 2016

@author: Jiajia
"""

import numpy as np
import random
from scipy.interpolate import interp1d
from scipy.interpolate import splev, splrep

#measured time; measured data; computed time; sample propotion; number of times of sampling
def cal_rsq_interp(mtime,mdata,ctime,prop,iter_num):  
    try:
        #intial r squares array to store results
        rsq_quadratic = []
        rsq_cubic = []
        rsq_spline_k2 = []
        rsq_spline_k3 = []
        
        #repeat sampling
        for i in range(iter_num):
            m = len(mtime)
            n = round(prop * m)
            sample_loc = np.sort(random.sample(range(m),n))
            mtime_sp = mtime[sample_loc]
            mdata_sp = mdata[sample_loc]
            
            #quadratic
            quadratic_interp = interp1d(mtime_sp,mdata_sp,kind='quadratic')
            quadratic_computes = quadratic_interp(ctime)
            
            #cubic
            cubic_interp = interp1d(mtime_sp,mdata_sp,kind='cubic')
            cubic_computes = cubic_interp(ctime)
            
            #B-Spline k=2
            bspline_k2_interp = splrep(mtime,mdata,w=None, xb=None, xe=None, k=2, task=0, s=None, t=None, full_output=0, per=0, quiet=1)
            bspline_k2_computes = splev(ctime,bspline_k2_interp)
            
            #B-Spline k=3
            bspline_k3_interp = splrep(mtime,mdata,w=None, xb=None, xe=None, k=3, task=0, s=None, t=None, full_output=0, per=0, quiet=1)
            bspline_k3_computes = splev(ctime,bspline_k3_interp)
          
            #fetching loci for non sampled data
            non_sample_loc = np.setdiff1d(range(m),sample_loc)
            
            #test set, origin data
            test_ori = mdata[non_sample_loc]
            mtime_nsp = mtime[non_sample_loc]
            
            #fetching loci for computed data
            ctime_list = ctime.tolist()
            computed_loc = []
            for i in range(len(mtime_nsp)):
                computed_loc.append(ctime_list.index(mtime_nsp[i]))
            
            #test set, computed data
            test_computed_quadratic = quadratic_computes[computed_loc]
            test_computed_cubic = cubic_computes[computed_loc]
            test_computed_bspline_k2 = bspline_k2_computes[computed_loc]
            test_computed_bspline_k3 = bspline_k3_computes[computed_loc]
            
            #calculate r square
            rsq_quadratic.append(sum(np.square(test_ori - test_computed_quadratic)))
            rsq_cubic.append(sum(np.square(test_ori - test_computed_cubic)))
            rsq_spline_k2.append(sum(np.square(test_ori - test_computed_bspline_k2)))
            rsq_spline_k3.append(sum(np.square(test_ori - test_computed_bspline_k3)))
        
        #print average output of sampling
        print('Average root square for quadratic regression: %d' % np.average(rsq_quadratic))
        print('Average root square for cubic regression: %d' % np.average(rsq_cubic))
        print('Average root square for B-Spline k2 regression: %d' % np.average(rsq_spline_k2))
        print('Average root square for B-Spline k3 regression: %d' % np.average(rsq_spline_k3))
            
    except AttributeError:
        break


