import numpy as np
import pandas as pd
# from scipy import stats
#from scipy.interpolate import interp1d
from scipy.interpolate import splev, splrep
from scipy import io as spio

df_bna = pd.read_csv("/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Bna_exp_fpkm.csv",header=None,index_col=0)

df_bna.columns = ["0h","0.5h","1h","2h","3h","4h","5h","6h","8h","10h","12h","14h","16h","20h","24h","28h","32h","36h","40h","44h","48h","56h","64h","72h"]

# tmp = np.ravel(df_bna)
# stats.scoreatpercentile(tmp[tmp>0],99)

df_bna_filter = df_bna[(df_bna.sum(axis=1)>24) & (df_bna.sum(axis=1)<48000)]

df_bna_filter.to_csv('/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Bna_exp_fpkm_filter.csv',encoding='utf-8')

# interp1d with step 0.5h
measured_time = np.array([0,0.5,1,2,3,4,5,6,8,10,12,14,16,20,24,28,32,36,40,44,48,56,64,72])
computed_time = np.arange(0,72.5,0.5)

computed_results = pd.DataFrame(index=df_bna_filter.index, columns=computed_time)
computed_results = computed_results.fillna(0) 

for index, row in df_bna_filter.iterrows():
    measures = row
    # cubic_interp = interp1d(measured_time,measures,kind='quadratic')
    # computes = cubic_interp(computed_time)

    bspline_interp = splrep(measured_time,measures,w=None, xb=None, xe=None, k=2, task=0, s=None, t=None, full_output=0, per=0, quiet=1)
    computes = splev(computed_time,bspline_interp)
    
    computed_results.set_value(index=index,col=computed_time,value=computes)
    
# save results
spio.savemat('/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Bna_exp_fpkm.mat',{'computed_results': computed_results, 'df_bna_filter': df_bna_filter, 'df_bna': df_bna})

# write csv
computed_results.to_csv('/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Bna_exp_fpkm_filter_interp_spline_k2.csv',encoding='utf-8')

