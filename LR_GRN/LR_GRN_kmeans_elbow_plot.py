import pandas as pd
import numpy as np
import random
from scipy.spatial.distance import cdist, pdist
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

df_bna = pd.read_csv("/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Bna_exp_fpkm_filter_interp_spline_k2.csv",header='infer',index_col=0)

# random sample 1000
sample_loc = np.sort(random.sample(range(len(df_bna)),1000))
df_test = df_bna.iloc[sample_loc]

# transform to matrix
df_test_mat = df_test.as_matrix()

K = range(1,50)

KM = [KMeans(n_clusters=k).fit(df_test_mat) for k in K]
centroids = [k.cluster_centers_ for k in KM]

D_k = [cdist(df_test_mat, cent, 'euclidean') for cent in centroids]
cIdx = [np.argmin(D,axis=1) for D in D_k]
dist = [np.min(D,axis=1) for D in D_k]
avgWithinSS = [sum(d)/df_test_mat.shape[0] for d in dist]

# Total with-in sum of square
wcss = [sum(d**2) for d in dist]
tss = sum(pdist(df_test_mat)**2)/df_test_mat.shape[0]
bss = tss-wcss
perVarExp = bss/tss*100

# number of K you choose-1
kIdx = 20-1

# elbow curve
fig = plt.figure(2, figsize=(15 , 10))
ax = fig.add_subplot(111)
ax.plot(K, avgWithinSS, 'b*-')
ax.plot(K[kIdx], avgWithinSS[kIdx], marker='o', markersize=12, 
        markeredgewidth=2, markeredgecolor='r', markerfacecolor='None')
plt.grid(True)
plt.xlabel('Number of clusters', fontsize=18)
plt.ylabel('Average within-cluster sum of squares', fontsize=18)
plt.title('Elbow for KMeans clustering', fontsize=24)
fig.savefig('/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Elbow_KMeans_soq.pdf',format='pdf')

fig = plt.figure(1, figsize=(15 , 10))
ax = fig.add_subplot(111)
ax.plot(K, perVarExp, 'b*-')
ax.plot(K[kIdx], perVarExp[kIdx], marker='o', markersize=12,
        markeredgewidth=2, markeredgecolor='r', markerfacecolor='None')
plt.grid(True)
plt.xlabel('Number of clusters', fontsize=18)
plt.ylabel('Percentage of variance explained', fontsize=18)
plt.title('Elbow for KMeans clustering', fontsize=24)
fig.savefig('/Users/Jiajia/WORK/projects/LR_GN/diceseq_exp/Elbow_KMeans_pve.pdf',format='pdf')
