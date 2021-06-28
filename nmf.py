import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as spr
from numpy import genfromtxt
from sklearn.decomposition import NMF
from sklearn.cluster import KMeans
from matplotlib.pyplot import figure
import nimfa
import csv
import copy


fname='new_nn_Freqs_matrix_Pop56_5Locus_15082017.csv'
fname_pop_names='Pop_Name_5locus_56_15082017.csv'
outfile='nmf.csv'

def normalize_dat(haplotype_freq,method):
    # normalization for the input
    # 1 is to keep the full size
    #2 is the square root
    # 3 is fully normalized to 1
    hap_freq_norm=copy.deepcopy(haplotype_freq)
    sum_col=haplotype_freq.sum(axis=0)
    if (method == 2):
        sum_col=np.sqrt(sum_col);
    if (method == 3):
        sum_col=sum_col*0+1;
    print(sum_col)
    for col in range(haplotype_freq.shape[1]):
        hap_freq_norm[:,col]=haplotype_freq[:,col]/sum_col[col]
    return hap_freq_norm
        
def nmf_comp(Num_Comp,haplotype_freq,norm_method,comp_method='NIMFA'):
    # run the nmf, we can use either nimfa or the sklearn nmf
    # the output is a dictionary with W and H
    hap_freq_norm=normalize_dat(haplotype_freq,norm_method)
    if (comp_method=='NIMFA'):
        nmf = nimfa.Nmf(hap_freq_norm, max_iter=200, rank=Num_Comp, update='euclidean', objective='fro')
        nmf_fit = nmf()
        W = nmf_fit.basis()
        H = nmf_fit.coef()
    if (comp_method=='NMF'):
        model = NMF(n_components=Num_Comp, init='random', random_state=0,solver='cd')
        W = model.fit_transform(hap_freq_norm)
        H = model.components_
    return {'W': W, 'H': H}


def plot_nmf(nmf1,name_list,clustered=True):
    # plot W to visualize the resutls
    W=nmf1['W']
    H=nmf1['H']
    width = 0.8    # the width of the bars: can also be len(x) sequence
    N_Bars=H.shape[1]
    ind = np.arange(N_Bars)    # the x locations for the groups
    kmeans = KMeans(n_clusters=H.shape[0]).fit(H.transpose()) 
    H1=H
    ordered_list=name_list # I keep this array of names and then will update the order if needed
    if (clustered):
        labels=kmeans.labels_ 
        order=np.argsort(labels)
        H1=H1[:,order]
        ordered_list = [name_list[i] for i in order]
    for col in range(H.shape[1]):
        H1[:,col]=H1[:,col]/sum(H1[:,col])
        H1[:,col]=np.cumsum(H1[:,col])
    figure(num=None, figsize=(16, 10), dpi=80, facecolor='w', edgecolor='k')

    for row in reversed(range(H.shape[0])):
        plt.barh(ind, H1[row,:], edgecolor='white')
    plt.yticks(np.arange(H.shape[1]), ordered_list) # add ticks to plot
    # Custom x axis
    #plt.xticks(r, names)
    plt.ylabel("Populations")
 
    # Show graphic
    plt.show()

def save_nmf(nmf1,name_list,f_output):
    # save only W as CSV
    W=nmf1['W']
    H=nmf1['H']
    H1=H
    ordered_list=name_list # I keep this array of names and then will update the order if needed
    for col in range(H.shape[1]):
        H1[:,col]=H1[:,col]/sum(H1[:,col])
    np.savetxt(f_output, H1, delimiter=',')


# In[55]:


with open(fname_pop_names, 'rt') as f:
    reader = csv.reader(f)
    name_list = list(reader)
name_list1 = [val for sublist in name_list for val in sublist]

haplotype_freq=genfromtxt(fname, delimiter=',')
Num_Comp=8


# In[58]:


nmf1=nmf_comp(Num_Comp,haplotype_freq,norm_method=3,comp_method='NMF')
save_nmf(nmf1,name_list,outfile)
plot_nmf(nmf1,name_list1,clustered=True)


# In[59]:


plot_nmf(nmf1,name_list1,clustered=True)

