import collections,pickle,ast,os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

####################################################################
"""
Helper functions
    -analyzeSensitivity

Plot Functions
    - heat
    - PCA_plot
    - trajectory (pick one that looks converged)
    - plotBaseCaseIndividualParams
    - plotCostFunc

"""


####################################################################

def analyzeSensitivity():
    def analyzeOneThing(path,burnin=0):
        with open(path,'r') as f: raw_output = f.read().split('\n')[1:]
        output = np.array(map(lambda x: np.array(ast.literal_eval(x)),raw_output))
        if  burnin: output=output[burnin:,:]
        ones = np.ones((42,))
        s    = np.std(output,axis=0)
        ns  = np.divide(ones, s, out=np.zeros_like(s), where=s!=0)
        return output,ns

    with open('12.9.17_beta1000_corrected_costs.txt','r') as f:
        raw_costs = np.array(f.read().split('\n')[1:])[2000:] #ASSUME 2000 burnin
        #no burnin
#       raw_costs = np.array(f.read().split('\n')[1:])
    base_output,base_ns  = analyzeOneThing('12.9.17_beta1000_corrected_inputs.txt', 2000)
    f_output,f_ns        = analyzeOneThing('12.9.17_beta1000_forces.txt',2000)
    e_output,e_ns        = analyzeOneThing('12.9.17_beta1000_energies.txt',2000)

    return {'result' : base_ns,'cost': raw_costs,'trajs': base_output
                ,'result_forceonly':f_ns,'result_energyonly':e_ns
                ,'trajs_forceonly':f_ns,'trajs_energyonly':e_ns}

param_labels = ['$C H \phi_{0}$', '$C H A_{1}$', '$C H A_{2}$', '$C H A_{3}$'
    , '$C H A_{4}$', '$C C \phi_{0}$', '$C C A_{1}$', '$C C A_{2}$', '$C C A_{3}$'
    , '$C C A_{4}$', '$H H \phi_{0}$', '$H H A_{1}$', '$H H A_{2}$', '$H H A_{3}$'
    , '$H H A_{4}$', '$C C sss h_{0}$', '$C C sss B_{1}$', '$C C sss B_{2}$'
    , '$C C sps h_{0}$', '$C C sps B_{1}$', '$C C sps B_{2}$', '$C C pps h_{0}$'
    , '$C C pps B_{1}$', '$C C pps B_{2}$', '$C C ppp h_{0}$', '$C C ppp B_{1}$'
    , '$C C ppp B_{2}$', '$H C sss h_{0}$', '$H C sss B_{1}$', '$H C sss B_{2}$'
    , '$H C sps h_{0}$', '$H C sps B_{1}$', '$H C sps B_{2}$', '$H H sss h_{0}$'
    , '$H H sss B_{1}$', '$H H sss B_{2}$', '$H \epsilon_{s}$', '$H \epsilon_{p}$'
    , '$H U$', '$C \epsilon_{s}$', '$C \epsilon_{p}$', '$C U$']

noH_param_labels = ['$C H \phi_{0}$', '$C H A_{1}$', '$C H A_{2}$', '$C H A_{3}$'
    , '$C H A_{4}$', '$C C \phi_{0}$', '$C C A_{1}$', '$C C A_{2}$', '$C C A_{3}$'
    , '$C C A_{4}$', '$C C sss h_{0}$', '$C C sss B_{1}$', '$C C sss B_{2}$'
    , '$C C sps h_{0}$', '$C C sps B_{1}$', '$C C sps B_{2}$', '$C C pps h_{0}$'
    , '$C C pps B_{1}$', '$C C pps B_{2}$', '$C C ppp h_{0}$', '$C C ppp B_{1}$'
    , '$C C ppp B_{2}$', '$H C sss h_{0}$', '$H C sss B_{1}$', '$H C sss B_{2}$'
    , '$H C sps h_{0}$', '$H C sps B_{1}$', '$H C sps B_{2}$','$H \epsilon_{s}$'
    , '$H \epsilon_{p}$', '$H U$', '$C \epsilon_{s}$', '$C \epsilon_{p}$', '$C U$']

def noHanalyzeSensitivity():
    def noHanalyzeOneThing(path,burnin=0):
        with open(path,'r') as f: raw_output = f.read().split('\n')[1:]
        output = np.array(map(lambda x: np.array(ast.literal_eval(x)),raw_output))
        noH_output = np.delete(output, [10, 11, 12, 13, 14, 33, 34, 35], axis = 1)
        if burnin: noH_output = noH_output[burnin:,:]
        ones = np.ones((34,))
        s = np.std(noH_output, axis=0)
        ns = np.divide(ones, s, out=np.zeros_like(s), where=s!=0)
        return output, ns
    
    base_output, base_ns = noHanalyzeOneThing('12.9.17_beta1000_corrected_inputs.txt',2000)
    return {'result' : base_ns, 'trajs': base_output}

def heat():
    from scipy import stats
    f,ax = plt.subplots(nrows=1,ncols=1)

    trajs   = analyzeSensitivity()['trajs']
    a = np.zeros((42, 42))
    for i in range(42):
        for j in range(i,42):
            a[i,j] = stats.linregress(trajs[:,i],trajs[:,j]).rvalue
            a[j,i] = stats.stats.pearsonr(trajs[:,i],trajs[:,j])[0]
    store = np.sort(a, axis=None)
    print store[1721]
    print store[1720]
    plt.imshow(a, cmap='hot', interpolation='nearest')
    plt.colorbar()
    ax.set_xticks(np.arange(42)); ax.set_xticklabels(param_labels)
    ax.set_yticks(np.arange(42)); ax.set_yticklabels(param_labels)
    plt.xticks(fontsize=4,rotation=45)
    plt.yticks(fontsize=4,rotation=45)
    ax.xaxis.label.set_size(4)
    ax.set_xlabel('Parameter', fontsize = 10);ax.set_ylabel('Parameter',fontsize = 10)

#   plt.show()

def trajectory(i):
    """
    i < 1: plot cost trajectory
    otherwise: plot trajectory of parameter i (starting at 1)
    """
    f,ax = plt.subplots(nrows=1,ncols=1)

    if i<1: traj= analyzeSensitivity()['cost']
    else:   traj = analyzeSensitivity()['trajs'][:,i-1]
    ax.plot(np.arange(len(traj)),traj)
#   ax.set_xlabel('Time');ax.set_ylabel(param_labels[i])
    ax.set_xlabel('Time', fontsize = 20); ax.set_ylabel('Cost function (energy + forces)', fontsize = 20)
    plt.show()

def noHplotBaseCaseIndividualParams():
    result = noHanalyzeSensitivity()['result']
    f, ax = plt.subplots(nrows=1,ncols=1)
    ax.set_xticks(np.arange(1,35))
    ax.set_xticklabels(noH_param_labels)
    for tick in ax.get_xticklabels(): tick.set_rotation(45)
    for i in range(34):
        ax.bar(i+1,result[i],width=0.5,color = 'red',label='_nolegend_')
    plt.show()

def  plotBaseCaseIndividualParams():
    result   = analyzeSensitivity()['result']
    f,ax = plt.subplots(nrows=1,ncols=1)
    ax.set_xticks(np.arange(1,43))
    ax.set_xticklabels(param_labels)
    for tick in ax.get_xticklabels(): tick.set_rotation(45)
    for i in range(42):
        ax.bar(i+1,result[i],width=0.5,color = 'red',label='_nolegend_')
    ax.set_xlabel('Parameter');ax.set_ylabel(r'$\frac{1}{\sigma}$', fontsize = 20)
    plt.show()


def plotCostFunc(normalize=False):
    results = [ analyzeSensitivity()['result']
                 ,analyzeSensitivity()['result_forceonly']
                 ,analyzeSensitivity()['result_energyonly']]

    if normalize: map(lambda x: x/np.sum(x),results)
    f,ax = plt.subplots(nrows=1,ncols=1)
    ax.set_xticks(np.arange(1,43))
    ax.set_xticklabels(param_labels)
    for tick in ax.get_xticklabels(): tick.set_rotation(45)
    colors = ['red','blue','green']
    leg   = ['Energy and Forces','Forces','Energies']
    shift  = [.666,1,1.333]
    for i,r in enumerate(results):
        for j in range(42):
            ax.bar(j+shift[i],r[j],width=0.3,color = colors[i],label=leg[i] if j==0 else None)
    ax.legend()
    ax.set_xlabel('Parameter', fontsize = 15);ax.set_ylabel(r'$\frac{1}{\sigma}$', fontsize = 20)
    ax.set_title('Quantification of cost function', fontsize = 20)
    plt.show()


def PCA_plot(plot_variance=False):
    f,ax = plt.subplots(nrows=1,ncols=1)

    trajs   = analyzeSensitivity()['trajs']
    pca       = PCA(n_components=42,whiten=True)
    pca.fit(trajs)
    if plot_variance:
        x         = np.arange(1,43,1)
        y         =  pca.explained_variance_ratio_*100
        ax.semilogy(x[:-1],y[:-1],linestyle='None', marker='o')
        ax.set_title('PCA Analysis of TB Parameter Space',fontsize=18)
        ax.set_xlabel('Number of principal components',fontsize=18)
        ax.set_ylabel('Fraction of variance explained',fontsize=18)
        plt.show()
    else:
        pc1 = sorted(zip(map(abs,pca.components_[0]),param_labels))
        for i in range(42):
            ax.bar(i-.3,abs(pca.components_[0][i]),width=0.2,color = 'red',label='1st' if i==0 else None)
            ax.bar(i,abs(pca.components_[1][i]),width=0.2,color = 'blue',label='2nd' if i==0 else None)
            ax.bar(i+.3,abs(pca.components_[2][i]),width=0.2,color = 'green',label='3rd' if i==0 else None)
        ax.set_xticks(np.arange(42))
        ax.set_xticklabels(param_labels)       
        for tick in ax.get_xticklabels(): tick.set_rotation(45)
        ax.legend()
        ax.set_xlabel('Parameter',fontsize=15);ax.set_ylabel('Magnitude of Principal Component Loading',fontsize = 12)
        plt.show()

def noH_PCA_plot(plot_variance=False):
    f,ax = plt.subplots(nrows=1,ncols=1)
    trajs  = noHanalyzeSensitivity()['trajs']
    pca    = PCA(n_components=34, whiten=True)
    pca.fit(trajs)
    if plot_variance:
        x       = np.arange(1, 35, 1)
        y       = pca.explained_variance_ratio_*100
        ax.semilogy(x[:-1],y[:-1],linestyle='None', marker='o')
        ax.set_xlabel('Number of principal components',fontsize=18)
        ax.set_ylabel('Fraction of variance explained',fontsize=18)
    else:
        pc1 = sorted(zip(map(abs,pca.components_[0]),noH_param_labels))
        for i in range(34):
            ax.bar(i-.3, abs(pca.components_[0][i]),width=0.2,color = 'red',label='1st' if i==0 else None)
            ax.bar(i,abs(pca.components_[1][i]),width=0.2,color = 'blue',label='2nd' if i==0 else None)
            ax.bar(i+.3,abs(pca.components_[2][i]),width=0.2,color = 'green',label='3rd' if i==0 else None)
    ax.set_xticks(np.arange(34))
    ax.set_xticklabels(noH_param_labels)
    for tick in ax.get_xticklabels(): tick.set_rotation(45)
    ax.legend()
    ax.set_xlabel('Parameter');ax.set_ylabel('Magnitude of Principal Component Loading')
    plt.show()


if __name__ == '__main__':
    PCA_plot()
#   trajectory(0)
#   noH_PCA_plot()
#   noHplotBaseCaseIndividualParams()
#   analyzeSensitivity()
#   heat()
#   plotCostFunc()
    #plotCostFunc()#plotBaseCaseIndividualParams()#heat()#plotBaseCaseIndividualParams()# trajectory(0) # #analyzeSensitivity()#
