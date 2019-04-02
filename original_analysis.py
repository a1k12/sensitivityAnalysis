import collections,pickle,ast,os
import numpy as np
import matplotlib
#matplotlib.use('GTKAgg')
import matplotlib.pyplot     as plt
####################################################################
####################################################################
####################################################################

def rec_dd(): return collections.defaultdict(rec_dd)

def analyzeSensitivity(folder):

    resultDict = collections.defaultdict(rec_dd)

    for d in os.listdir('.'):
        if 'DS' not in d and 'costs' not in d and 'result' not in d:
            print 'd ',d
            i,j,k,l = map(int,d[:-4].split('_'))
            with open(d,'r') as f:
                raw_output = f.read().split('\n')[1:]

            output = np.array(map(lambda x: np.array(ast.literal_eval(x)),raw_output))
            m = np.abs(np.mean(output,axis=0))
            s = np.std(output,axis=0)
            ns = np.divide(s, m, out=np.zeros_like(s), where=m!=0)

            resultDict[i][j][k][l] = ns

    with open('result.pckl','wb') as f: pickle.dump(resultDict,f)



def    plotBaseCaseIndividualParams():
    with open('result.pckl','rb') as f: result = pickle.load(f)
    f,ax = plt.subplots(nrows=1,ncols=1)

    basecase = result[1][1][2][0]
    for i in range(len(basecase)):
        ax.bar(i,basecase[i],width=0.5,color = 'red',label='_nolegend_')
    plt.show()

def    plotAlphaDependence():
    with open('result.pckl','rb') as f: result = pickle.load(f)
    f,ax = plt.subplots(nrows=1,ncols=1)

    c1 = result[1][1][0][0]
    c2 = result[1][1][1][0]
    c3 = result[1][1][2][0] # base case
    for i in range(len(c1)):
        ax.bar(i-.33,c1[i],width=0.33,color = 'red')
        ax.bar(i,c2[i],width=0.33,color = 'blue')
        ax.bar(i+.33,c3[i],width=0.33,color = 'green')

    legend = ax.legend((r'$\alpha  = 0.05$',r'$\alpha  = 0.01$',r'$\alpha  = 0.005$'),loc='upper right', fancybox=True)
    ax.set_title(r'$\alpha $ dependence',fontsize=18)
    ax.set_xlabel('Parameters',fontsize=18)
    ax.set_ylabel(r'$\frac{\sigma}{\mu}$', fontsize=18)
    plt.show()

def   plotCostFuncDependence():
    with open('result.pckl','rb') as f: result = pickle.load(f)
    f,ax = plt.subplots(nrows=1,ncols=1)

    c1 = result[0][1][2][0]
    c2 = result[1][1][2][0] # base case
    c3 = result[2][1][2][0]
    for i in range(len(c1)):
        ax.bar(i-.33,c1[i],width=0.33,color = 'red',label='1')
        ax.bar(i,c2[i],width=0.33,color = 'blue',label='2')
        ax.bar(i+.33,c3[i],width=0.33,color = 'green',label='3')

    legend = ax.legend(('Forces','Energy+Forces','Energy'),loc='upper right', fancybox=True)
    legend.draggable()
    ax.set_title(r'Cost function dependence',fontsize=18)
    ax.set_xlabel('Parameters',fontsize=18)
    ax.set_ylabel(r'$\frac{\sigma}{\mu}$',fontsize=18)

    plt.show()

def  plotBetaDependence():
    with open('result.pckl','rb') as f: result = pickle.load(f)
    f,ax = plt.subplots(nrows=1,ncols=1)

    c1 = result[1][1][2][0] #base case
    c2 = result[1][1][2][1]
    for i in range(len(c1)):
        ax.bar(i-.25,c1[i],width=0.5,color = 'red',label='1')
        ax.bar(i+.25,c2[i],width=0.5,color = 'blue',label='2')

    legend = ax.legend((r'$\beta = 1$',r'$\beta = 10$'),loc='upper right',fancybox=True)
    legend.draggable()
    ax.set_title(r'$\beta $ dependence',fontsize=18)
    ax.set_xlabel('Parameters',fontsize=18)
    ax.set_ylabel(r'$\frac{\sigma}{\mu}$',fontsize=18)

    plt.show()


def    plotJumpFuncDependence():
    with open('result.pckl','rb') as f: result = pickle.load(f)
    f,ax = plt.subplots(nrows=1,ncols=1)

    j1 = result[1][1][2][0] # base case
    j2 = result[1][0][2][0]

    for i in range(len(j1)):
        ax.bar(i-.25,j1[i],width=0.5,color = 'red',label='_nolegend_')
        ax.bar(i+.25,j2[i],width=0.5,color = 'blue',label='_nolegend_')
    plt.show()
def PCA_plot():
    f,ax = plt.subplots(nrows=1,ncols=1)

    from sklearn.decomposition import PCA
    with open('1_1_2_0.txt','r') as f: raw_output = f.read().split('\n')[1:]
    output = np.array(map(lambda x: np.array(ast.literal_eval(x)),raw_output))
    print '\nraw output ',output
    s         = (1101,1)
    ones      = np.ones(s)
    raw_mean  = np.mean(output,axis=0)
    raw_sd    = np.std(output,axis=0)
    m         = np.dot(ones,[raw_mean])
    sd        = np.dot(ones,[raw_sd])

    output    = output - m
    pca_input = np.nan_to_num(np.divide(output, sd))

    pca       = PCA(n_components=42)#,whiten=True)
    pca.fit(pca_input)
    x         = np.arange(1,43,1)
    y         =  pca.explained_variance_/np.sum(pca.explained_variance_)*100

    ax.semilogy(x[:-1],y[:-1],linestyle='None', marker='o')
    ax.set_title('PCA Analysis of TB Parameter Space',fontsize=18)
    ax.set_xlabel('Number of principal components',fontsize=18)
    ax.set_ylabel('Fraction of variance explained',fontsize=18)
    plt.show()
if __name__ == '__main__': PCA_plot()
