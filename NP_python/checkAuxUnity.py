import matplotlib.pyplot as plt

import numpy as np

import csv

import os

from glob import iglob

from os.path import getctime

import sys



def checkAux(path):

    if os.path.isdir(path):

        fileNamePattern=path+'/**/*position.txt'

        LatestFile = max(iglob(fileNamePattern,recursive=True),key=getctime)
        print(LatestFile)
    else:

        LatestFile=path

    f=open(LatestFile)



    #reader = csv.reader(f,delimiter=',')

    #x=list(reader)

    #header = f.readline()

    all_Y = np.genfromtxt(f,delimiter="\t", dtype='<f8',skip_header=0)



    timediff = np.diff(all_Y[:,0])

    bin_edges = np.linspace(start=10, stop=30,num=20 + 1, endpoint=True)

    n, bins, patches=plt.hist(timediff,bins=bin_edges/1000,alpha=0.7, rwidth=0.85)   

    m=n.max()/2

    plt.text(0.02, m, r'$\mu={:.3f}$'.format(np.mean(timediff)))



    plt.grid(axis='y', alpha=0.75)

    plt.xlabel('[s]')

    plt.ylabel('Frequency')

    plt.title('Frametimes, n: '+str(all_Y.shape[0]))

    



    #headers=header.split('\t')

    #headers = headers[1:]

    nlines= all_Y.shape[1]-1 #don't want to plot time

    fig, ax = plt.subplots()

    step=1.5

    for i in range(nlines):

        dat=[sub_x[i+1] for sub_x in all_Y]

        dat=np.array(dat,dtype='float')

        dat = (dat - dat.min())/np.ptp(dat)

        ax.plot(dat+step*i)



    #labels = [i for i in headers]

    ax.set_yticks(np.arange(0,nlines)*step+.25)

    #ax.set_yticklabels(labels)

    for t, l in zip(ax.get_yticklabels(), ax.lines):

        t.set_color(l.get_color())

    plt.title(LatestFile.split('\\')[-1])

    plt.show()









if __name__=='__main__':

    if len(sys.argv)==1:

        path = r'D:/temp'

    else:

        path = sys.argv[1]

    checkAux(path)