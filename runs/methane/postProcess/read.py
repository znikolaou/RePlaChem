from scipy.io import FortranFile
import numpy as np
intVar='i4'
dblVar='f8'

def readSpeciesAndTargets(src):
    fl=src+'/speciesAndTargets.txt'
    f=open(fl,'r')
    f.readline()
    ncase=int(f.readline())
    f.readline()
    ndata=int(f.readline())
    f.readline()
    nspec=int(f.readline())
    f.readline()
    spec=[]
    for i in range(nspec):
        spec.append(f.readline())
    f.readline()
    ntrg=int(f.readline())
    f.readline()
    trgSpec=[]
    for i in range(ntrg):
        trgSpec.append(f.readline())
    f.close()
    
    return ncase,ndata,nspec,ntrg,spec,trgSpec

def read(src,i,j,nspec,ntrg):
    fl=src+'/case'+str(i).zfill(3)+'data'+str(j).zfill(8)+'.dat'
    f=FortranFile(fl,'r')
    data=f.read_record(intVar, intVar, \
            np.dtype((dblVar,(nspec,ntrg))) )

    ntrg=data[0]
    nspec=data[1]
    oic=data[2]
    f.close()
    
    return oic



