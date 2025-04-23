from read import *
src='../output/'

ncase,ndata,nspec,ntrg,spec,trgSpec=readSpeciesAndTargets(src)

print('Species:',nspec)
for item in spec: print(item)
print('Target species:',ntrg)
for item in trgSpec: print(item)

print('OICs:')
for i in range(1,ncase+1):
    for j in range(1,ndata+1):
        oic=read(src,i,j,nspec,ntrg)
        for itrg in range(ntrg):
            print('Target:', trgSpec[itrg])
            for jspec in range(nspec):
                print(spec[jspec],' ',oic[jspec][itrg])




