import sys, os, math
import numpy as np
import ExtendedTSC
import traak_plots
from sklearn import cluster
import matplotlib.pyplot as plt

#np.set_printoptions(threshold=np.nan)

# stepsize in ps
stepsize = 500
save = False
outdirname = None

# generate list of input dat files:
indir = os.path.abspath('../measurements/20170508_basic')
infiles = []
for infile in os.listdir(indir):
    if infile.endswith('.dat'):
        infiles.append(infile)
infiles = [os.path.join(indir, i) for i in infiles]

datalist = []
for filename in infiles:
    datalist.append(ExtendedTSC.ExtendedTSC(filename))
for elem in datalist:
    elem.simplify_indexing()

if save:
    plotdir = os.path.join(indir,outdirname)
    os.mkdir(plotdir)
    os.chdir(plotdir)

# prep clustering on all chi1/chi2 points by building one big array
idxref = np.empty((len(datalist),2), dtype=int)
for idx,elem in enumerate(datalist):
    if idx == 0:
        idxref[idx,0] = np.shape(elem.access['W262_chi1A'])[1]
        idxref[idx,1] = np.shape(elem.access['W262_chi1B'])[1] + idxref[idx,0]
    else:
        idxref[idx,0] = np.shape(elem.access['W262_chi1A'])[1] + idxref[idx-1,1]
        idxref[idx,1] = np.shape(elem.access['W262_chi1B'])[1] + idxref[idx,0]
print idxref[-1,1]

alldihedral = np.empty((idxref[-1,1],2), dtype=float)
print np.shape(alldihedral)

for idx,elem in enumerate(datalist):
    if idx == 0:
        alldihedral[0:idxref[idx,0],0] = np.copy(elem.access['W262_chi1A'])
        alldihedral[0:idxref[idx,0],1] = np.copy(elem.access['W262_chi2A'])
        alldihedral[idxref[idx,0]:idxref[idx,1],0] = np.copy(elem.access['W262_chi1B'])
        alldihedral[idxref[idx,0]:idxref[idx,1],1] = np.copy(elem.access['W262_chi2B'])
    else:
        alldihedral[idxref[idx-1,1]:idxref[idx,0],0] = np.copy(elem.access['W262_chi1A'])
        alldihedral[idxref[idx-1,1]:idxref[idx,0],1] = np.copy(elem.access['W262_chi2A'])
        alldihedral[idxref[idx,0]:idxref[idx,1],0] = np.copy(elem.access['W262_chi1B'])
        alldihedral[idxref[idx,0]:idxref[idx,1],1] = np.copy(elem.access['W262_chi2B'])

for idx,n in np.ndenumerate(alldihedral[:,0]):
    if n < -2.0:
        alldihedral[idx,0] += math.pi*2
alldihedral = np.rad2deg(alldihedral)

k=7
kmeans = cluster.KMeans(n_clusters=k, n_init=50)
kmeans.fit(alldihedral)
labels = kmeans.labels_
centroids = kmeans.cluster_centers_

# A
print np.shape(labels)
for idx,elem in enumerate(datalist):
    print elem.trajname
    plotname = elem.trajname.split('/')[-1].split('.')[0]+'_A'
    if idx == 0:
        targetidxA = (0,idxref[idx,0])
    else:
        targetidxA = (idxref[idx-1,1],idxref[idx,0])
    # print targetidxA
    y1 = [
    ('fenA', elem.access['fenestrationA'][0,:]),
    ('expA', elem.access['expansionA'][0,:])
    ]
    y2 = [
    ('cluster', labels[targetidxA[0]:targetidxA[1]])
    ]
    traak_plots.timeseries_dual_yscale(plotname,y1,y2,elem.access['time'][0,:]/1000.0,y2_limits=[0,8])

# B
for idx,elem in enumerate(datalist):
    print elem.trajname
    plotname = elem.trajname.split('/')[-1].split('.')[0]+'_B'
    targetidxB = (idxref[idx,0],idxref[idx,1])
    # print targetidxA
    y1 = [
    ('fenB', elem.access['fenestrationB'][0,:]),
    ('expB', elem.access['expansionB'][0,:])
    ]
    y2 = [
    ('cluster', labels[targetidxB[0]:targetidxB[1]])
    ]
    traak_plots.timeseries_dual_yscale(plotname,y1,y2,elem.access['time'][0,:]/1000.0,y2_limits=[0,8])
    
# print np.shape(labels)
# too_high = alldihedral > 180
# alldihedral[too_high] -= 360
# too_low = alldihedral < -180
# alldihedral[too_low] += 360

# f, ax = plt.subplots()
# for i in range(k):
#     ds = alldihedral[np.where(labels==i)]
#     ax.plot(ds[:,0],ds[:,1],'o')
#     lines = ax.plot(centroids[i,0],centroids[i,1],'kx')
#     plt.setp(lines,ms=15.0)
# ax.set_xlim([-180,180])
# ax.set_ylim([-180,180])
# ax.set_xlabel('Chi1')
# ax.set_ylabel('Chi2')
# plt.show()

# for elem in datalist:
#     print elem.trajname
#     plotname = elem.trajname.split('/')[-1].split('.')[0]
#     data1 = [
#     ('fenA', elem.access['fenestrationA'][0,:]),
#     ('fenB', elem.access['fenestrationB'][0,:])
#     ]
#     data2 = [
#     ('expA', elem.access['expansionA'][0,:]),
#     ('expB', elem.access['expansionB'][0,:])
#     ]
#     traak_plots.timeseries_dual_yscale(plotname, data1, data2, elem.access['time'][0,:]/1000.0)
