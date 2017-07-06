import sys, os, math
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.gridspec import GridSpec
from core.TimeseriesCore import TimeseriesCore as tsc
from sklearn import cluster
from itertools import cycle

data_dir = '/Users/anatale/school/UCSF/Grabe_lab/data/traj_features/50ps/state_metrics'

data_names = [
'traakWT_full_npt',
'traakWT_s1s3_npt',
'traakWT_s2s4_npt',
'traakG124I_full_npt',
'traakG124I_s1s3_npt',
'traakG124I_s2s4_npt',
'traakTM4_npt'
]

all_datasets = {}
reader = tsc(verbose=False)
for key in data_names:
    all_datasets[key] = reader._data_reader(os.path.join(data_dir, '%s_state_metrics_all_frames.dat' % key))

all_datadicts = {}
for key in data_names:
    all_datadicts[key] = all_datasets[key].dataset_to_dict()

total_frames = 0
idxs = {}
for key in data_names:
    start_idx = total_frames
    total_frames += all_datasets[key].get_length()
    idxs[key] = (start_idx,total_frames)

for key in data_names:
    print key,idxs[key]

def get_all(subkey):
    target = np.empty((1,0), dtype=float)
    for key in data_names:
        target = np.concatenate([target, all_datadicts[key][subkey]], axis=1)
    return target

wingA = get_all('wingA')
wingB = get_all('wingB')

wingmax = np.amax([np.amax(wingA),np.amax(wingB)])
wingmin = np.amin([np.amin(wingA),np.amin(wingB)])

fenA = get_all('fenA')
fenB = get_all('fenB')

fenmax = np.amax([np.amax(fenA),np.amax(fenB)])
fenmin = np.amin([np.amin(fenA),np.amin(fenB)])

# for key in data_names:
    # f, ax = plt.subplots()
    # ax.scatter(wingA[0,idxs[key][0]:idxs[key][1]],fenA[0,idxs[key][0]:idxs[key][1]], color='b', s=2, edgecolor='none', alpha=0.8)
    # ax.scatter(wingB[0,idxs[key][0]:idxs[key][1]],fenB[0,idxs[key][0]:idxs[key][1]], color='r', s=2, edgecolor='none', alpha=0.8)
    # ax.set_xlim(wingmin,wingmax)
    # ax.set_ylim(fenmin,fenmax)
    # plt.title(key)
    # plt.show()
# for key in data_names:
#     a = wingA[0,idxs[key][0]:idxs[key][1]]
#     b = wingB[0,idxs[key][0]:idxs[key][1]]
#     f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
#     ax1.hist(a, 100, normed=1)
#     ax2.hist(b, 100, normed=1)
#     ax1.set_xlim(wingmin,wingmax)
#     plt.title(key)
#     plt.show()
# for key in data_names:
#     a = fenA[0,idxs[key][0]:idxs[key][1]]
#     b = fenB[0,idxs[key][0]:idxs[key][1]]
#     f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
#     ax1.hist(a, 100, normed=1)
#     ax2.hist(b, 100, normed=1)
#     ax1.set_xlim(fenmin,fenmax)
#     plt.title(key)
#     plt.show()

target = np.concatenate([np.concatenate([wingA,fenA], axis=0),np.concatenate([wingB,fenB], axis=0)], axis=1)
target = target.T
print np.shape(target)

# kmeans clustering
k=6
kmeans = cluster.KMeans(n_clusters=k, n_init=50, random_state=87)
kmeans.fit(target)
labels = kmeans.labels_
centroids = kmeans.cluster_centers_

# f, ax = plt.subplots()
# colors=['violet','limegreen','silver','tomato','seagreen','blue','blueviolet','palevioletred']
# color = cycle(colors).next
# for i in range(k):
#     ds = target[np.where(labels==i)]
#     ax.scatter(ds[:,0],ds[:,1], c=color(), s=1, edgecolor='none')
#     lines = ax.plot(centroids[i,0],centroids[i,1],'kx')
#     plt.setp(lines,ms=15.0)
#     ax.set_xlim(wingmin,wingmax)
#     ax.set_ylim(fenmin,fenmax)
# plt.show()

aaa,bbb = np.split(labels, 2, axis=0)

for key in data_names:
    stdtime = all_datadicts[key]['time'].flatten()
    f, ax = plt.subplots()
    ax.scatter(stdtime, aaa[idxs[key][0]:idxs[key][1]]+0.1, color='red', s=3, edgecolor='none')
    ax.scatter(stdtime, bbb[idxs[key][0]:idxs[key][1]]-0.1, color='blue', s=3, edgecolor='none')
    ax.set_xlim(stdtime[0],stdtime[-1])
    ax.set_ylim(-1,6)
    plt.title(key)
    plt.show()
