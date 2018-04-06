import pyETS
import libETS.data_model as dm
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import patches as mpatches
from matplotlib import collections


class Limits(object):
    def __init__(self):
        self._xmin = 1e20
        self._xmax = -1e20
        self._ymin = 1e20
        self._ymax = -1e20

    def include (self,x,y,r=0.):
        self._xmin = min(self._xmin,x-r)
        self._xmax = max(self._xmax,x+r)
        self._ymin = min(self._ymin,y-r)
        self._ymax = max(self._ymax,y+r)

    @property
    def xmin(self):
        return self._xmin

    @property
    def xmax(self):
        return self._xmax

    @property
    def ymin(self):
        return self._ymin

    @property
    def ymax(self):
        return self._ymax

def plot_cobra(c, patches, limits, facecolor='none', edgecolor='black', plot_dot=False):
    x, y, r = c.center.real, c.center.imag, c.innerLinkLength + c.outerLinkLength
    limits.include(x, y, r)
    circle = mpatches.Circle((x, y), r, facecolor=facecolor, edgecolor=edgecolor, lw=1.)
    patches.append(circle)
    if plot_dot:
        x, y, r = c.dotcenter.real, c.dotcenter.imag, c.rdot
        circle = mpatches.Circle((x, y), r, facecolor=facecolor, edgecolor=edgecolor, lw=1.)
        patches.append(circle)

def plot(cbr):
    def plotNodes(ax, nodes, limits):
        # targets
        patches = []
        for i, c in enumerate(nodes):
            plot_cobra(c, patches, limits, plot_dot=True)

        collection = collections.PatchCollection(patches, match_original=True)
        ax.add_collection(collection)

    limits=Limits()
    fig = plt.figure(figsize=[15, 15])
    ax = plt.subplot(111)  # note we must use plt.subplots, not plt.subplot
    plotNodes(ax, cbr, limits)
    ax.set_ylim([limits.ymin, limits.ymax])
    ax.set_xlim([limits.xmin, limits.xmax])
    plt.show()

def plot_assignment(cobras,targets,tpos,res):
    xmin, xmax, ymin, ymax = 1e3,-1e3,1e3,-1e3
    fig = plt.figure(figsize=[15, 15])
    ax = plt.subplot(111)  # note we must use plt.subplots, not plt.subplot
    patches = []
    limits=Limits()
    for tidx, cidx in res.items():
        c = cobras[cidx]
        tp = tpos[tidx]
        t = targets[tidx]
        color = "red" if isinstance(t, dm.ScienceTarget) else "green"
        plot_cobra(c, patches, limits, edgecolor=color)
        x, y = c.center.real, c.center.imag
        tx, ty = tp.real, tp.imag
        line = mlines.Line2D([x, tx], [y, ty], lw=1, color= color)
        ax.add_line(line)
    collection = collections.PatchCollection(patches, match_original=True)
    ax.add_collection(collection)
    ax.set_ylim([limits.ymin, limits.ymax])
    ax.set_xlim([limits.xmin, limits.xmax])
    plt.show()

# define locations of the input files
catalog_path = "/home/martin/codes/ets_fiber_assigner/pfs_target_list"
fscience_targets = catalog_path+"/pfs_preliminary_target_cosmology.dat"
#fscience_targets = catalog_path+"/pfs_preliminary_target_galaxy.dat"
fcal_stars       = catalog_path+"/pfs_preliminary_target_cosmology_fcstars.dat"
fsky_pos         = catalog_path+"/pfs_preliminary_target_cosmology_sky.dat"

# read all targets into a sigle list, giving them their proper types
tgt = dm.readScientificFromFile(fscience_targets, "sci")
tgt += dm.readCalibrationFromFile(fcal_stars, "cal")
tgt += dm.readCalibrationFromFile(fsky_pos, "sky")

# get a complete, idealized focal plane configuration
cobras = dm.getFullFocalPlane()

# build reduced Cobra list to speed up calculation
cobras = [c for c in cobras if abs(c.center) < 30]


#plot(cobras)
# point the telescope at the center of all science targets
raTel, decTel = dm.telescopeRaDecFromFile(fscience_targets)
posang = 0.
otime = "2016-04-03T08:00:00Z"
telescopes=[]
nvisit = 20
for _ in range(nvisit):
    telescopes.append(dm.Telescope(cobras, 1., raTel+np.random.normal()*1e-2, decTel+np.random.normal()*1e-2, posang, otime))
    print (telescopes[-1]._ra, telescopes[-1]._dec)

tpos = [tele.get_fp_positions(tgt) for tele in telescopes]

classdict = {}
classdict["sci_P1"] = {"nonObservationCost": 100, "partialObservationCost": 1e9}
classdict["sci_P2"] = {"nonObservationCost": 90, "partialObservationCost": 1e9}
classdict["sci_P3"] = {"nonObservationCost": 80, "partialObservationCost": 1e9}
classdict["sci_P4"] = {"nonObservationCost": 70, "partialObservationCost": 1e9}
classdict["sci_P5"] = {"nonObservationCost": 60, "partialObservationCost": 1e9}
classdict["sci_P6"] = {"nonObservationCost": 50, "partialObservationCost": 1e9}
classdict["sci_P7"] = {"nonObservationCost": 40, "partialObservationCost": 1e9}
classdict["sky"] = {"numRequired": 2}
classdict["cal"] = {"numRequired": 1}

res = dm.observeWithNetflow(telescopes[0].Cobras, tgt, tpos, classdict, 300.)
print (len(res))
for vis, tp in zip(res,tpos):
    plot_assignment(cobras, tgt, tp, vis)

