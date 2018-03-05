import abc
import numpy as np
import pycconv
import pyETS


class Cobra(object):
    def __init__ (self, ID, center, dotcenter, rdot, li, lo):
        self._ID = str(ID)
        self._center = complex(center)
        self._dotcenter = complex(dotcenter)
        self._rdot = float(rdot)
        self._innerLinkLength = float(li)
        self._outerLinkLength = float(lo)

    @property
    def ID(self):
        """ID of the object : str"""
        return self._ID

    @property
    def center(self):
        """position of Cobra center in the focal plane : pos"""
        return self._center

    @property
    def dotcenter(self):
        """position of dot center in the focal plane : pos"""
        return self._dotcenter

    @property
    def rdot(self):
        """dot radius : float"""
        return self._rdot

    @property
    def innerLinkLength(self):
        """length of the inner link : float"""
        return self._innerLinkLength

    @property
    def outerLinkLength(self):
        """length of the outer link : float"""
        return self._outerLinkLength


class Telescope(object):
    def __init__(self, Cobras, collisionRadius, ra, dec, posang, time):
        self._Cobras = tuple(Cobras)
        self._cobraCollisionRadius = float(collisionRadius)
        self._ra = float(ra)
        self._dec = float(dec)
        self._posang = float(posang)
        self._time = str(time)

    @property
    def Cobras(self):
        """return all Cobras : tuple(Cobra) or list(Cobra)"""
        return self._Cobras

#    @property
#    def Spokes(self):
#        """return geometry of spokes fiducials (TBD)"""

    @property
    def cobraCollisionRadius(self):
        """the radius of a Cobra tip : float
        This is used for collision detection. It could also become a property
        of each individual Cobra, but apparently this quantity is not very
        variable.
        """
        return self._cobraCollisionRadius

    def select_visible_targets(self, tgt):
        for t in tgt:
            t.calc_position(self._ra, self._dec, self._posang, self._time)
        pos = [t.position for t in tgt]
        cbr = []
        for c in self._Cobras:
            cbr.append([c.center, c.innerLinkLength, c.outerLinkLength,
                        c.dotcenter, c.rdot])
        tmp = pyETS.getVis(pos, cbr)
        return [tgt[k] for k in tmp.keys()]

    def subtract_obs_time(self, tgt, obs, tvisit):
        res = []
        for k in obs.keys():
            tgt[k].reduce_time(tvisit)
        for t in tgt:
            if t.obs_time > 0.:
                res.append(t)
        return res

    def observeWithETS(self, tgt, nvisit, tvisit, assigner):
        tgt = self.select_visible_targets(tgt)

        cbr = []
        for c in self._Cobras:
            cbr.append([c.center, c.innerLinkLength, c.outerLinkLength,
                        c.dotcenter, c.rdot])

        res = []
        for i in range(nvisit):
            pos = [t.position for t in tgt]
            pri = [t.priority for t in tgt]
            time = [t.obs_time for t in tgt]
            tmp = pyETS.getObs(pos, time, pri, cbr, assigner)
            d = {}
            for key, value in tmp.items():
                d[tgt[key].ID] = self._Cobras[value].ID
            res.append(d)
            tgt = self.subtract_obs_time(tgt, tmp, tvisit)
            print(len(tgt))
        return res

    def observeWithNetflow(self, tgt, nvisit):
        tgt = self.select_visible_targets(tgt)


class Target(object):
    def __init__(self, ID, ra, dec):
        self._ID = str(ID)
        self._ra = float(ra)
        self._dec = float(dec)
        self._position = None

    @property
    def ID(self):
        """ID of the object : str"""
        return self._ID

    @property
    def ra(self):
        """the rectascension : float"""

    @property
    def dec(self):
        """the declination : float"""

    def calc_position(self, raTel, decTel, posang, time):
        self._position = pycconv.cconv([self._ra], [self._dec],
                                       raTel,decTel,posang,time)[0]

    @property
    def position(self):
        """the position in the focal plane : pos"""
        return self._position


class ScienceTarget(Target):
    def __init__(self, ID, ra, dec, obs_time, pri):
        super(ScienceTarget, self).__init__(ID, ra, dec)
        self._obs_time = float(obs_time)
        self._pri = int(pri)

    @property
    def priority(self):
        """target priority : int"""
        return self._pri

    @property
    def nonObservationCost(self):
        return self._non_obs_ost

    @property
    def partialObservationCost(self):
        return self._partial_obs_cost

    @property
    def obs_time(self):
        """required observation time : float"""
        return self._obs_time

    def reduce_time(self, dt):
        self._obs_time -= float(dt)


class CalibTarget(Target):
    @abc.abstractmethod
    def numRequired():
        """ returns the required number of targets for this target class"""
        pass


class SkyCalibTarget(CalibTarget):
# must implement numRequired()
    pass


class StarCalibTarget(CalibTarget):
# must implement numRequired()
    pass

# ... maybe other calibration classes

def telescopeRaDecFromFile(file):
    with open(file) as f:
        ras=[]
        decs=[]
        ll=f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                ra,dec = (float(tt[1]), float(tt[2]))
                ras.append(ra)
                decs.append(dec)
    return float(np.average(ras)), float(np.average(decs))


def readScientificFromFile(file):
    with open(file) as f:
        res = []
        ll=f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_,ra,dec,tm,pri = (str(tt[0]), float(tt[1]), float(tt[2]),
                                       float(tt[3]), int(tt[4]))
                res.append(ScienceTarget(id_, ra, dec, tm, pri))
    return res


def readCalibrationFromFile(file, cls):
    with open(file) as f:
        res = []
        ll=f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_,ra,dec = (str(tt[0]), float(tt[1]), float(tt[2]))
                res.append(cls(id_, ra, dec))
    return res


def getFocalPlane():
    cobras = pyETS.getAllCobras()
    res = []
    for i in range(len(cobras)):
        ID = "C{}".format(i)
        res.append(Cobra(ID, cobras[i][0], cobras[i][3], cobras[i][4],
                         cobras[i][1], cobras[i][2]))
    return res


catalog_path = "/home/martin/codes/ets_fiber_assigner/pfs_target_list"
fscience_targets = catalog_path+"/pfs_preliminary_target_cosmology.dat"
fcal_stars       = catalog_path+"/pfs_preliminary_target_cosmology_fcstars.dat"
fsky_pos         = catalog_path+"/pfs_preliminary_target_cosmology_sky.dat"

tgt = readScientificFromFile(fscience_targets)
#tgt += readCalibrationFromFile(fcal_stars, StarCalibTarget)
#tgt += readCalibrationFromFile(fsky_pos, SkyCalibTarget)

cobras = getFocalPlane()

# build reduced Cobra list to speed up calculation
cobras = [c for c in cobras if abs(c.center) < 20]

# point the telescope at the center of all science targets
raTel, decTel = telescopeRaDecFromFile(fscience_targets)
posang = 0.
time = "2016-04-03T08:00:00Z"
telescope = Telescope(cobras, 1., raTel, decTel, posang, time)

# select only targets that are visible by the active cobras
res=telescope.observeWithETS(tgt, 20, 300, "draining")
#res=telescope.observeWithNetflow(tgt, 1)
print(res)
