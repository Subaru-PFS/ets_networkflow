import abc
import numpy as np
import pycconv
import pyETS
from collections import OrderedDict, defaultdict
import time
import pulp

def build_network(cobras, targets, nvisits, tvisit):
    Cv_i = defaultdict(list)  # Cobra visit inflows
    Tv_o = defaultdict(list)  # Target visit outflows
    Tv_i = defaultdict(list)  # Target visit inflows
    T_o = defaultdict(list)  # Target outflows (only science targets)
    T_i = defaultdict(list)  # Target inflows (only science targets)
    CTCv_o = defaultdict(list)  # Calibration Target visit outflows
    STC_o = defaultdict(list)  # Science Target outflows
    prob = pulp.LpProblem("problem", pulp.LpMinimize)
    cost = pulp.LpVariable("cost", 0)
    #cost = 0.

    nreqvisit = []
    for t in targets:
        if isinstance(t, ScienceTarget):
            nreqvisit.append(int(t.obs_time/tvisit))
        else:
            nreqvisit.append(0)

    # determine visibilities
    pos = [t.position for t in targets]
    cbr = [[c.center, c.innerLinkLength, c.outerLinkLength,
                    c.dotcenter, c.rdot] for c in cobras]

    vis = pyETS.getVis(pos, cbr)

    def newvar(lo, hi):
        newvar._varcount+=1
        print("var",newvar._varcount)
        return pulp.LpVariable("v{}".format(newvar._varcount), lo, hi, cat=pulp.LpInteger)
    newvar._varcount = 0

    vis = [vis for _ in range(nvisits)]  # just replicate for now
    for ivis in range(nvisits):
        for tidx, val in vis[ivis].items():
            tgt = targets[tidx]
            if isinstance(tgt, ScienceTarget):
                f = newvar(0,1)#pulp.LpVariable("T{}_Tv{}_v{}".format(tidx, tidx, ivis), 0, 1, cat=pulp.LpInteger)
                T_o[tidx].append(f)
                Tv_i[(tidx,ivis)].append(f)
                if len(T_o[tidx])==1:  # freshly created
                    f = newvar(0,1)#pulp.LpVariable("STC{}_T{}".format(tidx, tidx), 0, 1, cat=pulp.LpInteger)
                    T_i[tidx].append(f)
                    STC_o[type(tgt)].append(f)
                    if len(STC_o[type(tgt)])==1:  # freshly created
                        f = newvar(0,None)#pulp.LpVariable("STC{}_SINK".format(tidx), 0, None, cat=pulp.LpInteger)
                        STC_o[type(tgt)].append(f)
                        cost += f*tgt.nonObservationCost
                    f = newvar(0,None)#pulp.LpVariable("T{}_SINK".format(tidx), 0, None, cat=pulp.LpInteger)
                    T_o[tidx].append(f)
                    cost += f*tgt.partialObservationCost
            elif isinstance(tgt, CalibTarget):
                f = newvar(0,1)#pulp.LpVariable("CTCv{}_Tv{}_v{}".format(tidx, tidx, ivis), 0, 1, cat=pulp.LpInteger)
                Tv_i[(tidx,ivis)].append(f)
                CTCv_o[(type(tgt),ivis)].append(f)
            for cidx in val:
                f = newvar(0,1)#pulp.LpVariable("Tv{}_Cv{}_v{}".format(tidx, cidx, ivis), 0, 1, cat=pulp.LpInteger)
                Cv_i[(cidx,ivis)].append(f)
                Tv_o[(tidx,ivis)].append((f, cidx))

    prob += cost
    # every Cobra can observe at most one target per visit
    for inflow in Cv_i.values():
        prob += pulp.lpSum([f for f in inflow]) <= 1

    # every calibration target class must be observed a minimum number of times
    # every visit
    for key, value in CTCv_o.items():
        prob += pulp.lpSum([v for v in value]) >= key[0].numRequired()

    # inflow and outflow at every Tv node must be balanced
    for key, ival in Tv_i.items():
        oval = Tv_o[key]
        prob += pulp.lpSum([v for v in ival]+[-v[0] for v in oval]) == 0

    # inflow and outflow at every T node must be balanced
    for key, ival in T_i.items():
        oval = T_o[key]
        nvisits = nreqvisit[key]
        prob += pulp.lpSum([nvisits*v for v in ival]+[-v for v in oval]) == 0

    # Science targets must be either observed or go to the sink
    for key, val in STC_o.items():
        prob += pulp.lpSum([v for v in val]) == 1000000

    #print(prob)
    status = prob.solve(pulp.COIN_CMD(msg=1, keepFiles=0,
                                      maxSeconds=1000,
                                          threads=2, dual=10.))
    for k1, v1 in Tv_o.items():
        for i2 in v1:
            visited = pulp.value(i2[0]) > 0
            if visited:
                tidx, ivis = k1
                cidx = i2[1]
                print("visit:{} tgt {} observed by cobra {}".format(ivis, tidx, cidx))
    #print(prob)
    exit()

class Cobra(object):
    """An object holding all relevant information describing a single Cobra
    positioner.
    This includes center position, black dot position, black dot radius, and
    inner and outer link lengths.
    All lengths in mm, positions are stored as complex numbers.
    """
    def __init__(self, ID, center, dotcenter, rdot, li, lo):
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
    """An object describing a telescope configuration to be used for observing
    a target field. Includes a list of Cobras, telescope RA/Dec and position
    angle and an observation time according to ISO8601 UTC.
    The minimum allowed distance between two Cobra tips is also stored here.
    Based on this information, and given a list of targets, this object can
    compute assignment strategies using different algorithms (currently network
    flow and ETs approaches).
    """
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
        return res

    def observeWithNetflow(self, tgt, nvisit, tvisit):
        for t in tgt:
            t.calc_position(self._ra, self._dec, self._posang, self._time)
        #tgt = self.select_visible_targets(tgt)
        build_network (self._Cobras, tgt, nvisit, tvisit)


class Target(object):
    """Base class for all target types observable with PFS. All targets are
    initialized with RA/Dec and an ID string. From RA/Dec, the target can
    determine its position on the focal plane, once also the telescope attitude
    is known. """
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
                                       raTel, decTel, posang, time)[0]

    @property
    def position(self):
        """the position in the focal plane : pos"""
        return self._position


class ScienceTarget(Target):
    """Derived from the Target class, with the additional attributes priority
    and observation time.

    All different types of ScienceTarget need to be derived from this class."
    """
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
        if self._pri == 1:
            return 1000.
        elif self._pri == 2:
            return 100.
        elif self._pri == 3:
            return 10.
        else:
            return 1.

    @property
    def partialObservationCost(self):
        return 1e9

    @property
    def obs_time(self):
        """required observation time : float"""
        return self._obs_time

    def reduce_time(self, dt):
        self._obs_time -= float(dt)


class CalibTarget(Target):
    """Derived from the Target class.
    Needs to define a class variable called 'nonObservationCost'
    and a property called 'classname'"""
    @abc.abstractmethod
    def numRequired():
        """returns the required number of targets for this target class"""
        pass


def telescopeRaDecFromFile(file):
    with open(file) as f:
        ras = []
        decs = []
        ll = f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                ra, dec = (float(tt[1]), float(tt[2]))
                ras.append(ra)
                decs.append(dec)
    return float(np.average(ras)), float(np.average(decs))


def readScientificFromFile(file):
    with open(file) as f:
        res = []
        ll = f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_, ra, dec, tm, pri = (
                    str(tt[0]), float(tt[1]), float(tt[2]),
                    float(tt[3]), int(tt[4]))
                res.append(ScienceTarget(id_, ra, dec, tm, pri))
    return res


def readCalibrationFromFile(file, cls):
    with open(file) as f:
        res = []
        ll = f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_, ra, dec = (str(tt[0]), float(tt[1]), float(tt[2]))
                res.append(cls(id_, ra, dec))
    return res


def getFullFocalPlane():
    cobras = pyETS.getAllCobras()
    res = []
    for i in range(len(cobras)):
        ID = "{}".format(i)
        res.append(Cobra(ID, cobras[i][0], cobras[i][3], cobras[i][4],
                         cobras[i][1], cobras[i][2]))
    return res
