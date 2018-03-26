import abc
import numpy as np
import pycconv
import pyETS
from collections import OrderedDict
from pfs_netflow.survey_plan import buildSurveyPlan
from pfs_netflow.lp import buildLPProblem, computeStats, solve
import time
import pulp
from pfs_netflow.plotting import plotSurveyPlan, plotFocalPlane
import pfs_netflow.datamodel as dm
from collections import defaultdict

def build_network(cobras, targets, nvisits, tvisit):
    Cv_i = defaultdict(list)  # Cobra visit inflows
    Tv_o = defaultdict(list)  # Target visit outflows
    Tv_i = defaultdict(list)  # Target visit outflows
    T_o = defaultdict(list)  # Target visit outflows
    T_i = defaultdict(list)  # Target visit outflows
    CTCv_o = defaultdict(list)  # Target visit outflows
    STC_o = defaultdict(list)  # Target visit outflows
    prob = pulp.LpProblem("problem", pulp.LpMinimize)

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

    vis = [vis for _ in range(nvisits)]  # just replicate for now
    for tidx, val in vis[0].items():
        tgt = targets[tidx]
        if isinstance(tgt, ScienceTarget):
            for i in range(nvisits):
                f = pulp.LpVariable("T{}_Tv{}_v{}".format(tidx, tidx, i), 0, 1, cat=pulp.LpInteger)
                T_o[tidx].append(f)
                Tv_i[(tidx,i)].append(f)
            f = pulp.LpVariable("STC{}_T{}".format(tidx, tidx), 0, 1, cat=pulp.LpInteger)
            T_i[tidx].append(f)
            STC_o[type(tgt)].append(f)
        elif isinstance(tgt, CalibTarget):
            for i in range(nvisits):
                f = pulp.LpVariable("CTCv{}_Tv{}_v{}".format(tidx, tidx, i), 0, 1, cat=pulp.LpInteger)
                Tv_i[(tidx,i)].append(f)
                CTCv_o[(type(tgt),i)].append(f)
        for cidx in val:
            for i in range(nvisits):
                f = pulp.LpVariable("Tv{}_Cv{}_v{}".format(tidx, cidx, i), 0, 1, cat=pulp.LpInteger)
                Cv_i[(cidx,i)].append(f)
                Tv_o[(tidx,i)].append(f)

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
        prob += pulp.lpSum([v for v in ival]+[-v for v in oval]) == 0

    # inflow and outflow at every T node must be balanced
    for key, ival in T_i.items():
        oval = T_o[key]
        nvisits = nreqvisit[key]
        prob += pulp.lpSum([nvisits*v for v in ival]+[-v for v in oval]) == 0

    print(prob)

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
        tgt = self.select_visible_targets(tgt)
        build_network (self._Cobras, tgt, nvisit, tvisit)
        # build dictionaries
        # supply_dict:
        # - for all scientific targets: infinity
        # - for all types of calibration targets: take number from class
        # cost_dict:
        # - scientific targets: get from object
        # -
        xcobras = OrderedDict()
        for i, c in enumerate(self._Cobras):
            xcobras[c.ID] = [np.real(c.center), np.imag(c.center)]

        targets = OrderedDict()
        for i, t in enumerate(tgt):
            targets[t.ID] = [np.real(t.position), np.imag(t.position)]

        nreqvisit = []
        for t in tgt:
            if isinstance(t, ScienceTarget):
                nreqvisit.append(int(t.obs_time/tvisit))
            else:
                nreqvisit.append(0)

        # determine visibilities
        pos = [t.position for t in tgt]
        cbr = []
        for c in self._Cobras:
            cbr.append([c.center, c.innerLinkLength, c.outerLinkLength,
                        c.dotcenter, c.rdot])
        vis = pyETS.getVis(pos, cbr)

        visibilities = OrderedDict()
        for key, val in vis.items():
            tid = tgt[key].ID
            cc = [self._Cobras[c].ID for c in val]
            visibilities[tid] = cc

        class_dict = {}
        cost_dict = {}
        supply_dict = {}
        for t in tgt:
            if isinstance(t, ScienceTarget):
                cls = "sci_P{}".format(t.priority)
                class_dict[t.ID] = cls
                supply_dict[cls] = np.inf
                cost_dict[cls] = (t.nonObservationCost,
                                  t.partialObservationCost)
            elif isinstance(t, CalibTarget):
                cls = t.classname()
                class_dict[t.ID] = cls
                supply_dict[cls] = t.numRequired()
                cost_dict[cls] = t.nonObservationCost
            else:
                raise TypeError
        A = 0.
        cost_dict['cobra_move'] = lambda d: d*A

        g = buildSurveyPlan(xcobras, targets, nreqvisit, visibilities,
                            class_dict, cost_dict, supply_dict, nvisit, 500.,
                            [0., 0.])
        print("Building LP problem ...")
        start_time = time.time()
        prob, flows, cost = buildLPProblem(g, cat='Integer')
        # prob, flows, cost = buildLPProblem(g, cat='Continuous')
        time_to_build = time.time() - start_time
        print("Time to build model: {:.4e} s".format(time_to_build))
        # Solve problem!
        print("Solving LP problem ...")
        start_time = time.time()

        status = solve(prob, maxSeconds=100)  # , solver="GUROBI")

        def setflows(g, flows):
            for a in g.arcs.values():
                k = '{}={}'.format(a.startnode.id, a.endnode.id)
                if k in flows:
                    a.flow = pulp.value(flows[k])
        setflows(g, flows)

        res = []
        for i in range(nvisit):
            res.append({})
        for a in g.arcs.values():
            n1, n2 = a.startnode, a.endnode
            if a.flow > 0.01:
                if type(n2) == dm.CobraVisit and type(n1) == dm.TargetVisit:
                    visit = n2.visit
                    cobraID = n2.cobra.id[2:]
                    targetID = n1.target.id[2:]
                    res[visit][targetID] = cobraID
                elif type(n2) == dm.CobraVisit and type(n1) == dm.CalTarget:
                    visit = n2.visit
                    cobraID = n2.cobra.id
                    cobraID = cobraID[2:]
                    vstr = "_v{}".format(visit)
                    targetID = n1.id[2:-len(vstr)]
                    res[visit][targetID] = cobraID

        time_to_solve = time.time() - start_time
        print("Solve status is [{}].".format(pulp.LpStatus[status]))
        print("Time to solve: {:.4e} s".format(time_to_solve))

        stats = computeStats(g, flows, cost)

        NSciTargets = 0
        for t in tgt:
            if isinstance(t, ScienceTarget):
                NSciTargets += 1
        print("{} = {}".format('Value of cost function', pulp.value(stats.cost)))
        print("[{}] out of {} science targets get observed.".format(int(stats.NSciObs), NSciTargets))
        print("For {} out of these all required exposures got allocated.".format(stats.NSciComplete))
        print("{} targets get sent down the overflow arc.".format(stats.Noverflow))
        print("{} out of {} cobras observed a target in one or more exposures.".format(stats.Ncobras_used, len(self._Cobras)))
        print("{} cobras observed a target in all exposures.".format(stats.Ncobras_fully_used))
        return res


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
