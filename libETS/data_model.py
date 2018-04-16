import numpy as np
import pycconv
import pyETS
from collections import OrderedDict, defaultdict


def _get_visibility(cobras, tpos):
    cbr = [[c.center, c.innerLinkLength, c.outerLinkLength, c.dotcenter,
            c.rdot] for c in cobras]

    return pyETS.getVis(tpos, cbr)


def _build_network(cobras, targets, tpos, classdict, tvisit, vis_cost=None, gurobi=False):
    Cv_i = defaultdict(list)  # Cobra visit inflows
    Tv_o = defaultdict(list)  # Target visit outflows
    Tv_i = defaultdict(list)  # Target visit inflows
    T_o = defaultdict(list)  # Target outflows (only science targets)
    T_i = defaultdict(list)  # Target inflows (only science targets)
    CTCv_o = defaultdict(list)  # Calibration Target class visit outflows
    STC_o = defaultdict(list)  # Science Target outflows
    if gurobi:
        import gurobipy as gbp
        prob = gbp.Model("problem")
        lpSum = gbp.quicksum
        def add_constraint(problem, constraint):
            problem.addConstr(constraint)
        def varValue(var):
            return var.X
    else:
        import pulp
        prob = pulp.LpProblem("problem", pulp.LpMinimize)
        cost = pulp.LpVariable("cost", 0)
        lpSum = pulp.lpSum
        def add_constraint(problem, constraint):
            problem += constraint
        def varValue(var):
            return pulp.value(var)

    nreqvisit = []
    for t in targets:
        if isinstance(t, ScienceTarget):
            nreqvisit.append(int(t.obs_time/tvisit))
        else:
            nreqvisit.append(0)

    vis = [_get_visibility(cobras, tp) for tp in tpos]
    nvisits = len(vis)

    if vis_cost is None:
        vis_cost = [0.]*nvisits

    def newvar(lo, hi):
        newvar._varcount += 1
        if gurobi:
            if lo == 0 and hi == 1:
               return prob.addVar(vtype=gbp.GRB.BINARY, name="v{}".format(newvar._varcount))
            else:
               return prob.addVar(vtype=gbp.GRB.INTEGER, name="v{}".format(newvar._varcount), lb=lo, ub=hi)
        else:
            return pulp.LpVariable("v{}".format(newvar._varcount), lo, hi,
                                   cat=pulp.LpInteger)
    newvar._varcount = 0

    # define LP variables

    # create nodes for every visit and calibration target class
    for key, value in classdict.items():
        if value["calib"]:
            for ivis in range(nvisits):
                f = newvar(0,100000000)
                CTCv_o[(key, ivis)].append(f)
                cost += f*value["nonObservationCost"]

    for ivis in range(nvisits):
        for tidx, val in vis[ivis].items():
            tgt = targets[tidx]
            TC = tgt.targetclass
            Class = classdict[TC]
            if isinstance(tgt, ScienceTarget):
                # Target node to target visit node
                f = newvar(0, 1)
                T_o[tidx].append(f)
                Tv_i[(tidx, ivis)].append(f)
                if len(T_o[tidx]) == 1:  # freshly created
                    # Science Target class node to target node
                    f = newvar(0, 1)
                    T_i[tidx].append(f)
                    STC_o[TC].append(f)
                    if len(STC_o[TC]) == 1:  # freshly created
                        # Science Target class node to sink
                        f = newvar(0, None)
                        STC_o[TC].append(f)
                        cost += f*Class["nonObservationCost"]
                    # Science Target node to sink
                    f = newvar(0, None)
                    T_o[tidx].append(f)
                    cost += f*Class["partialObservationCost"]
            elif isinstance(tgt, CalibTarget):
                # Calibration Target class node to target visit node
                f = newvar(0, 1)
                Tv_i[(tidx, ivis)].append(f)
                CTCv_o[(TC, ivis)].append(f)
            for cidx in val:
                # target visit node to cobra visit node
                f = newvar(0, 1)
                Cv_i[(cidx, ivis)].append(f)
                Tv_o[(tidx, ivis)].append((f, cidx))
                cost += f*vis_cost[ivis]

    # Cost function
    prob += cost

    # Constraints

    # every Cobra can observe at most one target per visit
    for inflow in Cv_i.values():
        add_constraint(prob, lpSum([f for f in inflow]) <= 1)

    # every calibration target class must be observed a minimum number of times
    # every visit
    for key, value in CTCv_o.items():
        add_constraint(prob, lpSum([v for v in value]) >= \
                classdict[key[0]]["numRequired"])

    # inflow and outflow at every Tv node must be balanced
    for key, ival in Tv_i.items():
        oval = Tv_o[key]
        add_constraint(prob, lpSum([v for v in ival]+[-v[0] for v in oval]) == 0)

    # inflow and outflow at every T node must be balanced
    for key, ival in T_i.items():
        oval = T_o[key]
        nvis = nreqvisit[key]
        add_constraint(prob, lpSum([nvis*v for v in ival]+[-v for v in oval]) == 0)

    # Science targets must be either observed or go to the sink
    for key, val in STC_o.items():
        add_constraint(prob, lpSum([v for v in val]) == len(val)-1)

    if gurobi:
        prob.ModelSense = 1  # minimize
        model.optimize()
    else:
        status = prob.solve(pulp.COIN_CMD(msg=1, keepFiles=0, maxSeconds=100,
                                          threads=1, dual=10.))

    res = [{} for _ in range(nvisits)]
    for k1, v1 in Tv_o.items():
        for i2 in v1:
            visited = varValue(i2[0]) > 0
            if visited:
                tidx, ivis = k1
                cidx = i2[1]
                res[ivis][tidx] = cidx
    return res

def observeWithNetflow(cbr, tgt, tpos, classdict, tvisit, vis_cost=None):
    return _build_network(cbr, tgt, tpos, classdict, tvisit, vis_cost)


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

    @property
    def cobraCollisionRadius(self):
        """the radius of a Cobra tip : float
        This is used for collision detection. It could also become a property
        of each individual Cobra, but apparently this quantity is not very
        variable.
        """
        return self._cobraCollisionRadius

    def subtract_obs_time(self, tgt, obs, tvisit):
        res = []
        for k in obs.keys():
            tgt[k].reduce_time(tvisit)
        for t in tgt:
            if t.obs_time > 0.:
                res.append(t)
        return res

    def get_fp_positions(self, tgt):
        return [t.fp_position(self._ra, self._dec, self._posang, self._time)
                for t in tgt]

#    def observeWithNetflow(self, tgt, classdict, nvisit, tvisit, vis_cost=None):
#        for t in tgt:
#            t.calc_position(self._ra, self._dec, self._posang, self._time)
#        return _build_network(self._Cobras, tgt, classdict, nvisit, tvisit, vis_cost)


class Target(object):
    """Base class for all target types observable with PFS. All targets are
    initialized with RA/Dec and an ID string. From RA/Dec, the target can
    determine its position on the focal plane, once also the telescope attitude
    is known. """
    def __init__(self, ID, ra, dec, targetclass):
        self._ID = str(ID)
        self._ra = float(ra)
        self._dec = float(dec)
        self._targetclass = targetclass

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

    def fp_position(self, raTel, decTel, posang, time):
        return pycconv.cconv([self._ra], [self._dec],
                             raTel, decTel, posang, time)[0]

    @property
    def targetclass(self):
        """string representation of the target's class"""
        return self._targetclass


class ScienceTarget(Target):
    """Derived from the Target class, with the additional attributes priority
    and observation time.

    All different types of ScienceTarget need to be derived from this class."
    """
    def __init__(self, ID, ra, dec, obs_time, pri, prefix):
        super(ScienceTarget, self).__init__(ID, ra, dec,
                                            "{}_P{}".format(prefix, pri))
        self._obs_time = float(obs_time)
        self._pri = int(pri)

    @property
    def obs_time(self):
        """required observation time : float"""
        return self._obs_time

    def reduce_time(self, dt):
        self._obs_time -= float(dt)


class CalibTarget(Target):
    """Derived from the Target class."""


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


def readScientificFromFile(file, prefix):
    with open(file) as f:
        res = []
        ll = f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_, ra, dec, tm, pri = (
                    str(tt[0]), float(tt[1]), float(tt[2]),
                    float(tt[3]), int(tt[4]))
                res.append(ScienceTarget(id_, ra, dec, tm, pri, prefix))
    return res


def readCalibrationFromFile(file, targetclass):
    with open(file) as f:
        res = []
        ll = f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_, ra, dec = (str(tt[0]), float(tt[1]), float(tt[2]))
                res.append(CalibTarget(id_, ra, dec, targetclass))
    return res


def getFullFocalPlane():
    cobras = pyETS.getAllCobras()
    res = []
    for i in range(len(cobras)):
        ID = "{}".format(i)
        res.append(Cobra(ID, cobras[i][0], cobras[i][3], cobras[i][4],
                         cobras[i][1], cobras[i][2]))
    return res
