from __future__ import print_function
from __future__ import absolute_import
import pulp
from collections import namedtuple
from collections import OrderedDict
import numpy as np
import time
from . import datamodel as dm


def solve(prob, maxSeconds=5, solver='COIN_CMD'):
    if solver == 'COIN_CMD':
        status = prob.solve(pulp.COIN_CMD(msg=1, keepFiles=1,
                                          maxSeconds=maxSeconds,
                                          threads=6, dual=10.))

    elif solver == 'GUROBI':
        status = prob.solve(pulp.GUROBI(msg=1))
    else:
        print("ERROR: Unkown solver {}.".format(solver))
        return None
    return status


def buildLPProblem(g, name="MinCostFlowTest", cat='Integer'):
    start_time = time.time()

    prob = pulp.LpProblem(name, pulp.LpMinimize)
    flows = {}

    def addFlow(flows, n, l, u=None):
        f = pulp.LpVariable(n, l, u, cat=cat)
        flows[f.name] = f

    # add flow variables for target class to target arcs
    for tcid, tc in g.sciTargetClasses.items():
        for arc in tc.outarcs:
            if arc.endnode.id == "SINK":
                addFlow(flows, arc.id, 0, 1e6)
            else:
                addFlow(flows, arc.id, 0, 1)  # capacity of one
                # target to target visit arcs
                for arc2 in arc.endnode.outarcs:
                    if arc2.endnode.id == "SINK":
                        addFlow(flows, arc2.id, 0, 1e6)
                    else:
                        addFlow(flows, arc2.id, 0, 1)  # capacity of one

    # add flow variables for target visit to cobra visit arcs
    for aid, a in g.arcs.items():
        if type(a) == dm.TargetVisitToCobraVisitArc:
            addFlow(flows, a.id, 0, 1)

    # add flow variables for cobra visit to cobra arcs
    for cid, c in g.cobras.items():
        for arc in c.inarcs:
            addFlow(flows, arc.id, 0, 1)

    # Now add constraints: At every intermediate node, inflow (* gain) = outflow
    for tcid, tc in g.sciTargetClasses.items():
        # for now set supply equal to number of targets in that target
        # class, i.e. Ideally we get them all observed.
        # If it is sufficient to observe a subset (i.e. N out of M) then this needs to be modified.
        S = tc.supply
        prob += pulp.lpSum([flows[a.id] for a in tc.outarcs]) == S

    # target nodes
    for tid, t in g.sciTargets.items():
        prob += pulp.lpSum( [ flows[a.id] for a in t.inarcs]) * t.gain == \
            pulp.lpSum([ flows[a.id] for a in t.outarcs])

    # target visit nodes
    for tv in g.targetVisits.values():
            prob += pulp.lpSum([flows[a.id] for a in tv.inarcs]) == \
                sum([flows[a.id] for a in tv.outarcs])

    # cobra visit nodes
    for cvid, cv in g.cobraVisits.items():
        prob += pulp.lpSum([flows[a.id] for a in cv.inarcs]) == \
            pulp.lpSum([flows[a.id] for a in cv.outarcs])

    # for calibration targets
    INCLUDE_CALIB = True
    if INCLUDE_CALIB:
        # add flow variables for calib. target class to target arcs
        for tcid, tc in g.calTargetClasses.items():
            for arc in tc.outarcs:
                if arc.endnode.id == "SINK":
                    addFlow(flows, arc.id, 0, 1e6)
                else:
                    addFlow(flows, arc.id, 0, 1)  # capacity of one

        # add flow for overflow arcs from calibb. target nodes
        #  mF: eventually marry with loop above, keep separate for readibility now.
        # for tid,t in g.calTargets.iteritems():
        #    addFlow(flows, "{}=SINK".format(tid),0,1e6)

        # Now add constraints: At every intermediate node, inflow (* gain) = outflow
        for tcid, tc in g.calTargetClasses.items():
            # for now set supply equal to number of targets in that target
            # class, i.e. Ideally we get them all observed.
            # If it is sufficient to observe a subset (i.e. N out of M) then this needs to be modified.
            S = tc.supply
            prob += pulp.lpSum([flows[a.id] for a in tc.outarcs]) == S

        # target nodes
        for tid, t in g.calTargets.items():
            prob += pulp.lpSum([flows[a.id] for a in t.inarcs]) * t.gain == \
                pulp.lpSum([flows[a.id] for a in t.outarcs])

    print("Building cost equation ...")
    # attribute cost to the flows along the overflow arcs
    # Cost occurs in two ways for now, either by
    # flow occuring from a targetClass to the sink node (targets do not get observed at all)
    # or flow from a trget node to the sink node (target was only partially observed)
    cost = pulp.LpVariable("cost", 0)

    prob += cost == pulp.lpSum([a.cost * flows[a.id] for a in g.overflowArcs.values()])\
        + pulp.lpSum([a.cost * flows[a.id] for a in g.targetToTargetVisitArcs.values()]) \
        + pulp.lpSum([a.cost * flows[a.id] for a in g.targetVisitToCobraVisitArcs.values()])

    # This sets the cost as objective function for the optimisation.
    prob += cost
    time_to_finish = time.time() - start_time
    print(" Time to completion: {:.2f} s".format(time_to_finish))
    return prob, flows, cost

