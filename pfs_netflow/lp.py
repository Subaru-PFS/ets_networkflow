from __future__ import print_function
from __future__ import absolute_import
import pulp
from collections import namedtuple
from collections import OrderedDict
import numpy as np
import time
from . import datamodel as dm
from pfs_netflow.utils import pp 
from gurobipy import quicksum
from pfs_netflow.survey_plan import compute_collision_pairs, compute_collision_flow_pairs

from pfs_netflow import ics_interface
from pfs_netflow import survey_plan


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
                # target to target pointing arcs
                for arc2 in arc.endnode.outarcs:
                    if arc2.endnode.id == "SINK":
                        addFlow(flows, arc2.id, 0, 1e6)
                    else:
                        addFlow(flows, arc2.id, 0, 1)  # capacity of one

    # add flow variables for target pointing to cobra pointing arcs
    for aid, a in g.arcs.items():
        if type(a) == dm.TargetPointingToCobraPointingArc:
            addFlow(flows, a.id, 0, 1)

    # add flow variables for cobra pointing to cobra arcs
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

    # target pointing nodes
    for tv in g.targetPointings.values():
            prob += pulp.lpSum([flows[a.id] for a in tv.inarcs]) == \
                sum([flows[a.id] for a in tv.outarcs])

    # cobra pointing nodes
    for cvid, cv in g.cobraPointings.items():
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
        + pulp.lpSum([a.cost * flows[a.id] for a in g.targetToTargetPointingArcs.values()]) \
        + pulp.lpSum([a.cost * flows[a.id] for a in g.targetPointingToCobraPointingArcs.values()])

    # This sets the cost as objective function for the optimisation.
    prob += cost
    time_to_finish = time.time() - start_time
    print(" Time to completion: {:.2f} s".format(time_to_finish))
    return prob, flows, cost


def buildLPProblemGRB(g, name="MinCostFlowTest"):
    """
    Build the LP problem using guropipy's native interface.
    """
    from gurobipy import Model, quicksum, GRB
    
    NCobras = len(g.cobras)
    NSciTargets = len(g.sciTargets)
    NCalTargets = len(g.calTargets)
    NPOINTINGS = len( g.pointings )

    summary = ""
    summary += pp("Building LP problem ...")
    summary += pp("NPOINTINGS = {}".format(NPOINTINGS))
    summary += pp("Searching optimal strategy to observe in ")
    summary += pp(" {} pointings".format(NPOINTINGS))
    summary += pp(" {} science targets".format(NSciTargets))
    summary += pp(" {} calib. targets".format(int(NCalTargets/NPOINTINGS) ))
    summary += pp(" {} cobras".format(NCobras))
    summary += pp("num nodes: {}".format(len(g.nodes)))
    summary += pp("num edges: {}".format(len(g.arcs)))
    
    start_time = time.time()

    m = Model(name)

    flows = {}

    def addFlow(m, flows, n, l=0, u=GRB.INFINITY, vtype=GRB.INTEGER):
        f = m.addVar(vtype=vtype, name=n, lb=l, ub=u)
        flows[n] = f

    # add flow variables for target class to target arcs
    for tcid, tc in g.sciTargetClasses.items():
        for arc in tc.outarcs:
            if arc.endnode.id == "SINK":
                addFlow(m, flows, arc.id, 0, 1e6)
            else:
                addFlow(m, flows, arc.id, vtype=GRB.BINARY)  # capacity of one
                # target to target pointing arcs
                for arc2 in arc.endnode.outarcs:
                    if arc2.endnode.id == "SINK":
                        addFlow(m, flows, arc2.id, 0, 1e6)
                    else:
                        addFlow(m, flows, arc2.id, vtype=GRB.BINARY)  # capacity of one
                        
    # add flow variables for target pointing to cobra pointing arcs
    for aid, a in g.arcs.items():
        if type(a) == dm.TargetPointingToCobraPointingArc:
            addFlow(m, flows, a.id, vtype=GRB.BINARY)

    # add flow variables for cobra pointing to cobra arcs
    for cid, c in g.cobras.items():
        for arc in c.inarcs:
            addFlow(m, flows, arc.id, vtype=GRB.BINARY)
            
            
    # Now add constraints: At every intermediate node, inflow (* gain) = outflow
    for tcid, tc in g.sciTargetClasses.items():
        # for now set supply equal to number of targets in that target
        # class, i.e. Ideally we get them all observed.
        # If it is sufficient to observe a subset (i.e. N out of M) then this needs to be modified.
        S = tc.supply
        m.addConstr( quicksum([flows[a.id] for a in tc.outarcs]) == S )

    # target nodes
    for tid, t in g.sciTargets.items():
        m.addConstr( quicksum( [ flows[a.id] for a in t.inarcs]) * t.gain == \
                     quicksum( [ flows[a.id] for a in t.outarcs])\
                   )

    # target pointing nodes
    for tv in g.targetPointings.values():
            m.addConstr( quicksum([flows[a.id] for a in tv.inarcs]) == \
                         quicksum([flows[a.id] for a in tv.outarcs])\
                       )

            
    # cobra pointing nodes
    for cvid, cv in g.cobraPointings.items():
        m.addConstr( quicksum([flows[a.id] for a in cv.inarcs]) == \
                     quicksum([flows[a.id] for a in cv.outarcs])\
                   )

    # for calibration targets
    INCLUDE_CALIB = True
    if INCLUDE_CALIB:
        # add flow variables for calib. target class to target arcs
        for tcid, tc in g.calTargetClasses.items():
            for arc in tc.outarcs:
                if arc.endnode.id == "SINK":
                    addFlow(m, flows, arc.id, 0, 1e6)
                else:
                    addFlow(m, flows, arc.id, vtype=GRB.BINARY)  # capacity of one

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
        m.addConstr( quicksum([flows[a.id] for a in tc.outarcs]) == S )

    # target nodes
    for tid, t in g.calTargets.items():
        m.addConstr( quicksum([flows[a.id] for a in t.inarcs]) * t.gain == \
                     quicksum([flows[a.id] for a in t.outarcs])\
                   )
        
    print("Building cost equation ...")
    # attribute cost to the flows along the overflow arcs
    # Cost occurs in two ways for now, either by
    # flow occuring from a targetClass to the sink node (targets do not get observed at all)
    # or flow from a trget node to the sink node (target was only partially observed)
    cost = m.addVar(vtype=GRB.CONTINUOUS, name="cost")

    m.addConstr( cost == quicksum([a.cost * flows[a.id] for a in g.overflowArcs.values()])\
        + quicksum([a.cost * flows[a.id] for a in g.targetToTargetPointingArcs.values()]) \
        + quicksum([a.cost * flows[a.id] for a in g.targetPointingToCobraPointingArcs.values()])
               )
    
    # This sets the cost as objective function for the optimisation.
    m.setObjective(cost, GRB.MINIMIZE)
    time_to_finish = time.time() - start_time
    time_to_build = time.time() - start_time
    summary += pp("Time to build model: {:.4e} s".format(time_to_build))


    # bBldly add flows and costs to the survey plan graph.
    # reaally should be part of the LP problem, but gurobipy's Model
    # is immutable.
    g.flows = flows
    g.cost = cost
    
    return m


def add_fiber_fiber_collision_avoidance(g, model, pointings, target_fplane_pos):
    """
    Adds flow constrint equations to avoid fiber to fiber endpoint collision "a-priori",
    i.e. solutions which cause collisions won't be considered in the first place.
    """
    collision_pairs = compute_collision_pairs(pointings, target_fplane_pos) 

    flow_pairs = compute_collision_flow_pairs(g, collision_pairs)
    print("Adding {} collision avoidance constraints.".format(len(flow_pairs)))              
    for fp in flow_pairs:
        if fp[0] in g.flows and fp[1] in g.flows:
            model.addConstr( quicksum( [ g.flows[ fp[0] ], g.flows[ fp[1] ] ] ) <= 1. )
            

def add_fiber_elbow_collision_avoidance(g, model, ivisibilities, fiber_collision_radius = 1.0):
    """
    Adds flow constrint equations to avoid fiber to fiber AND fiber to elbow endpoint collision "a-priori",
    i.e. solutions which cause collisions won't be considered in the first place.
    1.074 is good radius to also avoid grazing collisions at link.
    """
    elbowPositions = ics_interface.getElbowPositions(g, ivisibilities)

    collision_flow_pairs = \
    ics_interface.compute_collision_flow_pairs(g, elbowPositions, fiber_collision_radius = fiber_collision_radius) 
    ics_interface.addCollisionFlowConstraints(model, g.flows, collision_flow_pairs)


def optimize(model, g, SOLVE_TIME_LIMIT, NUMER_OF_THREDS, MIP_GAP):
    # Solve problem!
    summary = ""
    summary += pp("Solving LP problem ...")
    start_time = time.time()

    model.Params.timelimit = SOLVE_TIME_LIMIT  # time limit in seconds
    model.Params.Threads = NUMER_OF_THREDS     # number of thread to use
    model.Params.MIPGap = MIP_GAP              # solution tolerance 

    model.optimize()

    time_to_solve = time.time() - start_time
    summary += pp("Solve status is [{}].".format( model.status ))
    summary += pp("Time to solve: {:.4e} s".format(time_to_solve))

    survey_plan.setflows(model, g)