from __future__ import print_function
from __future__ import absolute_import
import pulp
#from pulp import LpVariable, LpProblem, LpMinimize, GLPK, LpStatus, value, CPLEX, COIN, COIN_CMD, LOQO_CMD, GUROBI_CMD, GUROBI
from pulp import LpVariable, LpProblem, LpMinimize, GLPK, LpStatus, value, CPLEX, COIN, COIN_CMD, GUROBI_CMD, GUROBI
from numpy import random
from collections import namedtuple
from numpy import nan
import time

from . import datamodel as dm

def solve(prob, maxSeconds=5, solver='COIN_CMD'):
    #status = prob.solve()
    #status = prob.solve(GLPK(msg = 1, keepFiles=1))
    #status = prob.solve(COIN(msg = 1, keepFiles=1, maxSeconds=maxSeconds))
    if solver == 'COIN_CMD':
        status = prob.solve(COIN_CMD(msg = 1, keepFiles=1, maxSeconds=maxSeconds, threads=6, dual=10.))
    
    #status = prob.solve(LOQO_CMD(msg = 1, keepFiles=1, path='/Users/mxhf/ownCloudRZG/work/MPE/pfs/ETS/LOQO/loqo'))
    #status = prob.solve(CPLEX(msg = 1, keepFiles=1))
    #status = prob.solve(GUROBI_CMD(msg = 1, keepFiles=1))
    elif solver == 'GUROBI':
        status = prob.solve(GUROBI(msg = 1))
    else:
        print("ERROR: Unkown solver {}.".format(solver) )
        return None
    
    return status


def buildLPProblem(g, name="MinCostFlowTest", RSEP=5., cat='Integer'):
    start_time = time.time()

    prob = LpProblem(name, LpMinimize)
    flows = {}
    
    def addFlow(flows, n, l, u=None):
        f = LpVariable(n, l, u, cat=cat)
        #f = LpVariable(n, l, u, cat='continuous')
        flows[f.name] = f

    # add flow variables for target class to target arcs
    for tcid,tc in g.sciTargetClasses.items():
        for tid,t in tc.targets.items():
            addFlow(flows, "{}={}".format(tc.id,t.id), 0, 1) # capacity of one
    
    # add flow variables for target to target visits
    for tid,t in g.sciTargets.items():
        for visit in g.visits:
            addFlow(flows, r"{}={}_v{}".format(t.id,t.id,visit), 0, 1) # capacity of one
            #print r"{}={}_v{}".format(t.id,t.id,visit)

    # add flow variables for target visit to cobra visit arcs
    for aid,a in g.arcs.items():
        if type(a) == dm.TargetVisitToCobraVisitArc:
            addFlow(flows, "{}={}".format(a.startnode.id,a.endnode.id), 0,1)

    # add flow variables for cobra visit to cobra arcs
    for cid,c in g.cobras.items():
        for visit in g.visits:
            addFlow(flows, "{}_v{}={}".format(c.id,visit,c.id), 0, 1) 

    # add flow for overflow arc from targetClass node
    # mF: evetually marry with loop above, keep separate for readibility now.
    for tcid,tc in g.sciTargetClasses.items():
        addFlow(flows, "{}=SINK".format(tcid),0,1e6) 

    # add flow for overflow arcs from science target nodes
     # mF: evetually marry with loop above, keep separate for readibility now.
    for tid,t in g.sciTargets.items():
        addFlow(flows, "{}=SINK".format(tid),0,1e6) 
          
    # Now add constraints: At every intermediate node, inflow (* gain) = outflow
    for tcid,tc in g.sciTargetClasses.items():
        # for now set supply equal to number of targets in that target
        # class, i.e. Ideally we get them all observed.
        # If it is sufficient to observe a subset (i.e. N out of M) then this needs to be modified.
        S = tc.supply
        prob += pulp.lpSum([ flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in tc.outarcs]) == S
    
    # target nodes
    for tid,t in g.sciTargets.items():
        prob += pulp.lpSum( [ flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in t.inarcs]) * t.gain == \
            pulp.lpSum([ flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in t.outarcs])
    
    # target visit nodes
    for tv in g.targetVisits.values():
            prob += pulp.lpSum( [ flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in tv.inarcs]) == \
                sum([ flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in tv.outarcs])

    # cobra visit nodes
    for cvid,cv in g.cobraVisits.items():
        prob += pulp.lpSum( [ flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in cv.inarcs]) == \
            pulp.lpSum([ flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in cv.outarcs])
    
    ### for calibration targets
    INCLUDE_CALIB = True
    if INCLUDE_CALIB:
        # add flow variables for calib. target class to target arcs
        for tcid,tc in g.calTargetClasses.items():
            for tid,t in tc.targets.items():
                addFlow(flows, "{}={}".format(tc.id,t.id), 0, 1) # capacity of one
        
        # add flow for overflow arc from calib. targetClass node
        for tcid,tc in g.calTargetClasses.items():
            addFlow(flows, "{}=SINK".format(tcid),0,1e6) 
        
        # add flow for overflow arcs from calibb. target nodes
         # mF: eventually marry with loop above, keep separate for readibility now.
        #for tid,t in g.calTargets.iteritems():
        #    addFlow(flows, "{}=SINK".format(tid),0,1e6) 

        # Now add constraints: At every intermediate node, inflow (* gain) = outflow
        for tcid,tc in g.calTargetClasses.items():
            # for now set supply equal to number of targets in that target
            # class, i.e. Ideally we get them all observed.
            # If it is sufficient to observe a subset (i.e. N out of M) then this needs to be modified.
            S = tc.supply
            prob += pulp.lpSum([ flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in tc.outarcs]) == S
        
        # target nodes
        for tid,t in g.calTargets.items():
            prob += pulp.lpSum( [ flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in t.inarcs]) * t.gain == \
                pulp.lpSum([ flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in t.outarcs])
    ###
    


    print("Building cost equation ...")
    # attribute cost to the flows along the overflow arcs
    # Cost occures in two ways for now, either by 
    # flow occuring from a targetClass to the sink node (targets do not get observed at all)
    # or flow from a trget node to the sink node (target was only partially observed)
   
    cost = LpVariable("cost", 0) 
        

    prob += cost ==  pulp.lpSum([ a.cost * flows[a.id] for a in g.overflowArcs.values() ])\
        +  pulp.lpSum([ a.cost * flows[a.id] for a in g.targetToTargetVisitArcs.values() ]) \
        +  pulp.lpSum([ a.cost * flows[a.id] for a in g.targetVisitToCobraVisitArcs.values() ])
      

    # This sets the cost as objective function for the optimisation.
    prob += cost
    time_to_finish = time.time() - start_time
    print(" Time to completion: {:.2f} s".format(time_to_finish))
        
    return prob, flows, cost




def computeStats(g, flows, cost):
    Stats = namedtuple('Stats', \
        ['cost', 'NSciObs', 'NCalObs', 'NSciComplete', 'NCalComplete', 'Noverflow', 'Ncobras_used', 'Ncobras_fully_used'])

    NSciObs = 0
    NSciComplete = 0
    NCalObs = nan
    NCalComplete = nan
    NVISITS = len(g.visits)
    
    for t in g.sciTargets.values():
        NSciObs += value(sum([flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in t.inarcs]))
        NSciComplete += int(sum([ value(flows['{}={}'.format(a.startnode.id,a.endnode.id)]) for a in t.outarcs]) == t.gain) 

    #for t in g.sciTargets.itervalues():
    #    print t.id
    #    for a in t.outarcs:
    #        f = value(flows['{}={}'.format(a.startnode.id,a.endnode.id)])
    #        print " flow {}".format(f)
        
    Noverflow = 0
    for tcid,tc in g.sciTargetClasses.items():
        Noverflow += \
            value(flows['{}={}'.format(g.arcs["{}=SINK".format(tcid)].startnode.id,g.arcs["{}=SINK".format(tcid)].endnode.id)])

    Ncobras_used = 0
    Ncobras_fully_used = 0
    for c in g.cobras.values():
        v = value(sum([flows['{}={}'.format(a.startnode.id,a.endnode.id)] for a in c.inarcs]))
        Ncobras_used += int(v>0)
        Ncobras_fully_used += int(v == NVISITS)
        
    return Stats(value(cost), NSciObs, NCalObs, NSciComplete, NCalComplete, Noverflow, Ncobras_used, Ncobras_fully_used)

#stats = computeStats(g, flows, cost)
#print stats.NSciObs
#print stats.NSciComplete
