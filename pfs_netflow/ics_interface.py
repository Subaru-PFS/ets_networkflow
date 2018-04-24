import numpy as np

import time

from collections import OrderedDict

from scipy.spatial.distance import cdist

from gurobipy import quicksum

from ics.cobraOps.Bench import Bench
from ics.cobraOps.CobrasCalibrationProduct import CobrasCalibrationProduct
from ics.cobraOps.CollisionSimulator import CollisionSimulator
from ics.cobraOps.DistanceTargetSelector import DistanceTargetSelector
from ics.cobraOps.RandomTargetSelector import RandomTargetSelector
from ics.cobraOps.TargetGroup import TargetGroup


def getElbowPositions(g, ivisibilities):
    print("Computing potential cobra elbow positions ...")
    # Find elbow focal plane positions for all possible cobra to target assignments
    elbowPositions = OrderedDict()

    for pid in ivisibilities:
        # convert cobra positions to array of complex numbers
        cids = [ c.id for c in g.cobras.values()]
        cpos = np.array( [ complex(g.cobras[cid].getX(pid) , g.cobras[cid].getY(pid)) for cid in cids] )
        # and create bench object
        bench = Bench(cobraCenters=np.array( cpos ))
        print("Number of cobras:", bench.cobras.nCobras)
    
        elbowPositions[pid] = OrderedDict()

        for i, cid in enumerate(cids):
            if not cid[2:] in ivisibilities[pid]:
                continue

            ttid = ivisibilities[pid][cid[2:] ]
            ttid

            for tid in ttid:
                stid = "T_" + tid
                if not stid in g.targets:
                    # in case it is a calibration target
                    stid = "T_{}_v{}".format(tid, pid)
                    if not stid in g.targets:
                        print("Did not find target {}".format(stid))
                        continue
                t = g.targets[stid]
                x = t.getX(pid)
                y = t.getY(pid)


                targetPositionsForCobra1 = [complex(x,y)]
                ebp = bench.cobras.calculateCobraElbowPositions(cobraIndex=i, fiberPositions=targetPositionsForCobra1)[0]

                fpid = "T_{}_v{}={}_v{}".format(tid,pid,cid,pid)
                elbowPositions[pid][fpid] = ebp.real,ebp.imag
    return elbowPositions



def compute_collision_flow_pairs(g, elbowPositions, fiber_collision_radius = 1.):
    print("Finding fiber-fiber, fiber-elbow, and elbow-elbow collisions ...")
    collision_flow_pairs  = []

    for pid in g.pointings:

        # target x/y positions
        txx  =  [t.getX(pid) for tid, t in g.targets.items()] 
        tyy  =  [t.getY(pid) for tid, t in g.targets.items()] 


        tfid = []
        for tid, t in g.targets.items():

            if tid in g.calTargets:
                fid = g.calTargets[tid].inarcs[0].id
                tfid.append( fid )
            elif tid in g.sciTargets:
                fid = "{}={}_v{}".format(tid,tid,pid)
                tfid.append( fid )
            else:
                #continue # this target is not part of the problem, probably did not survive RMAX cut
                fid = "NAN"
                print("Error, target {} not part of the problem.".format(tid))


        # elbow x/y positions
        ebxx = [ ebp[0] for ebp in elbowPositions[pid].values()]
        ebyy = [ ebp[1] for ebp in elbowPositions[pid].values()]
        ebid = [ eid for eid in elbowPositions[pid]]

        xx = txx + ebxx
        yy = tyy + ebyy
        ID = tfid + ebid

        #N = len(ID)
        points = list( zip(xx,yy) )
        Y = cdist( points, points )

        # any target separation that is smaller than 2 x the collision radius will be flagged a s collision
        cc = Y <= (fiber_collision_radius*2.) 

        N = len(xx)
        ncoll = int( (np.sum(cc.flatten()) - N)/2. )

        print ("Pointing {}: Found  {:d} collision pairs.".format( pid, ncoll  ))

        # identify collision pairs
        # array of indices
        ii = np.arange( N )
        for i in range(cc.shape[0]):
            x1,y1 =  xx[i], yy[i]
            # only iterate over the indices that are colliding and the upper diagonal in the collision matrix
            jj = ii[ cc[i,:] * ii > i ] 
            for j in jj: 
                if cc[i,j]:
                    x2,y2 =  xx[j], yy[j]
                    collision_flow_pairs.append( (ID[i], ID[j]) )

    return collision_flow_pairs


def addCollisionFlowConstraints(model, flows, collision_flow_pairs):
    print("Adding flow constraints to model ...")
    
    i = 0
    for fp in collision_flow_pairs:
        if fp[0] in flows and fp[1] in flows:
            i += 1
            model.addConstr( quicksum( [ flows[ fp[0] ], flows[ fp[1] ] ] ) <= 1. , "{}_OR_{}".format(fp[0],fp[1]))
        elif not fp[0] in flows:
            print("No flow variable for {}".format(fp[0]))
        elif not fp[1] in flows:
            print("No flow variable for {}".format(fp[1]))
    print("Added {} collision flow constraints.".format(i))
    
    
def getBench(g, pid):
    # convert cobra positions to array of complex numbers
    cids = [  c.id for c in g.cobras.values() ]
    
    cpos = np.array( [ complex(g.cobras[cid].getX(pid) , g.cobras[cid].getY(pid)) for cid in cids] )
    # and create bench object
    bench = Bench(cobraCenters=np.array( cpos ))
    print("Number of cobras:", bench.cobras.nCobras)
    return bench


def getCollSim(g, pid, bench, PLOT=False, PLOTTRAJ=False, ZOOMIN=False):
    # initialize assigment list with Null
    # as required by the collision simulator code.
    assignments = OrderedDict()
    for c in g.cobras.values(): 
        assignments[c.id] = ('NULL',0.+0.j)

    # Now find which cobras have been assigned to which targets. 
    # Loop over all targetPointing to CobraPointing arcs.Filter for those that correspond to the current pointing.
    for a in filter(lambda x : x.pointing == pid, g.targetPointingToCobraPointingArcs.values()):
        if a.flow > 0.:
            t = a.startnode.target
            c = a.endnode.cobra 
            #assignments[c.id] = t.id
            assignments[c.id] = (t.id, complex(t.getX(pid), t.getY(pid)) )

    # convert cobra positions to array of complex numbers
    cids = [ cid for cid in assignments]
    cpos = np.array( [ complex(g.cobras[cid].getX(pid) , g.cobras[cid].getY(pid)) for cid in assignments] )
    # and create bench object
    #bench = Bench(cobraCenters=np.array( cpos ))
    print("Number of cobras:", bench.cobras.nCobras)

    # Create list of fiber to cobra assignments
    targetpos       = np.array([t[1] for t in assignments.values()])
    ids             = np.array([t[0] for t in assignments.values()])
    selectedTargets = TargetGroup(positions=targetpos, ids=ids)
    print("Number of  targets:", selectedTargets.nTargets)


    # Simulate an observation
    start = time.time()
    simulator = CollisionSimulator(bench, selectedTargets)
    simulator.run()
    print("Number of cobras involved in collisions:", simulator.nCollisions)
    print("Number of colliding cobras that are not also end collisions: ", simulator.nCollisions - simulator.nEndPointCollisions)
    print("Total simulation time (s):", time.time() - start)

    if PLOT:
        # Plot the simulation results
        __ = simulator.plotResults(paintFootprints=False, size=(10,10))

        W = 200
        xc = -75.
        yc =  20.
        xmin,xmax = -W/2.+xc,W/2.+xc
        ymin,ymax = -W/2.+yc,W/2.+yc

        #for i,(x,y) in enumerate(zip( np.real( targetpos ), np.imag( targetpos ) ) ):
        #    if xmin < x and x < xmax and ymin < y and y < ymax:
        #        plt.text(x,y,ids[i], size=12)

        if ZOOMIN:
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])



    if PLOTTRAJ:
        # Animate one of the trajectory collisions
        (problematicCobras,) = np.where(np.logical_and(simulator.collisions, simulator.endPointCollisions == False))

        if len(problematicCobras) > 0:
            simulator.animateCobraTrajectory(problematicCobras[0])
        else:
            print("No trajectory collisions detected...")

        # Pause the execution to have time to inspect the figures
        plotUtils.pauseExecution()
        
    return simulator, assignments, cids
    
