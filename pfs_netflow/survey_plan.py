from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import range
import numpy as np
from collections import OrderedDict
import pulp

from astropy.table import Table, Column

from scipy.spatial.distance import cdist
import time

from . import datamodel as dm


def buildSurveyPlan(cobras, targets, nreqv_dict, visibilities, class_dict,
                    cost_dict, supply_dict, RMAX, CENTER=(0.,0.), COBRAS=[]):
    """
    Builds a graph which represents a survey plan.
    

    Args:
        cobras (OrderedDict): A dictionary of cobra x,y in the focal plane. The key is the cobra ID.
        targets (OrderedDict): A dictionary of target x,y in the focal plane. The key is the target ID.
        nreqv_dict (OrderedDict): Number of required pointings per target. The key is the target ID, elements are int.
        visibilities (OrderedDict): Dictionary describing which target can be observed by which cobra. This is
                                    a dictionary of dictionaries, one for each pointing. 
        class_dict (OrderedDict): Dictionary to which specifies the target class per target ID.
        cost_dict (OrderedDict): Dictionary specifing the cost function (cost of non-observation and more).
        supply_dict (OrderedDict): Dictionary specifing the supply (number targets in each target class that must be observed).
        RMAX (float): Limits the problem to a specific maximum radius in the focal plane.
        CENTER (2-tuple of float): Limits the problem to a specific center in the focal plane.
        COBRAS (list of CID): Limits the problem to a specific list of cobras in the focal plane.
        

    Returns:
        datamodel.SurveyPlan: The survey graph.
    """
    
    print("buildSurveyPlan")


    # Check if a move distance dependant cost function was given. If not set to zero.
    if 'cobra_move' not in cost_dict:
        # no specific cost function for the cobra move distance is given. Set to zero cost.
        costOfD = lambda x : 0.
    # MR FIXME: "else" missing?
    costOfD = cost_dict['cobra_move']


    # generate a directed graph
    g = dm.SurveyPlan()
    g.pointings = list(visibilities)

    # Add global sink node
    T = dm.Sink()
    g.add_node(T)

    # Add nodes for the cobras and the cobra pointings and arcs between them
    for cid in cobras:
        x, y = cobras[cid]
        if (x-CENTER[0])**2 + (y-CENTER[1])**2 > RMAX**2:
            continue
        if COBRAS != [] and cid not in COBRAS:
            continue
        c = dm.Cobra(cid, fplane_pos=(x,y))
        g.add_node(c)
        # replicate node as many times as there are pointings
        for pid in visibilities:
            cv = dm.CobraPointing(cid=cid, cobra=c, pointing=pid)
            g.add_node(cv)
            g.add_arc(dm.CobraPointingToCobraArc(cv, c))

    # Add nodes for target classes
    for tc in np.unique(list(class_dict.values())):
        # add target-class super node
        if tc.startswith("sci_"):
            targetClass = dm.SciTargetClass(tc)
            targetClass.cost = cost_dict[tc][0]
            targetClass.cost_partial_compl = cost_dict[tc][1]
            targetClass.supply = supply_dict[tc]
            g.add_node(targetClass)

            # Add costly overflow arc to the target class
            g.add_arc(dm.OverflowArc(targetClass.cost, targetClass, T))
        elif tc.startswith("cal_") or tc.startswith("sky_"):
            for pid in g.pointings:
                targetClass = dm.CalTargetClass(tc, pid)
                targetClass.cost = cost_dict[tc]
                targetClass.supply = supply_dict[tc]
                g.add_node(targetClass)

                # Add costly overflow arc to the target class
                g.add_arc(dm.OverflowArc(targetClass.cost, targetClass, T))
        else:
            print("Error unknown target class {}".format(tc))


    # Add nodes for the targets and the target pointings and arcs between them
    for tid in targets:
        tc = class_dict[tid]
        fplane_positions = targets[tid]
        nrv = nreqv_dict[tid]
        
        txx = np.array( [fplane_positions[pid][0] for pid in fplane_positions] )
        tyy = np.array( [fplane_positions[pid][1] for pid in fplane_positions] )
        
        if ( (txx - CENTER[0])**2 + (tyy - CENTER[1])**2 > RMAX**2 ).all():
            continue
        if tc.startswith("sci_"):
            t = dm.SciTarget(tid, fplane_positions=fplane_positions)
            t.target_class = tc
            t.gain = nrv
            g.add_node(t)
            
            # Add as many TargetPointing nodes for this target as there are pointings
            for pid in visibilities:
                tv = dm.TargetPointing(tid, target=t, pointing=pid)
                g.add_node(tv)
                ttva = dm.TargetToTargetPointingArc(t, tv)
                # Here we assign the cost for the respective pointing
                # increasing the cost for later pointings encourages
                # earlier observation.
                ttva.cost = cost_dict["pointings"][pid]
                g.add_arc(ttva)

            # Add arc from target class to target
            targetClass = g.sciTargetClasses[dm.SciTargetClass.getID(tc)]
            targetClass.add_target(t)
            g.add_arc(dm.TargetClassToTargetArc(targetClass, t))

            # Add VERY costly overflow arc from targets to sink
            g.add_arc(dm.OverflowArc(targetClass.cost_partial_compl, t, T))

        elif tc.startswith("cal_") or tc.startswith("sky_"):
            for pid in visibilities:
                # Add as many nodes for each calibration target as there are pointings
                if tc.startswith("cal_"):
                    t = dm.StarCalTarget(tid, fplane_positions=fplane_positions, pointing=pid)
                else:
                    t = dm.SkyCalTarget(tid, fplane_positions=fplane_positions, pointing=pid)
                #t = dm.CalTarget(tid, fplane_positions=fplane_positions, pointing=pid)
                t.target_class = tc
                t.gain = 1
                g.add_node(t)

                # Add arc from target class to target
                targetClass = g.calTargetClasses[dm.CalTargetClass.getID(tc, pid)]
                targetClass.add_target(t)
                g.add_arc(dm.TargetClassToTargetArc(targetClass, t))

    # Calibration targets only need to be observed once.
    # Remember we replicate calibration targets such that there is one per
    # observation.
    for tid, t in g.calTargets.items():
        t.gain = 1.

    # Set supply for each science target class.
    for tcid, targetClass in g.sciTargetClasses.items():
        if targetClass.supply == np.inf:
            # Rather than using infinite supply (which woudl result in infinite cost)
            # set supply for the targetClasses to number of targets.
            # I.e. we want to get them all observed (if possible).
            targetClass.supply = len(targetClass.targets)

    # Set supply for each calibration target class.
    for tcid, targetClass in g.calTargetClasses.items():
        if targetClass.supply == np.inf:
            targetClass.supply = len(targetClass.targets)

    for pointing,pid in enumerate(visibilities):
        # Add arcs corresponding to the visibilities, i.e.
        # to which cobra can observe which target in which exposure
        for id in visibilities[pid]:
            tc = class_dict[id]
            # these are all the cobras that can reach the target
            cobra_ids = ["{}".format(c) for c in visibilities[pid][id]]

            # bail out if its none
            if cobra_ids == []:
                continue

            if tc.startswith("sci_"):
                tid = dm.SciTarget.getID(id)
                # bail out if we didn't include use this target
                if tid not in g.sciTargets:
                    continue
                # For science targets we need to add edges between cobra pointings and target pointings.
                # So for each pointing, link all cobras that can reach
                # a specific target to that target.

                tvid = dm.TargetPointing.getID(id, pid)
                #tvid = "{}_v{}".format(tid, pointing)
                tv = g.targetPointings[tvid]
                for cid in cobra_ids:
                    cid2 = dm.Cobra.getID(cid)
                    # bail out if we didn't include use this cobra
                    if cid2 not in g.cobras:
                        continue

                    cv = g.cobraPointings[dm.CobraPointing.getID(cid, pid)]
                    a = dm.TargetPointingToCobraPointingArc(tv, cv)
                    a.pointing = pid

                    cx, cy = g.cobras[cid2].fplane_pos[0], g.cobras[cid2].fplane_pos[1]
                    tx, ty = g.sciTargets[tid].fplane_positions[pid][0], g.sciTargets[tid].fplane_positions[pid][1]

                    d = np.abs(np.sqrt((tx-cx)**2. + (ty-cy)**2.))
                    a.cost = costOfD(d)
                    a.d = d
                    e = g.add_arc(a)

            if tc.startswith("cal_") or tc.startswith("sky_"):
                # For calibration targets we need to add edges between cobra pointings and target (not targetPointing).
                    for cid in cobra_ids:
                        cid2 = dm.Cobra.getID(cid)
                        tid = dm.CalTarget.getID(id, pid)

                        # bail out if we didn't include use this target
                        if tid not in g.calTargets:
                            continue

                        t = g.calTargets[tid]

                        # bail out if we didn't include use this cobra
                        if cid2 not in g.cobras:
                            continue

                        cv = g.cobraPointings[dm.CobraPointing.getID(cid, pid)]
                        a = dm.TargetPointingToCobraPointingArc(t, cv)

                        cx, cy = g.cobras[cid2].fplane_pos[0], g.cobras[cid2].fplane_pos[1]
                        tx, ty = g.calTargets[tid].fplane_positions[pid][0], g.calTargets[tid].fplane_positions[pid][1]

                        d = np.abs( np.sqrt((tx-cx)**2. + (ty-cy)**2.))
                        a.cost = costOfD(d)
                        a.d = d

                        g.add_arc(a)

    return g


def compute_collision_pairs(pointings, target_fplane_pos):
    """
    Compute which pairs of targets would lead to enpoint collision
    if assigned simultanously for observation.
    """

    # find collision pairs
    fiber_collision_radius = 1.

    # I am sure this code cound be massively optimized by subdeviding the focal plane
    # also we probably don't need to do this on a per-pointing basis.
    collision_pairs = OrderedDict()


    for pid,(pointing_RA,pointing_DEC) in pointings.items():
        
        start_time = time.time()
        
        targets = target_fplane_pos[pid]

        txx = np.array( [t[0] for tid, t in targets.items()] )
        tyy = np.array( [t[1] for tid, t in targets.items()] )
        ID  = [tid for tid in targets] 
        
        N = len(ID)
        points = list( zip(txx,tyy) )
        Y = cdist( points[:N], points[:N] )

        # any target separation that is smaller than 2 x the collision radius will be flagged a s collision
        cc = Y <= (fiber_collision_radius*2.) 
        ncoll = int( (np.sum(cc.flatten()) - N)/2. )

        print ("Pointing {}: Found  {:d} collision pairs.".format( pid, ncoll  ))

        # identify collision pairs
        collision_pairs[pid] = []
        # array of indices
        ii = np.arange(N)
        for i in range(cc.shape[0]):
            x1,y1 =  txx[i], tyy[i]
            # only iterate over the indeces that are colliding and the upper diagonal in the collision matrix
            jj = ii[ cc[i,:] * ii > i ] 
            for j in jj: 
                if cc[i,j]:
                    x2,y2 =  txx[j], tyy[j]
                    collision_pairs[pid].append([(ID[i],x1,y1),(ID[j],x2,y2)])

        time_to_finish = time.time() - start_time
        print(" Time to completion: {:.2f} s".format(time_to_finish))
        
    return collision_pairs
    

def compute_collision_flow_pairs(g, collision_pairs):
    """
     Identifiy which flow variables correspond to which collision pairs.

     This is not a nice piece of code, mostly because we need to treat science and
     calibrations targets differently
     What we do:
     Loop over all collision pairs (science - science, science - cal, cal - cal)
      then for each pointing
      look if they are actually part of the graph (in case we are dealing with a 
      subregion of the focal plane only we might ignore them)
       identify the input flow arc (ther can be only one) for each of the two targets in the pair
       add the flow pairs to a list
     loop over all flow pairs and add a constraint equation
    """ 
    flow_pairs = []

    for pid in collision_pairs:
        for cp in collision_pairs[pid]:

            for pointing in g.pointings:
                tid1 = cp[0][0]
                tid2 = cp[1][0]

                tvid1 = "T_{}_v{}".format(cp[0][0],pointing)
                tvid2 = "T_{}_v{}".format(cp[1][0],pointing)

                # science targets have targetPointing nodes
                # calibrations targets do not (there is a doublicate for each pointings)
                if tvid1 in g.calTargets:
                    f1id = g.calTargets[tvid1].inarcs[0].id
                elif tvid1 in g.targetPointings:
                    f1id = g.targetPointings[tvid1].inarcs[0].id
                else:
                    continue # this target is not part of the problem, probably did not survive RMAX cut

                if tvid2 in g.calTargets:
                    f2id = g.calTargets[tvid2].inarcs[0].id
                elif tvid2 in g.targetPointings:
                    f2id = g.targetPointings[tvid2].inarcs[0].id
                else:
                    continue # this target is not part of the problem, probably did not survive RMAX cut


                flow_pairs.append([f1id, f2id])
    return flow_pairs




def computeStats(g):
    stats = OrderedDict()

    NSciObs = 0
    NSciComplete = 0
    NCalObs = np.nan
    NCalComplete = np.nan
    NPOINTINGS = len(g.pointings)

    for t in g.sciTargets.values():
        NSciObs += sum([a.flow for a in t.inarcs])
        NSciComplete += int(sum([a.flow for a in t.outarcs]) == t.gain)
            
    Noverflow = 0
    for tcid, tc in g.sciTargetClasses.items():
        aid = '{}=SINK'.format(tcid)
        Noverflow += int( g.arcs[aid].flow )

    Ncobras_used = 0
    Ncobras_fully_used = 0
    for c in g.cobras.values():
        v = sum([a.flow for a in c.inarcs])
        Ncobras_used += int(v > 0)
        Ncobras_fully_used += int(v == NPOINTINGS)
        
    stats['NSciObs'] = NSciObs
    #stats['NCalObs'] = NCalObs
    stats['NSciComplete'] = NSciComplete
    #stats['NCalComplete'] = NCalComplete
    stats['Noverflow'] = Noverflow
    stats['Ncobras_used'] = Ncobras_used
    stats['Ncobras_fully_used'] =  Ncobras_fully_used
    

    compl = {}
    for stc in g.sciTargetClasses.values():
        _NObs = 0
        _NComplete = 0
        for t in stc.targets.values(): 
            _NObs += int( sum([a.flow for a in t.inarcs] ) )
            _NComplete += int(sum([ a.flow for a in t.outarcs]) == t.gain)
            
        compl[stc.ID] = {'total' : len(stc.targets), 'observed' : _NObs, 'completed' : _NComplete}


    stats['completion'] = compl

    return stats


    
def computeCompletion(g, class_dict, outfilename):
    """
    Compute how many science targets, calibration stars and sky positions
    were observed in each poitning. Compute also 
    cumulative number for the science targets
    and write results to "outfilename".
    """
    print("Computing target completion ...")


    filterArcs = lambda a : type(a) == dm.TargetToTargetPointingArc

    # consistency checks ....
    noutarcs = []
    for t in g.sciTargets.values():
        aa = list( filter( filterArcs, t.outarcs) )
        noutarcs.append(len(aa))
        npointings = np.unique( noutarcs )[0]

    if not len( np.unique( noutarcs ) ) == 1:
        print("Error the number of outarcs from each science target not is not equal. Found ", np.unique( noutarcs ))

    if not npointings == len(g.pointings):
        print("Error the number of outarcs from each science target not match len(g.pointings). Found ", np.unique( noutarcs ))

    # now compute per class and per pointing completion analysis
    tclass_completion  = OrderedDict()
    tclasses = np.unique( [t.target_class for t in g.sciTargets.values()] )

    # initilize structure to hold results
    for i, pid in enumerate(g.pointings): 
        tclass_completion[pid] = OrderedDict()
        for tclass in np.unique(tclasses):
            tclass_completion[pid][tclass] = 0

    # count cumulative number of target in each class 
    # that got completed
    for i, pid in enumerate(g.pointings): 
        print("   Pointing: ", pid)

        for t in g.sciTargets.values():
            tclass = t.target_class
            aa = list( filter( filterArcs, t.outarcs) )
            nobs = np.round( np.sum( [a.flow for a in aa[:i+1] ] ) )
            if nobs >= t.gain:
                tclass_completion[pid][tclass] += 1
       

    # now calculate number of observed calibration objects in each pointing
    caltclasses = np.unique( [t.target_class for t in g.calTargets.values()] )
    caltclass_obs  = OrderedDict()
    for pid in g.pointings:
        ptargets = list( filter( lambda t : t.pointing == pid,   g.calTargets.values() ) )
        caltclass_obs[pid] = OrderedDict()
        for ctc in caltclasses:
            #lambda t : class_dict[t.id.split("_")[1] == ctc
            ctctargets = list( filter(  lambda t : t.target_class == ctc ,   ptargets ) )
            nobs = int( sum([t.inarcs[0].flow for t in ctctargets]) )
            caltclass_obs[pid][ctc] = nobs

    # Assemble results in a convenient table
    names = ['N', 'pid'] + [ t for t in np.unique(tclasses)] + [t for t in np.unique(caltclasses)] 
    dtype = ['i4', 'S8'] + ['i4'] * len(np.unique(tclasses)) + ['i4'] * len(np.unique(caltclasses))
    t = Table(names=names, dtype=dtype)

    for i, pid in enumerate( g.pointings ):
        row = [i, pid]

        for tclass in np.unique(tclasses):
            row += [ tclass_completion[pid][tclass] ]
        for tclass in np.unique(caltclasses):
            row += [ caltclass_obs[pid][tclass] ]

        t.add_row(row)

    
    t.write(outfilename, overwrite=True, format="ascii.commented_header")
    print("Done, write {}.".format(outfilename) )
    
    return t
