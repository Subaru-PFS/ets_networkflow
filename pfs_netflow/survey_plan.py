from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import range
import numpy as np
from collections import OrderedDict
import pulp

from . import datamodel as dm


def buildSurveyPlan(cobras, targets, nreqv_dict, visibilities, class_dict,
                    cost_dict, supply_dict, RMAX, CENTER=(0.,0.), COBRAS=[]):
    """
    Builds a graph which represents a survey plan.
    

    Args:
        cobras (OrderedDict): A dictionary of cobra x,y in the focal plane. The key is the cobra ID.
        targets (OrderedDict): A dictionary of target x,y in the focal plane. The key is the target ID.
        nreqv_dict (OrderedDict): Number of required visits per target. The key is the target ID, elements are int.
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
    g.visits = list(visibilities)

    # Add global sink node
    T = dm.Sink()
    g.add_node(T)

    # Add nodes for the cobras and the cobra visits and arcs between them
    for cid in cobras:
        x, y = cobras[cid]
        if (x-CENTER[0])**2 + (y-CENTER[1])**2 > RMAX**2:
            continue
        if COBRAS != [] and cid not in COBRAS:
            continue
        c = dm.Cobra(cid, fplane_pos=(x,y))
        g.add_node(c)
        # replicate node as many times as there are visits
        for pid in visibilities:
            cv = dm.CobraVisit(cid=cid, cobra=c, visit=pid)
            g.add_node(cv)
            g.add_arc(dm.CobraVisitToCobraArc(cv, c))

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
            for visit in g.visits:
                targetClass = dm.CalTargetClass(tc, visit)
                targetClass.cost = cost_dict[tc]
                targetClass.supply = supply_dict[tc]
                g.add_node(targetClass)

                # Add costly overflow arc to the target class
                g.add_arc(dm.OverflowArc(targetClass.cost, targetClass, T))
        else:
            print("Error unknown target class {}".format(tc))


    # Add nodes for the targets and the target visits and arcs between them
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
            t.gain = nrv
            g.add_node(t)
            
            # Add as many TargetVisit nodes for this target as there are visits
            for pid in visibilities:
                tv = dm.TargetVisit(tid, target=t, visit=pid)
                g.add_node(tv)
                ttva = dm.TargetToTargetVisitArc(t, tv)
                # Here we assign the cost for the respective visit
                # increasing the cost for later visits encourages
                # earlier observation.
                ttva.cost = cost_dict["visits"][pid]
                g.add_arc(ttva)

            # Add arc from target class to target
            targetClass = g.sciTargetClasses[dm.SciTargetClass.getID(tc)]
            targetClass.add_target(t)
            g.add_arc(dm.TargetClassToTargetArc(targetClass, t))

            # Add VERY costly overflow arc from targets to sink
            g.add_arc(dm.OverflowArc(targetClass.cost_partial_compl, t, T))

        elif tc.startswith("cal_") or tc.startswith("sky_"):
            for pid in visibilities:
                # Add as many nodes for each calibration target as there are visits
                if tc.startswith("cal_"):
                    t = dm.StarCalTarget(tid, fplane_positions=fplane_positions, visit=pid)
                else:
                    t = dm.SkyCalTarget(tid, fplane_positions=fplane_positions, visit=pid)
                #t = dm.CalTarget(tid, fplane_positions=fplane_positions, visit=pid)
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

    for visit,pid in enumerate(visibilities):
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
                # For science targets we need to add edges between cobra visits and target visits.
                # So for each visit, link all cobras that can reach
                # a specific target to that target.

                tvid = dm.TargetVisit.getID(id, pid)
                #tvid = "{}_v{}".format(tid, visit)
                tv = g.targetVisits[tvid]
                for cid in cobra_ids:
                    cid2 = dm.Cobra.getID(cid)
                    # bail out if we didn't include use this cobra
                    if cid2 not in g.cobras:
                        continue

                    cv = g.cobraVisits[dm.CobraVisit.getID(cid, pid)]
                    a = dm.TargetVisitToCobraVisitArc(tv, cv)
                    a.visit = visit

                    cx, cy = g.cobras[cid2].fplane_pos[0], g.cobras[cid2].fplane_pos[1]
                    tx, ty = g.sciTargets[tid].fplane_positions[pid][0], g.sciTargets[tid].fplane_positions[pid][1]

                    d = np.abs(np.sqrt((tx-cx)**2. + (ty-cy)**2.))
                    a.cost = costOfD(d)
                    a.d = d
                    e = g.add_arc(a)

            if tc.startswith("cal_") or tc.startswith("sky_"):
                # For calibration targets we need to add edges between cobra visits and target (not targetVisit).
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

                        cv = g.cobraVisits[dm.CobraVisit.getID(cid, pid)]
                        a = dm.TargetVisitToCobraVisitArc(t, cv)

                        cx, cy = g.cobras[cid2].fplane_pos[0], g.cobras[cid2].fplane_pos[1]
                        tx, ty = g.calTargets[tid].fplane_positions[pid][0], g.calTargets[tid].fplane_positions[pid][1]

                        d = np.abs( np.sqrt((tx-cx)**2. + (ty-cy)**2.))
                        a.cost = costOfD(d)
                        a.d = d

                        g.add_arc(a)

    return g



def compute_collision_flow_pairs(g, collision_pairs):
    """
     Identifiy which flow variables correspond to which collision pairs.

     This is not a nice piece of code, mostly because we need to treat science and
     calibrations targets differently
     What we do:
     Loop over all collision pairs (science - science, science - cal, cal - cal)
      then for each visit
      look if they are actually part of the graph (in case we are dealing with a 
      subregaion of the focal plane only we might ignore them)
       identify the input flow arc (ther can be only one) for each of the two targets in the pair
       add the flow pairs to a list
     loop over all flow pairs and add a constraint equation
    """ 
    flow_pairs = []

    for pid in collision_pairs:
        for cp in collision_pairs[pid]:

            for visit in g.visits:
                tid1 = cp[0][0]
                tid2 = cp[1][0]

                tvid1 = "T_{}_v{}".format(cp[0][0],visit)
                tvid2 = "T_{}_v{}".format(cp[1][0],visit)

                # science targets have targetVisit nodes
                # calibrations targets do not (there ais a doublicate for each visits)
                if tvid1 in g.calTargets:
                    f1id = g.calTargets[tvid1].inarcs[0].id
                elif tvid1 in g.targetVisits:
                    f1id = g.targetVisits[tvid1].inarcs[0].id
                else:
                    continue # this target is not part of the problem, probably did not survive RMAX cut

                if tvid2 in g.calTargets:
                    f2id = g.calTargets[tvid2].inarcs[0].id
                elif tvid2 in g.targetVisits:
                    f2id = g.targetVisits[tvid2].inarcs[0].id
                else:
                    continue # this target is not part of the problem, probably did not survive RMAX cut


                flow_pairs.append([f1id, f2id])
    return flow_pairs



def computeStats(g, flows, cost):
    stats = OrderedDict()

    NSciObs = 0
    NSciComplete = 0
    NCalObs = np.nan
    NCalComplete = np.nan
    NVISITS = len(g.visits)

    for t in g.sciTargets.values():
        NSciObs += pulp.value(sum([flows[a.id] for a in t.inarcs]))
        NSciComplete += int(sum([pulp.value(flows[a.id]) for a in t.outarcs]) == t.gain)
            
    Noverflow = 0
    for tcid, tc in g.sciTargetClasses.items():
        Noverflow += int( pulp.value(flows['{}=SINK'.format(tcid)]) )

    Ncobras_used = 0
    Ncobras_fully_used = 0
    for c in g.cobras.values():
        v = pulp.value(sum([flows[a.id] for a in c.inarcs]))
        Ncobras_used += int(v > 0)
        Ncobras_fully_used += int(v == NVISITS)
        
    stats['cost'] = pulp.value(cost)
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
            _NObs += int( pulp.value(sum([flows[a.id] for a in t.inarcs])) )
            _NComplete += int(sum([pulp.value(flows[a.id]) for a in t.outarcs]) == t.gain)
            
        compl[stc.ID] = {'total' : len(stc.targets), 'observed' : _NObs, 'completed' : _NComplete}


    stats['completion'] = compl

    return stats



    
def computeTargetCompletion(g, visibilities, class_dict, outfilename):
    """
    Compute how many science targets, calibration stars and sky positions
    were observed in each poitning. Compute also 
    cumulative number for the science targets
    and write results to "outfilename".
    """
    print("Computing target completion ...")
    nsci_observed_total = 0
    NVISITS = len(g.visits)
    with open(outfilename, 'w') as f:
        s  = "# algorithm: {}\n".format("netflow")
        s += "# total number of vistis: {}\n".format(NVISITS)
        s += "{:5s} {:5s} {:5s} {:5s} {:7s} {:7s} {:7s} {:7s}\n"\
            .format("V", "nsci", "ncal", "nsky", "nsci_obs", "ncal_obs", "nsky_obs", "nsci_cum")
        f.write(s)

        for pid in g.visits:
            print("Pointing {}".format(pid))
            nsci_observable = 0
            ncal_observable = 0
            nsky_observable = 0

            for tid,v in visibilities[pid].items():
                if len(v) > 0:
                    if class_dict[tid][:3] == 'cal':
                        ncal_observable += 1
                    elif class_dict[tid][:3] == 'sky':
                        nsky_observable += 1
                    elif class_dict[tid][:3] == 'sci':
                        nsci_observable += 1
                
                
            #print("Pointing {}".format(pid))
            nsci = 0
            ncal = 0
            nsky = 0
            for a in g.arcs.values():
                n1,n2 = a.startnode,a.endnode

                if a.flow > 0.:
                    if type(n2) == dm.CobraVisit and n2.visit == pid:
                        if type(n1) == dm.TargetVisit:
                            nsci += 1
                            nsci_observed_total += 1
                        elif type(n1) == dm.SkyCalTarget:
                            nsky  += 1
                        elif type(n1) == dm.StarCalTarget:
                            ncal  += 1

            s = "{:5s} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d}\n".\
                format(pid, nsci, ncal, nsky, nsci_observable, ncal_observable, nsky_observable, nsci_observed_total)
            f.write(s)

            #print(" Observed {} science targets, {} calibration targets and {} sky positions.".format(nsci, ncal, nsky))
            #print(" Observed {} science targets in total.".format(nsci_observed_total))
    print("Done, write {}.".format(outfilename) )