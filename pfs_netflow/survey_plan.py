from numpy import unique, sqrt, abs
from numpy import inf

import datamodel as dm

def setflows(g,flows):
    for a in g.arcs.itervalues():
        k = '{}={}'.format(a.startnode.id,a.endnode.id)
        if flows.has_key(k):
            a.flow = value(flows[k])


def buildSurveyPlan(cobras, targets, nreqvisits, visibilities, class_dict, \
                    cost_dict, supply_dict, NVISITS, RMAX, CENTER, COBRAS = []):
    print("buildSurveyPlan " )
    
    if not 'visits' in cost_dict:
        # no specific cost for the visits is give, make them all zero cost
        cost_dict['visits'] = [0.] * NVISITS
    else:
        assert len(cost_dict['visits']) >= NVISITS, \
        "The length of the array in (cost_dict['visits']) (= {}) is less than NVISITS (= {})."\
            .format(len(cost_dict['visits']),NVISITS) 
    
    if not 'cobra_move' in cost_dict:
        # no specific cost function for the cobra move distance is given. Set to zero cost.
        costOfD = lambda x : 0.
    costOfD = cost_dict['cobra_move']
    
    TARGETS = []

    # generate a directed graph 
    g = dm.SurveyPlan()

    g.visits = range(NVISITS)

    # Add global sink node
    T = dm.Sink("SINK")
    g.add_node(T)

    # Add nodes for the cobras and the cobra visists
    #  and arcs between them
    for cid in cobras:  
        x,y = cobras[cid] 
        if (x-CENTER[0])**2 + (y-CENTER[1])**2 > RMAX**2:
            continue
        if COBRAS != [] and not cid in COBRAS:
            continue
        c = dm.Cobra(id="C_{}".format(cid),x=x,y=y)
        g.add_node(c)
        # replicate node as many times as there are visits
        for visit in range(NVISITS):
            cvid = "C_{}_v{}".format(cid,visit)
            cv = dm.CobraVisit(id=cvid,cobra=c,visit=visit)
            g.add_node(cv)
            cvca = dm.CobraVisitToCobraArc(cv, c)
            g.add_arc(cvca)
    
    # Add nodes for target classes
    for tc in unique( class_dict.values() ):
        # add target-class super node
        if tc.startswith("sci_"):
            targetClass = dm.SciTargetClass("TClass_{}".format(tc))
            # give target classes all the same cost fro now.
            # we will introduce a machanism to assign a
            # class specific cost later.
            targetClass.cost = cost_dict[tc][0] 
            targetClass.cost_partial_compl = cost_dict[tc][1]
            targetClass.supply = supply_dict[tc]
            g.add_node(targetClass)
        elif tc.startswith("cal_") or tc.startswith("sky_"):
            for visit in g.visits:
                targetClass = dm.CalTargetClass("TClass_{}_v{}".format(tc,visit))
                # give calibration targets a larger cost. 
                # Again, we will introduce a machanism to assign a
                # class specific cost later.
                targetClass.cost = cost_dict[tc]
                targetClass.supply = supply_dict[tc]
                g.add_node(targetClass)
        else:
            print("Error unknown target class {}".format(tc))
                
        
        # Add costly overflow arc to the target class
        targetClassOverFlowArc = dm.OverflowArc(targetClass.cost, targetClass, T)
        g.add_arc( targetClassOverFlowArc )
    
    # Add nodes for the targets and the target visists
    #  and arcs between them
    for tid,nrv in zip(targets, nreqvisits):
        tc = class_dict[tid]
        x,y = targets[tid]
        if (x-CENTER[0])**2 + (y-CENTER[1])**2 > RMAX**2:
            continue
        if TARGETS != [] and not tid in TARGETS:
            continue
        if tc.startswith("sci_"):
            t = dm.SciTarget(id="T_{}".format(tid),x=x,y=y)
            t.gain = nrv
            g.add_node( t )
            
            # Add as many TargetVisit nodes for this target as there are visits
            for visit in g.visits:
                tvid = "T_{}_v{}".format(tid,visit)
                tv = dm.TargetVisit(id=tvid,target=t,visit=visit)
                g.add_node( tv )
                ttva = dm.TargetToTargetVisitArc(t, tv)
                # Here we assign the cost for the respective visit
                # increasing the cost for later visits encourages
                # earlier observation.
                ttva.cost = cost_dict["visits"][visit]
                e = g.add_arc(ttva)

            # Add arc from target class to target
            targetClass = g.sciTargetClasses["TClass_{}".format(tc)]
            targetClass.add_target(t)
            g.add_arc( dm.TargetClassToTargetArc(targetClass, t) )
            
            # Add VERY costly overflow arc from targets to sink
            g.add_arc( dm.OverflowArc(targetClass.cost_partial_compl, t, T) )

        
        elif tc.startswith("cal_") or tc.startswith("sky_"):
            for visit in g.visits:
                # Add as many nodes for each calibration target as there are visits
                t = dm.CalTarget(id="T_{}_v{}".format(tid,visit),x=x,y=y,visit=visit)
                t.gain = 1
                g.add_node( t )
                
                # Add arc from target class to target
                targetClass = g.calTargetClasses["TClass_{}_v{}".format(tc,visit)]
                targetClass.add_target(t)
                g.add_arc( dm.TargetClassToTargetArc(targetClass, t) )  

            
    # Calibration targets only need to be observed once.
    # Remember we replicate calibration targets such that there is one per 
    # observation.
    for tid,t in g.calTargets.iteritems():
         t.gain = 1.

    # Set supply for each science target class. 
    for tcid,targetClass in g.sciTargetClasses.iteritems():
        if targetClass.supply == inf:
            # Rather than using infinite supply (which woudl result in infinite cost)
            # set supply for the targetClasses to number of targets. 
            # I.e. we want to get them all observed (if possible).
            targetClass.supply = len(targetClass.targets)
        
    # Set supply for each calibration target class. 
    for tcid,targetClass in g.calTargetClasses.iteritems():
        if targetClass.supply == inf:
            targetClass.supply = len(targetClass.targets)
        
        
    # Add arcs corresponding to the visibilies, i.e. 
    # to which cobra can observe which target in which exposure
    for id in visibilities:
        tc = class_dict[id]
        # these are all the cobras that can reach the target 
        cobra_ids = ["C_{}".format(c) for c in visibilities[id]]
        
        # bail out if its none
        if cobra_ids == []:
            continue
        
        if tc.startswith("sci_"):
            tid="T_{}".format(id)
            # bail out if we didn't include use this target
            if not g.sciTargets.has_key(tid):
                continue
            # For science targets we need to add edges between cobra visits and target visits.
            # So for each visit, link all cobras that can reach
            # a specific target to that target.
            for visit in g.visits:
                #print "x:", visit, cobra_ids
                tvid = "{}_v{}".format(tid,visit)
                tv = g.targetVisits[tvid]
                for cid in cobra_ids:
                    #print "cid:", cid, g.cobras.has_key(cid)
                    # bail out if we didn't include use this cobra
                    if not g.cobras.has_key(cid):
                        continue
                    cvid = "{}_v{}".format(cid,visit)
                    cv = g.cobraVisits[cvid]
                    a = dm.TargetVisitToCobraVisitArc(tv,cv)
                    
                    cx,cy = g.cobras[cid].x, g.cobras[cid].y
                    tx,ty = g.sciTargets[tid].x, g.sciTargets[tid].y
               
                    d = abs( sqrt( (tx-cx)**2. + (ty-cy)**2. ) )
                    a.cost = costOfD(d)
                    a.d = d
                    e = g.add_arc(a)


        
        if tc.startswith("cal_") or tc.startswith("sky_"):
            # For calibration targets we need to add edges between cobra visits and target (not visit).
            for visit in g.visits:
                for cid in cobra_ids:
                    tid = "T_{}_v{}".format(id,visit)
                    
                    # bail out if we didn't include use this target
                    if not g.calTargets.has_key(tid):
                        continue
                    
                    t = g.calTargets[tid]
                    
                    
                    # bail out if we didn't include use this cobra
                    if not g.cobras.has_key(cid):
                        continue
                    cvid = "{}_v{}".format(cid,visit)
                    
                    
                    cv = g.cobraVisits[cvid]
                    a = dm.TargetVisitToCobraVisitArc(t,cv)
                    
                    cx,cy = g.cobras[cid].x, g.cobras[cid].y
                    tx,ty = g.calTargets[tid].x, g.calTargets[tid].y
                    
                    d = abs( sqrt( (tx-cx)**2. + (ty-cy)**2. ) )
                    a.cost = costOfD(d)
                    a.d = d
                    
                    e = g.add_arc(a)
                
    return g