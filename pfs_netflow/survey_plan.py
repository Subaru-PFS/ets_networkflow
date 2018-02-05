from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import range
import numpy as np

from . import datamodel as dm


def buildSurveyPlan(cobras, targets, nreqvisits, visibilities, class_dict,
                    cost_dict, supply_dict, NVISITS, RMAX, CENTER, COBRAS=[]):
    print("buildSurveyPlan")

    # Check if a visit specific cost was given as part of the cost function.
    # Set all visits to zero (extra) cost otherwise.
    if 'visits' not in cost_dict:
        # no specific cost for the visits is give, make them all zero cost
        cost_dict['visits'] = [0.] * NVISITS
    else:
        assert len(cost_dict['visits']) >= NVISITS, \
        "The length of the array in (cost_dict['visits']) (= {}) is less than NVISITS (= {})."\
            .format(len(cost_dict['visits']), NVISITS)

    # Check if a move distance dependant cost function was given. If not set to zero.
    if 'cobra_move' not in cost_dict:
        # no specific cost function for the cobra move distance is given. Set to zero cost.
        costOfD = lambda x : 0.
    # MR FIXME: "else" missing?
    costOfD = cost_dict['cobra_move']

    # Check if visibility map was passed as list. If so the first entry will be used for the
    # first visit, the second for the second and so forth.
    if isinstance(visibilities, list):
        assert len(visibilities) >= NVISITS,\
            "The length of the visibilities list (= {}) is less than NVISITS (= {})."\
            .format(len(visibilities), NVISITS)
    else:
        visibilities = [visibilities] * NVISITS

    TARGETS = []
    # generate a directed graph
    g = dm.SurveyPlan()
    g.visits = list(range(NVISITS))

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
        c = dm.Cobra(cid, x=x, y=y)
        g.add_node(c)
        # replicate node as many times as there are visits
        for visit in range(NVISITS):
            cv = dm.CobraVisit(cid=cid, cobra=c, visit=visit)
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
    for tid, nrv in zip(targets, nreqvisits):
        tc = class_dict[tid]
        x, y = targets[tid]
        if (x-CENTER[0])**2 + (y-CENTER[1])**2 > RMAX**2:
            continue
        if TARGETS != [] and tid not in TARGETS:
            continue
        if tc.startswith("sci_"):
            t = dm.SciTarget(tid, x=x, y=y)
            t.gain = nrv
            g.add_node(t)

            # Add as many TargetVisit nodes for this target as there are visits
            for visit in g.visits:
                tv = dm.TargetVisit(tid, target=t, visit=visit)
                g.add_node(tv)
                ttva = dm.TargetToTargetVisitArc(t, tv)
                # Here we assign the cost for the respective visit
                # increasing the cost for later visits encourages
                # earlier observation.
                ttva.cost = cost_dict["visits"][visit]
                g.add_arc(ttva)

            # Add arc from target class to target
            targetClass = g.sciTargetClasses[dm.SciTargetClass.getID(tc)]
            targetClass.add_target(t)
            g.add_arc(dm.TargetClassToTargetArc(targetClass, t))

            # Add VERY costly overflow arc from targets to sink
            g.add_arc(dm.OverflowArc(targetClass.cost_partial_compl, t, T))

        elif tc.startswith("cal_") or tc.startswith("sky_"):
            for visit in g.visits:
                # Add as many nodes for each calibration target as there are visits
                t = dm.CalTarget(tid, x=x, y=y, visit=visit)
                t.gain = 1
                g.add_node(t)

                # Add arc from target class to target
                targetClass = g.calTargetClasses[dm.CalTargetClass.getID(tc, visit)]
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

    for visit in g.visits:
        # Add arcs corresponding to the visibilies, i.e.
        # to which cobra can observe which target in which exposure
        for id in visibilities[visit]:
            tc = class_dict[id]
            # these are all the cobras that can reach the target
            cobra_ids = ["{}".format(c) for c in visibilities[visit][id]]

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

                tvid = dm.TargetVisit.getID(id, visit)
                #tvid = "{}_v{}".format(tid, visit)
                tv = g.targetVisits[tvid]
                for cid in cobra_ids:
                    cid2 = dm.Cobra.getID(cid)
                    # bail out if we didn't include use this cobra
                    if cid2 not in g.cobras:
                        continue

                    cv = g.cobraVisits[dm.CobraVisit.getID(cid, visit)]
                    a = dm.TargetVisitToCobraVisitArc(tv, cv)
                    a.visit = visit

                    cx, cy = g.cobras[cid2].x, g.cobras[cid2].y
                    tx, ty = g.sciTargets[tid].x, g.sciTargets[tid].y

                    d = np.abs(np.sqrt((tx-cx)**2. + (ty-cy)**2.))
                    a.cost = costOfD(d)
                    a.d = d
                    e = g.add_arc(a)

            if tc.startswith("cal_") or tc.startswith("sky_"):
                # For calibration targets we need to add edges between cobra visits and target (not visit).
                    for cid in cobra_ids:
                        cid2 = dm.Cobra.getID(cid)
                        tid = dm.CalTarget.getID(id, visit)

                        # bail out if we didn't include use this target
                        if tid not in g.calTargets:
                            continue

                        t = g.calTargets[tid]

                        # bail out if we didn't include use this cobra
                        if cid2 not in g.cobras:
                            continue

                        cv = g.cobraVisits[dm.CobraVisit.getID(cid, visit)]
                        a = dm.TargetVisitToCobraVisitArc(t, cv)

                        cx, cy = g.cobras[cid2].x, g.cobras[cid2].y
                        tx, ty = g.calTargets[tid].x, g.calTargets[tid].y

                        d = np.abs( np.sqrt((tx-cx)**2. + (ty-cy)**2.))
                        a.cost = costOfD(d)
                        a.d = d

                        g.add_arc(a)

    return g
