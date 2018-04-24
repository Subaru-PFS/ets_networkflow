import numpy as np
from numpy import hstack
from numpy import cos, deg2rad, sqrt
from numpy import unique, array
from collections import OrderedDict
from numpy import random

def pp(s):
    print(s)
    return s + "\n"

def combineTargetLists(science_targets, cal_stars, sky_pos, DRMAX_SQ):
    ID         = hstack([science_targets['ID'],cal_stars['ID'],sky_pos['ID']])
    ra         = hstack([science_targets['RA'],cal_stars['RA'],sky_pos['RA']])
    dec        = hstack([science_targets['DEC'],cal_stars['DEC'],sky_pos['DEC']])
    exp_times  = hstack([science_targets['EXP_TIME'],cal_stars['EXP_TIME'],sky_pos['EXP_TIME']])
    priorities = hstack([science_targets['Priority'],cal_stars['Priority'],sky_pos['Priority']])

    # make up target classes from target type and priority
    types = ["sci"]*len(science_targets) + ["cal"]*len(cal_stars) + ["sky"]*len(sky_pos)
    class_dict = {}
    for id,t,p in zip(ID, types, priorities):
        class_dict[id] = '{}_P{:02d}'.format(t, p)
        
    # find field center and extent
    cRA  = np.median( science_targets['RA'] )  
    cDEC =  np.median( science_targets['DEC'] ) 
    RAmin,RAmax   = np.min( science_targets['RA'] ), np.max( science_targets['RA'] )
    DECmin,DECmax = np.min( science_targets['DEC'] ), np.max( science_targets['DEC'] )


    dra = (ra - cRA)*cos(deg2rad(cDEC))
    ddec = dec - cDEC
    ii = (dra**2. + ddec **2.) <= DRMAX_SQ

    ID = ID[ii].tolist()
    ra = ra[ii].tolist()
    dec = dec[ii].tolist()
    exp_times = exp_times[ii].tolist()
    priorities = priorities[ii].tolist()
    types = np.array(types)[ii].tolist()

    c = np.array([ class_dict[t][:3] for t in ID ])
    ii_sci = c == 'sci'
    ii_cal = c == 'cal'
    ii_sky = c == 'sky'

    return ID, ra, dec, exp_times, priorities, types, class_dict, ii_sci, ii_cal, ii_sky

def randomPriorities(ii_sci, high=1, low=4):
    """
    Create random priorities.
    """
    from numpy import random
    from numpy import array

    random.seed(42)

    N =  sum(ii_sci) 

    newpri = array( random.uniform(high, low, size=N) , dtype=int)
    priorities = np.array(priorities)

    priorities[ii_sci] = newpri
    priorities = priorities.tolist()
    
    return priorities

def randomExposures(ID, NPOINTINGS=10):
    """
    Create a uniform random number of required visits per target
    from 1 to NPOINTINGS.
    """
    nv = np.floor( random.uniform(NPOINTINGS+1, size=len(ID)) ) 

    # compute number of required visists from exposure times
    # and block length
    nreqv_dict = {}
    for id,t,nrv in zip(ID, types, nv):
        nreqv_dict[id] = int(nrv)

    print( "Required repointings", unique( [v for v in nreqv_dict.values()] ))
    
    return nreqv_dict 


def invert_vis(visibilities):
    """
    Invert visibility map, i.e. for each cobra list the observable targets
    """
    ivisibilities = OrderedDict()
    for pid in visibilities:
        cnt = 0
        ivisibilities[pid] = OrderedDict()
        for v in visibilities[pid]:
            t = v
            cc = visibilities[pid][v]

            for c in cc:
                if c in ivisibilities[pid]:
                    ivisibilities[pid][c].append(v)
                else:
                    ivisibilities[pid][c] = [v]

    return ivisibilities


def computeObservabityStats(visibilities, class_dict, ii_sci, ii_cal, ii_sky):
    ivisibilities = invert_vis(visibilities)           

    for pid in visibilities:
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
                #break


        print(" {} targets positions in total.".format(sum(ii_sci) ))
        print(" {} cal. targets in total.".format(sum(ii_cal) ))
        print(" {} sky positions in total.".format(sum(ii_sky) ))
        print(" {:6d} cobras have at least one target in reach.".format(len(ivisibilities[pid])))
        print(" {:6d} science targets\n {:6d} calibration targets,\n {:6d} sky positions \nare in reach of at least one cobra."\
              .format(nsci_observable, ncal_observable, nsky_observable))  
        print("")
     
def convert_to_newstyle_ets(ID, ets_cobras, ets_target_pos, visibility_maps):
    # obtain cobra centers in old ETS dictionary style
    cobras = OrderedDict()
    for i,c in enumerate(ets_cobras):
            x,y = np.real( ets_cobras[i][0] ), np.imag( ets_cobras[i][0] )
            cobras["{:d}".format(i)] = [x,y]

    # obtain dot centers in dictionary style
    dots = OrderedDict()
    for i,c in enumerate(ets_cobras):
            x,y = np.real( ets_cobras[i][-2] ), np.imag( ets_cobras[i][-2] )
            r = ets_cobras[i][-1]
            dots["{:d}".format(i)] = [x,y,r]

    # obtain target positions in old ETS dictionary style
    target_fplane_pos = OrderedDict()
    for pid in ets_target_pos:
        tt = OrderedDict()
        for j,c in enumerate(ets_target_pos[pid]):
                x,y = np.real(c),np.imag(c)
                tt[ID[j]] = [float(x),float(y)]
        target_fplane_pos[pid] = tt

    # obtain visibilities in old ETS dictionary style
    visibilities = OrderedDict()
    for pid in visibility_maps:
        vv = OrderedDict()   
        for v in visibility_maps[pid]:
            t = ID[v]
            cc = ["{:d}".format(c) for c in visibility_maps[pid][v]]
            vv[t] = cc
        visibilities[pid] = vv 

    # for per-pointing list of targets, build
    # single list of targets with multiple entries for the focal plane positions
    targets = OrderedDict()
    for pid in target_fplane_pos:
        tt = target_fplane_pos[pid] 
        for tid in tt:
            x,y = tt[tid]
            if not tid in targets:
                 targets[tid] = OrderedDict()
            targets[tid][pid] = (x,y)
        
    return cobras, dots, target_fplane_pos, visibilities, targets