import matplotlib.pyplot as plt

import matplotlib.lines as mlines
from matplotlib import patches as mpatches
from matplotlib import collections 
import numpy as np

import datamodel as dm

def plotSurveyPlan(g, name=""):
    
    def nlabel(ax, xy, text):
        y = xy[1] - 0.  # shift y-value for label so that it's below the artist
        x = xy[0] + 0.015  # shift x-value for label so that it's next to the artist
        ax.text(x, y, text, family='sans-serif', ha = 'left', va='center', fontsize=6)

    def plotNodes(ax, nodes, x, color, label, allnodes):
        # targets
        N = len(nodes)
        grid = [[x, 1. - (float(i)+1.)/float(N+1)] for i in range(N)]
        patches = []
        for i,(id,n) in enumerate(nodes.iteritems()):
            n.px,n.py = grid[i]
            circle = mpatches.Circle(grid[i],.01, facecolor=color, edgecolor='grey', lw=1.)
            patches.append(circle)
            nlabel(ax, grid[i], id)
            
        collection = collections.PatchCollection(patches, match_original=True)
        
        ax.add_collection(collection)
        allnodes += nodes.keys()
        plt.text(x,0.,label, ha = 'center', va='top')
     
    def plotArcs(ax, g, allnodes):
        for a in g.arcs.itervalues():
            n1,n2 = a.startnode,a.endnode
            if n1.id in allnodes and n2.id in allnodes:
                # add a line
                x, y = np.array([[n1.px,n2.px], [n1.py,n2.py]])
                lw = a.flow * 2 + 1.
                line = mlines.Line2D(x , y , lw=lw, alpha=.3,zorder=0)
                ax.add_line(line)

    
    cmap=plt.cm.hsv
              
    if True:
        fig = plt.figure(figsize=[15,15])
        # science targets
        ax = plt.subplot(111) # note we must use plt.subplots, not plt.subplot
        allnodes = [] # keep track of nodes actually contained in the plot
        plotNodes(ax, nodes=g.cobras, x=0.9, label="cobras", color="#19967d", allnodes=allnodes)
        plotNodes(ax, nodes=g.cobraVisits, x=0.7, label="cobras\nvisits", color="#99d5ca", allnodes=allnodes)
        plotNodes(ax, nodes=g.targetVisits, x=0.5, label="target\nvisits", color="#ffbbb1", allnodes=allnodes)
        plotNodes(ax, nodes=g.sciTargets, x=0.3, label="science\ntargets", color="#ff8f80", allnodes=allnodes)
        plotNodes(ax, nodes=g.sciTargetClasses, x=0.1, label="science\ntarget\nclasses", color="#c92d39", allnodes=allnodes)
        plotArcs(ax, g, allnodes)
        ax.set_ylim([-.0,1.0])
        ax.set_xlim([-.0,1.0])
        plt.axis('off')
        plt.savefig("{}_sci.pdf".format(name))
            
    if True:
        fig = plt.figure(figsize=[15,15])
        # calibration targets
        ax = plt.subplot(111) # note we must use plt.subplots, not plt.subplot
        allnodes = []
        plotNodes(ax, nodes=g.cobras, x=0.9, label="cobras", color="#19967d", allnodes=allnodes)
        plotNodes(ax, nodes=g.cobraVisits, x=0.7, label="cobras\nvisits",color="#99d5ca", allnodes=allnodes)
        plotNodes(ax, nodes=g.calTargets, x=0.4, label="calib\ntargets",color="#ffeca9", allnodes=allnodes)
        plotNodes(ax, nodes=g.calTargetClasses, x=0.2, label="calib.\ntarget\nclasses",color="#ef8d22", allnodes=allnodes)
        plotArcs(ax, g, allnodes)
        ax.set_ylim([-.0,1.0])
        ax.set_xlim([-.0,1.0])
        plt.axis('off')
        
        plt.savefig("{}_cal.pdf".format(name))
        
        
        
import matplotlib.pyplot as plt

import matplotlib.lines as mlines
from matplotlib import patches as mpatches
from matplotlib import collections 
import numpy as np

def plotFocalPlane(g, visit, summary="", XC=0., YC=0., W=400., name=""):
    
    def nlabel(ax, xy, text):
        y = xy[1] - 0.  # shift y-value for label so that it's below the artist
        x = xy[0] + 0.015  # shift x-value for label so that it's next to the artist
        ax.text(x, y, text, family='sans-serif', ha = 'left', va='center', fontsize=6)

    def plotNodes(ax, nodes, color, label, allnodes):
        # targets
        N = len(nodes)
        patches = []
        for i,(id,n) in enumerate(nodes.iteritems()):
            n.px,n.py = n.x,n.y
            circle = mpatches.Circle((n.x,n.y),.75, facecolor=color, edgecolor='grey', lw=1.)
            patches.append(circle)
            
        collection = collections.PatchCollection(patches, match_original=True)
        ax.add_collection(collection)
        allnodes += nodes.keys()
     
    def plotArcs(ax, g, allnodes, visit):
        for a in g.arcs.itervalues():
            n1,n2 = a.startnode,a.endnode
            if a.flow > 0.:
                alpha = 1.
            else:
                alpha = .2
                        
            if type(n2) == dm.CobraVisit and n2.visit == visit and type(n1) == dm.TargetVisit:
                    x, y = np.array([[n1.target.x,n2.cobra.x], [n1.target.y,n2.cobra.y]])

                    line = mlines.Line2D(x , y, lw=1, alpha=alpha, zorder=10)
                    ax.add_line(line)
            if type(n2) == dm.CobraVisit and n2.visit == visit and type(n1) == dm.CalTarget:
                    x, y = np.array([[n1.x,n2.cobra.x], [n1.y,n2.cobra.y]])
                    
                    line = mlines.Line2D(x , y, lw=1, alpha=alpha, zorder=10)
                    ax.add_line(line)                    

    
    cmap=plt.cm.hsv
              

    fig = plt.figure(figsize=[15,15])
    # science targets
    ax = plt.subplot(111) # note we must use plt.subplots, not plt.subplot
    allnodes = [] # keep track of nodes actually contained in the plot
    plotNodes(ax, nodes=g.cobras, label="cobras", color="#19967d", allnodes=allnodes)
    plotNodes(ax, nodes=g.sciTargets, label="science\ntargets", color="#ff8f80", allnodes=allnodes)
    #plotNodes(ax, nodes=g.calTargets, label="calib\ntargets",color="#ffeca9", allnodes=allnodes)
    plotArcs(ax, g, allnodes, visit)

    ax.set_ylim([YC-W/2.,YC+W/2.])
    ax.set_xlim([XC-W/2.,XC+W/2.])
    #plt.axis('off')
    plt.text(-0.1,-0.1, summary, ha='left', va='bottom', transform=ax.transAxes, fontsize=8)
    plt.text(0.5,1.0, "visit {}/{}".format(visit+1,len(g.visits)), ha='center', va='top', transform=ax.transAxes)
    plt.savefig("{}_visit{}_fp.pdf".format(name,visit))


#for a in g.arcs.itervalues():
#    k = '{}={}'.format(a.startnode.id,a.endnode.id)
#    if flows.has_key(k):
#        a.flow = value(flows[k])

#for v in range(NVISITS):
#    plotFocalPlane(g, v, summary, XC=-5., YC=90, W=110)

def printSurveyPlan(g):
    print("Visits: {}".format(g.visits))
    for cid in g.cobras:
        print("Cobra {}".format(cid))
    for cvid in g.cobraVisits:
        print("CobraVisit {}".format(cvid))
    for tid in g.sciTargetClasses:
        print("SciTargetClass {}".format(tid))
    for tid in g.calTargetClasses:
        print("CalTargetClass {}".format(tid))
    for tid in g.calTargets:
        print("CalTarget {}".format(tid))
    for tid in g.sciTargets:
        print("SciTarget {}".format(tid))
    for tvid in g.targetVisits:
        print("TargetVisit {}".format(tvid))

        
def printCalSurveyPlan(g):
    allnodes = []
    
    print("Visits: {}".format(g.visits))
    for tid in g.calTargetClasses:
        print("CalTargetClass {}".format(tid))
    allnodes += g.calTargetClasses.keys()
    
    for tid in g.calTargets:
        print("CalTarget {}".format(tid))
    allnodes += g.calTargets.keys()
    
    for aid,a in g.arcs.iteritems():
        n1,n2 = a.startnode,a.endnode
        if n1.id in allnodes and n2.id in allnodes:
                
            print("Arc {}".format(aid))

        
def printSciSurveyPlan(g):
    allnodes = []
    
    print("Visits: {}".format(g.visits))
    for tid in g.sciTargetClasses:
        print("SciTargetClass {}".format(tid))
    allnodes += g.sciTargetClasses.keys()
    
    for tid in g.sciTargets:
        print("SciTarget {}".format(tid))
    allnodes += g.sciTargets.keys()
    
    for tvid in g.targetVisits:
        print("TargetVisit {}".format(tvid))
    allnodes += g.targetVisits.keys()
    
    for aid,a in g.arcs.iteritems():
        n1,n2 = a.startnode,a.endnode
        if n1.id in allnodes and n2.id in allnodes:
                
            print("Arc {}".format(aid))

    
#printSciSurveyPlan(g)
# Plot targets
from numpy import array
from matplotlib import pyplot as plt


def plot_targets_in_fplane(targets, class_dict, RMAX = 15., CENTER = [-5.,90.] , xlim=[-200.,200.],ylim=[-200.,200.]):
    f = plt.figure(figsize=[5,5])

    sci = filter(lambda t : class_dict[t].startswith("sci"), targets )
    cal = filter(lambda t : class_dict[t].startswith("cal"), targets )
    sky = filter(lambda t : class_dict[t].startswith("sky"), targets )

    x = array( [targets[t][0] for t in sci])
    y = array( [targets[t][1] for t in sci])
    ii = ( (x-CENTER[0])**2 + (y-CENTER[1])**2 ) <= RMAX**2

    plt.plot(x[ii],y[ii],'x', ms=1, color='k')
    plt.plot(x[~ii],y[~ii],'x', ms=1, alpha=0.3, color='k')

    x = array( [targets[t][0] for t in cal])
    y = array( [targets[t][1] for t in cal])
    ii = ( (x-CENTER[0])**2 + (y-CENTER[1])**2 ) <= RMAX**2
    plt.plot(x[ii],y[ii],'+', ms=5, color='b')
    plt.plot(x[~ii],y[~ii],'x', ms=5, alpha=0.3, color='b')

    x = array( [targets[t][0] for t in sky])
    y = array( [targets[t][1] for t in sky])
    ii = ( (x-CENTER[0])**2 + (y-CENTER[1])**2 ) <= RMAX**2
    plt.plot(x[ii],y[ii],'+', ms=3, color='g')
    plt.plot(x[~ii],y[~ii],'x', ms=3, alpha=0.3, color='g')


    plt.ylabel('Y [mm]')
    plt.xlabel('X [mm]')
    plt.axis("equal")

    plt.xlim(xlim)
    plt.ylim(ylim)
