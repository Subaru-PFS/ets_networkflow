from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import range
import matplotlib.pyplot as plt

import matplotlib.lines as mlines
from matplotlib import patches as mpatches
from matplotlib import collections
import numpy as np

from . import datamodel as dm


def plotSurveyPlan(g, name="", PLOTSCI=True, PLOTCAL=True, LABELFLOWS=False, ALTLABELS=None):

    def nlabel(ax, xy, text):
        y = xy[1] - 0.  # shift y-value for label so that it's below the artist
        x = xy[0] + 0.015  # shift x-value for label so that it's next to the artist
        ax.text(x, y, text, family='sans-serif', ha='left', va='center', fontsize=6)

    def plotNodes(ax, nodes, x, color, label, allnodes):
        # targets
        N = len(nodes)
        plt.text(x, 0., label, ha='center', va='top')
        if N == 0:
            return
        grid = [[x, 1. - (float(i)+1.)/float(N+1)] for i in range(N)]
        patches = []
        for i, (id, n) in enumerate(nodes.items()):
            n.px, n.py = grid[i]
            
            circle = mpatches.Circle(grid[i], .01, facecolor=color, edgecolor='grey', lw=1.)
            patches.append(circle)
            nlabel(ax, grid[i], id)

        collection = collections.PatchCollection(patches, match_original=True)

        ax.add_collection(collection)
        allnodes += list(nodes.keys())
        

    def plotArcs(ax, g, allnodes, LABELFLOWS, ALTLABELS):
        for a in g.arcs.values():
            n1, n2 = a.startnode, a.endnode
            if n1.id in allnodes and n2.id in allnodes:
                # add a line
                x, y = np.array([[n1.px, n2.px], [n1.py, n2.py]])

                lw = a.flow * 2 + 1.
                c = 'k'
                rflw = a.flow
                if rflw < 1e-6:
                    rflw = 0.
                if rflw > 0.:
                    c = 'blue'
                line = mlines.Line2D(x, y, lw=lw, alpha=.3, zorder=0, c=c)
                if LABELFLOWS:
                    l = a.id
                    if ALTLABELS is not None:
                        l = ALTLABELS[a.id]
                    ax.text((n1.px+n2.px)/2., (n1.py+n2.py)/2., "{}:{}".format(l, rflw), size=8)
                ax.add_line(line)

    if PLOTSCI:
        fig = plt.figure(figsize=[15, 15])
        # science targets
        ax = plt.subplot(111)  # note we must use plt.subplots, not plt.subplot
        allnodes = []  # keep track of nodes actually contained in the plot
        plotNodes(ax, nodes=g.cobras, x=0.9, label="cobras", color="#19967d", allnodes=allnodes)
        plotNodes(ax, nodes=g.cobraPointings, x=0.7, label="cobras\npointings", color="#99d5ca", allnodes=allnodes)
        plotNodes(ax, nodes=g.targetPointings, x=0.5, label="target\npointings", color="#ffbbb1", allnodes=allnodes)
        plotNodes(ax, nodes=g.sciTargets, x=0.3, label="science\ntargets", color="#ff8f80", allnodes=allnodes)
        plotNodes(ax, nodes=g.sciTargetClasses, x=0.1, label="science\ntarget\nclasses", color="#c92d39", allnodes=allnodes)
        plotArcs(ax, g, allnodes, LABELFLOWS, ALTLABELS)
        ax.set_ylim([-.0, 1.0])
        ax.set_xlim([-.0, 1.0])
        plt.axis('off')
        plt.savefig("{}_sci.pdf".format(name))

    if PLOTCAL:
        fig = plt.figure(figsize=[15, 15])
        # calibration targets
        ax = plt.subplot(111)  # note we must use plt.subplots, not plt.subplot
        allnodes = []
        plotNodes(ax, nodes=g.cobras, x=0.9, label="cobras", color="#19967d", allnodes=allnodes)
        plotNodes(ax, nodes=g.cobraPointings, x=0.7, label="cobras\npointings", color="#99d5ca", allnodes=allnodes)
        if len(g.calTargets) > 0:
            plotNodes(ax, nodes=g.calTargets, x=0.4, label="calib\ntargets", color="#ffeca9", allnodes=allnodes)
        if len(g.calTargetClasses) > 0:
            plotNodes(ax, nodes=g.calTargetClasses, x=0.2, label="calib.\ntarget\nclasses", color="#ef8d22", allnodes=allnodes)
        plotArcs(ax, g, allnodes, LABELFLOWS, ALTLABELS)
        ax.set_ylim([-.0, 1.0])
        ax.set_xlim([-.0, 1.0])
        plt.axis('off')

        plt.savefig("{}_cal.pdf".format(name))


def plotFocalPlane(g, pid, summary="", XC=0., YC=0., W=400., name="", figsize=[15,15]):

    def nlabel(ax, xy, text):
        y = xy[1] - 0.  # shift y-value for label so that it's below the artist
        x = xy[0] + 0.015  # shift x-value for label so that it's next to the artist
        ax.text(x, y, text, family='sans-serif', ha='left', va='center', fontsize=6)

    def plotNodes(ax, nodes, color, label, allnodes):
        # targets
        patches = []
        for i, (id, n) in enumerate(nodes.items()):
            n.px, n.py = n.getX(pid), n.getY(pid)
            circle = mpatches.Circle((n.getX(pid), n.getY(pid)), .75, facecolor=color, edgecolor='grey', lw=1.)
            patches.append(circle)

        collection = collections.PatchCollection(patches, match_original=True)
        ax.add_collection(collection)
        allnodes += list(nodes.keys())

    def plotArcs(ax, g, allnodes, pid):
        for a in g.arcs.values():
            n1, n2 = a.startnode, a.endnode
            if a.flow > 0.:
                alpha = 1.
            else:
                alpha = .2

            if type(n2) == dm.CobraPointing and n2.pointing == pid and type(n1) == dm.TargetPointing:
                    x, y = np.array([[n1.target.getX(pid), n2.cobra.getX(pid)], [n1.target.getY(pid), n2.cobra.getY(pid)]])

                    line = mlines.Line2D(x, y, lw=1, alpha=alpha, zorder=10)
                    ax.add_line(line)
            if type(n2) == dm.CobraPointing and n2.pointing == pid and type(n1) == dm.CalTarget:
                    x, y = np.array([[n1.getX(pid), n2.cobra.getX(pid)], [n1.getY(pid), n2.cobra.getY(pid)]])

                    line = mlines.Line2D(x, y, lw=1, alpha=alpha, zorder=10)
                    ax.add_line(line)

    fig = plt.figure(figsize=figsize)
    # science targets
    ax = plt.subplot(111)  # note we must use plt.subplots, not plt.subplot
    allnodes = []  # keep track of nodes actually contained in the plot
    plotNodes(ax, nodes=g.cobras, label="cobras", color="#19967d", allnodes=allnodes)
    plotNodes(ax, nodes=g.sciTargets, label="science\ntargets", color="#ff8f80", allnodes=allnodes)
    plotNodes(ax, nodes=g.calTargets, label="calib\ntargets",color="#ffeca9", allnodes=allnodes)
    plotArcs(ax, g, allnodes, pid)

    ax.set_ylim([YC-W/2., YC+W/2.])
    ax.set_xlim([XC-W/2., XC+W/2.])
    #plt.axis('off')
    plt.text(-0.1, -0.1, summary, ha='left', va='bottom', transform=ax.transAxes, fontsize=8)
    plt.text(0.5, 1.0, "pointing {}".format(pid), ha='center', va='top', transform=ax.transAxes)
    plt.savefig("{}_pointing{}_fp.pdf".format(name, pid))

    
def plotTargetDistribution(ra, dec, types, pointings, target_fplane_pos, class_dict):
    """
    Plots the distribution of targets on the sky and in the focal plane for all dither positions.
    """

    # plot targets on sky
    f = plt.figure(figsize=[15,15])
    plt.title("Sky")
    ax = plt.subplot(111)
    ax.set_facecolor((.95,.95,1.))
    ii_sci = list( map( lambda x : x.startswith('sci') , types ) )
    ii_sky = list( map( lambda x : x.startswith('sky') , types ) )
    ii_cal = list( map( lambda x : x.startswith('cal') , types ) )
        
    plt.plot(np.array(ra)[ii_sky],np.array(dec)[ii_sky],'b.' , label='sky', ms=2)
    plt.plot(np.array(ra)[ii_cal],np.array(dec)[ii_cal],'ro' , label='cal. star', ms=2)
    plt.plot(np.array(ra)[ii_sci],np.array(dec)[ii_sci],'.', ms=1, label='science')
    #plt.axis('equal')
    l = plt.legend()
    l.draw_frame(False)
    
    for p in pointings.values():
        plt.plot(p[0], p[1],'yx')
    plt.xlabel("RA [Deg]")
    plt.ylabel("DEC [Deg]")


    f = plt.figure(figsize=[17,17])
    M = int(np.ceil( np.sqrt(len(pointings)) ))
    # plot targets in focal plane
    i = 0
    for pid,(pointing_RA,pointing_DEC) in pointings.items():
        i += 1


        ax = plt.subplot(M,M,i)
        plt.title("focal plane Pointing {}".format(pid))
        ax.set_facecolor((.95,.95,.95))

        targets = target_fplane_pos[pid]
        tclasses = [class_dict[tid] for tid in target_fplane_pos[pid]]

        
        txx = np.array( [t[0] for tid, t in targets.items()] )
        tyy = np.array( [t[1] for tid, t in targets.items()] )
        
        RMAX = 300.
        jj = txx**2. + tyy**2. < RMAX**2.

        _ii_sci = list( map( lambda x : x.startswith('sci') , tclasses ) )
        _ii_sky = list( map( lambda x : x.startswith('sky') , tclasses ) )
        _ii_cal = list( map( lambda x : x.startswith('cal') , tclasses ) )

        plt.plot(txx[_ii_sky * jj],tyy[_ii_sky * jj],'b.' , label='sky', ms=2)
        plt.plot(txx[_ii_cal * jj],tyy[_ii_cal * jj],'ro' , label='cal. star', ms=2)
        plt.plot(txx[_ii_sci * jj],tyy[_ii_sci * jj],'.', ms=1, label='science')


        plt.axis('equal')
        l = plt.legend()
        l.draw_frame(False)
        plt.xlabel("x [mm]")
        plt.ylabel("y [mm]")
    f.tight_layout()
    
    
    

def filterFoV(xx,yy, maxx=20., maxy=20.):
    ii  = xx >= -maxx
    ii *= xx <=  maxx
    ii *= yy >= -maxy
    ii *= yy <=  maxy
    return ii

def my_circle_scatter(axes, x_array, y_array, radius, **kwargs):
    for x, y, r in zip(x_array, y_array, radius):
        circle = plt.Circle((x,y), radius=r, edgecolor='k', facecolor='grey', **kwargs)
        axes.add_patch(circle)
    return True

def plotFP(cobras,dots,targets,maxx=20.,maxy=20.):
    """
    This plots the focal plane, with all the target positions at all
    dither positions. Can take very long, shoudl only be done with a subset of the total 
    targets.
    """
    f = plt.figure(figsize=[10,10])
    ax = plt.subplot(1,1,1)
    
    # rearrange dot and cobra information for easy plotting
    dxx = np.array( [d[0] for did, d in dots.items()] )
    dyy = np.array( [d[1] for did, d in dots.items()] )
    drr = np.array( [d[2] for did, d in dots.items()] )
    ciddd = np.array( [did for did, d in dots.items()] )

    cxx = np.array( [c[0] for cid, c in cobras.items()] )
    cyy = np.array( [c[1] for cid, c in cobras.items()] )
    ciddd = np.array( [cid for cid, t in cobras.items()] )
    
    # plot cobras
    ii = filterFoV(cxx,cyy,maxx=maxx,maxy=maxy)
    plt.plot(cxx[ii],cyy[ii],'.', c='orange')
    ii = filterFoV(cxx,cyy )
    for x,y,l in zip(cxx[ii],cyy[ii],ciddd[ii]):
        plt.text(x,y,l)

    # plot dots
    ii = filterFoV(dxx,dyy,maxx=maxx,maxy=maxy)
    my_circle_scatter(ax, dxx[ii],dyy[ii], drr[ii] )
    

    N = len(target_fplane_pos)
    for i,pid in enumerate(target_fplane_pos):
        alpha = 1. - float(i)/float( N )

        # rearrange target information for easy plotting
        txx = np.array( [t[0] for tid, t in target_fplane_pos[pid].items()] )
        tyy = np.array( [t[1] for tid, t in target_fplane_pos[pid].items()] )
        tiddd = np.array( [tid for tid, t in target_fplane_pos[pid].items()] )

        # plot targets
        ii = filterFoV(txx,tyy,maxx=maxx,maxy=maxy)
        plt.plot(txx[ii],tyy[ii],'bo', alpha=alpha, label='pointing {}'.format(pid))

        for x,y,l in zip(txx[ii],tyy[ii],tiddd[ii]):
            plt.text(x,y,l, size=8)

        # plot visibility arcs
        if True:
            for tid,cc in visibilities[pid].items():
                tx,ty = target_fplane_pos[pid][tid]
                for c in cc:
                    cx,cy = cobras[c]
                    #print c
                    #plt.text(cx,cy,"{}".format(c) )
                    plt.plot([cx,tx],[cy,ty],'k-', alpha=alpha)


        
            
       #plt.text(zip(txx[ii],tyy[ii]),tiddd[ii], size=6)   

    plt.axis('equal')
    plt.xlim([-maxx,maxx])
    plt.ylim([-maxy,maxy])
    plt.legend()
    
    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
