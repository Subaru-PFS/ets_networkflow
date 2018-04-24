from builtins import object
# # Datamodel for the flow network

from collections import OrderedDict
from numpy import inf, nan

#class FocalPlanePos(object):
#    def __init__(self, x,y):
#        self.x = x
#        self.y = y

class NetworkElement(object):
    def __init__(self):
        pass


# Node types ...
class Node(NetworkElement):
    def __init__(self, id):
        super(Node, self).__init__()
        self.id = id
        self.capacity = inf
        self.gain = 1
        self.inarcs = []
        self.outarcs = []
        self.supply = 0
        self.fplane_pos = (nan,nan)
        
    def getX(self, pid):
        return self.fplane_pos[0]
        
    def getY(self, pid):
        return self.fplane_pos[1]

class Sink(Node):
    def __init__(self):
        super(Sink, self).__init__("SINK")


class Cobra(Node):
    @staticmethod
    def getID(cid):
        return "C_{}".format(cid)

    def __init__(self, cid, fplane_pos):
        super(Cobra, self).__init__(self.getID(cid))
        self.fplane_pos = fplane_pos
        
    def getX(self, pid=None):
        """
        pid is only a fake parameter here to give the
        cobra's getX() the same signature as the target's 
        getX(). Other that for the targets, the cobras
        x,y don't cahnge from pointing to pointings.
        """
        return self.fplane_pos[0]
        
    def getY(self, pid=None):
        """
        pid is only a fake parameter here to give the
        cobra's getY() the same signature as the target's 
        getY(). Other that for the targets, the cobras
        x,y don't cahnge from pointing to pointings.
        """
        return self.fplane_pos[1]


class CobraPointing(Node):
    @staticmethod
    def getID(cid, pointing):
        return "C_{}_v{}".format(cid, pointing)

    def __init__(self, cid, cobra, pointing):
        super(CobraPointing, self).__init__(self.getID(cid, pointing))
        self.pointing = pointing
        self.cobra = cobra
      

class Network(object):
    """
    Network flow master class.
    """
    def __init__(self):
        self.nodes = {}
        self.arcs = {}

    def add_node(self, node):
        self.nodes[node.id] = node

    def add_arc(self, arc):
        self.arcs[arc.id] = arc
        self.nodes[arc.startnode.id].outarcs.append(arc)
        self.nodes[arc.endnode.id].inarcs.append(arc)


class SurveyPlan(Network):
    def __init__(self):
        super(SurveyPlan, self).__init__()
        # These are just handy for quick access to a certain type of
        # node or arc nodes
        self.cobras = OrderedDict()
        self.cobraPointings = OrderedDict()
        self.targets = {}
        self.calTargets = OrderedDict()
        self.sciTargets = OrderedDict()
        self.targetPointings = OrderedDict()
        self.calTargetClasses = OrderedDict()
        self.sciTargetClasses = OrderedDict()
        self.sinks = OrderedDict()
        # arcs
        self.targetClassToTargetArcs = OrderedDict()
        self.targetToTargetPointingArcs = OrderedDict()
        self.targetPointingToCobraPointingArcs = OrderedDict()
        self.cobraPointingToCobraArcs = OrderedDict()
        self.cobraToSinkArcs = OrderedDict()
        self.cobraPointingToCobraArcs = OrderedDict()
        self.overflowArcs = OrderedDict()

    def add_node(self, node):
        super(SurveyPlan, self).add_node(node)
        if type(node) == Cobra:
            self.cobras[node.id] = node
        elif type(node) == CobraPointing:
            self.cobraPointings[node.id] = node
        elif type(node) == CalTarget or type(node) == StarCalTarget or type(node) == SkyCalTarget:
            self.calTargets[node.id] = node
            self.targets[node.id] = node
        elif type(node) == SciTarget:
            self.sciTargets[node.id] = node
            self.targets[node.id] = node
        elif type(node) == TargetPointing:
            self.targets[node.target.id].targetPointings[node.pointing] = node
            self.targetPointings[node.id] = node
        elif type(node) == TargetClass:
            self.targetClasses[node.id] = node
        elif type(node) == CalTargetClass:
            self.calTargetClasses[node.id] = node
        elif type(node) == SciTargetClass:
            self.sciTargetClasses[node.id] = node
        elif type(node) == Sink:
            self.sinks[node.id] = node
        else:
            raise Exception("Unknown node type.")

    def get_overflowarcs(self):
        return [a for a in self.arcs if a == OverflowArc]

    def add_arc(self, arc):
        super(SurveyPlan, self).add_arc(arc)
        if type(arc) == TargetClassToTargetArc:
            self.targetClassToTargetArcs[arc.id] = arc
        elif type(arc) == TargetToTargetPointingArc:
            self.targetToTargetPointingArcs[arc.id] = arc
        elif type(arc) == TargetPointingToCobraPointingArc:
            self.targetPointingToCobraPointingArcs[arc.id] = arc
        elif type(arc) == CobraPointingToCobraArc:
            self.cobraPointingToCobraArcs[arc.id] = arc
        elif type(arc) == CobraToSinkArc:
            self.cobraToSinkArcs[arc.id] = arc
        elif type(arc) == CobraPointingToCobraArc:
            self.cobraPointingToCobraArcs[arc.id] = arc
        elif type(arc) == OverflowArc:
            self.overflowArcs[arc.id] = arc
        else:
            raise Exception("Unknown arc type: {}".format(type(arc)))


# Arc types ...
class Arc(NetworkElement):
    def __init__(self, startnode, endnode):
        super(Arc, self).__init__()
        self.id = "{}={}".format(startnode.id, endnode.id)
        self.startnode = startnode
        self.endnode = endnode
        self.capacity = 1
        self.cost = 0
        self.flow = 0


class TargetClassToTargetArc(Arc):
    def __init__(self, startnode, endnode):
        super(TargetClassToTargetArc, self).__init__(startnode, endnode)


class TargetToTargetPointingArc(Arc):
    def __init__(self, startnode, endnode):
        super(TargetToTargetPointingArc, self).__init__(startnode, endnode)


class TargetPointingToCobraPointingArc(Arc):
    def __init__(self, startnode, endnode):
        super(TargetPointingToCobraPointingArc, self).__init__(startnode, endnode)
        self.d = None  # Distance from target to cobra in mm.
        self.pointing = None


class CobraPointingToCobraArc(Arc):
    def __init__(self, startnode, endnode):
        super(CobraPointingToCobraArc, self).__init__(startnode, endnode)


class CobraToSinkArc(Arc):
    def __init__(self, startnode, endnode):
        super(CobraToSinkArc, self).__init__(startnode, endnode)
        self.cost = 0
        self.capacity = inf


class OverflowArc(Arc):
    def __init__(self, cost, startnode, endnode):
        super(OverflowArc, self).__init__(startnode, endnode)
        self.cost = cost
        self.capacity = inf


class TargetClass(Node):
    def __init__(self, id):
        super(TargetClass, self).__init__(id)
        self.ID = id
        self.capacity = inf
        self.targets = OrderedDict()
        # Cost of non-observation
        self.cost = None
        # How many objects of this class need to be observed.
        # For science object we probably generally want as many as we can?
        # For calibration objects only a subset N out of M may be required.
        self.num_required = inf

    def add_target(self, t):
        self.targets[t.id] = t


class CalTargetClass(TargetClass):
    @staticmethod
    def getID(tc, pointing):
        return "TClass_{}_v{}".format(tc, pointing)

    def __init__(self, tc, pointing):
        super(CalTargetClass, self).__init__(self.getID(tc, pointing))
        self.pointing = None


class SciTargetClass(TargetClass):
    @staticmethod
    def getID(tc):
        return "TClass_{}".format(tc)

    def __init__(self, tc):
        super(SciTargetClass, self).__init__(self.getID(tc))
        # Describes cost of partial completion.
        # Typically this should be a higher cost
        # than not observing the target at all.
        self.cost_partial_compl = None


class Target(Node):
    def __init__(self, id, fplane_positions, gain=3):
        super(Target, self).__init__(id)
        self.fplane_positions = fplane_positions
        self.gain = gain  # number of required exposures
        self.collision_group = None
        self.targetPointings = {}
        self.target_class = None

    def getX(self, pid):
        return self.fplane_positions[pid][0]
        
    def getY(self, pid):
        return self.fplane_positions[pid][1]
    
class SciTarget(Target):
    @staticmethod
    def getID(tid):
        return "T_{}".format(tid)

    def __init__(self, tid, fplane_positions, gain=3):
        super(SciTarget, self).__init__(self.getID(tid), fplane_positions, gain)


class CalTarget(Target):
    """
    Calibration targets are different in the sense that they can be reobserved
    an arbitrary number of times while their required number of pointings
    is always one. A reobservation will therefore never be enforced.
    """
    @staticmethod
    def getID(tid, pointing):
        return "T_{}_v{}".format(tid, pointing)

    def __init__(self, tid, fplane_positions, pointing, gain=3):
        super(CalTarget, self).__init__(self.getID(tid, pointing), fplane_positions, gain)
        self.pointing = pointing

class SkyCalTarget(CalTarget):
    """
    Really the same functionality as CalTarget, but help later discrimination.
    """
    def __init__(self, tid, fplane_positions, pointing, gain=3):
        super(SkyCalTarget, self).__init__(tid, fplane_positions, pointing, gain)
        self.pointing = pointing

class StarCalTarget(CalTarget):
    """
    Really the same functionality as CalTarget, but help later discrimination.
    """
    def __init__(self, tid, fplane_positions, pointing, gain=3):
        super(StarCalTarget, self).__init__(tid, fplane_positions, pointing, gain)
        self.pointing = pointing


class TargetPointing(Node):
    @staticmethod
    def getID(tid, pointing):
        return "T_{}_v{}".format(tid, pointing)

    def __init__(self, tid, target, pointing):
        super(TargetPointing, self).__init__(self.getID(tid, pointing))
        self.pointing = pointing
        self.target = target
