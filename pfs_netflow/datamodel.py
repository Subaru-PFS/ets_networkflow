from builtins import object
# # Datamodel for the flow network

from collections import OrderedDict
from numpy import inf


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


class Sink(Node):
    def __init__(self, id):
        super(Sink, self).__init__(id)


class Cobra(Node):
    def __init__(self, id, x, y):
        super(Cobra, self).__init__(id)
        self.x = x
        self.y = y


class CobraVisit(Node):
    def __init__(self, id, cobra, visit):
        super(CobraVisit, self).__init__(id)
        self.visit = visit
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
        self.cobraVisits = OrderedDict()
        self.targets = {}
        self.calTargets = OrderedDict()
        self.sciTargets = OrderedDict()
        self.targetVisits = OrderedDict()
        self.calTargetClasses = OrderedDict()
        self.sciTargetClasses = OrderedDict()
        self.sinks = OrderedDict()
        # arcs
        self.targetClassToTargetArcs = OrderedDict()
        self.targetToTargetVisitArcs = OrderedDict()
        self.targetVisitToCobraVisitArcs = OrderedDict()
        self.cobraVisitToCobraArcs = OrderedDict()
        self.cobraToSinkArcs = OrderedDict()
        self.cobraVisitToCobraArcs = OrderedDict()
        self.overflowArcs = OrderedDict()

    def add_node(self, node):
        super(SurveyPlan, self).add_node(node)
        if type(node) == Cobra:
            self.cobras[node.id] = node
        elif type(node) == CobraVisit:
            self.cobraVisits[node.id] = node
        elif type(node) == CalTarget:
            self.calTargets[node.id] = node
            self.targets[node.id] = node
        elif type(node) == SciTarget:
            self.sciTargets[node.id] = node
            self.targets[node.id] = node
        elif type(node) == TargetVisit:
            self.targets[node.target.id].targetVisits[node.visit] = node
            self.targetVisits[node.id] = node
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
        elif type(arc) == TargetToTargetVisitArc:
            self.targetToTargetVisitArcs[arc.id] = arc
        elif type(arc) == TargetVisitToCobraVisitArc:
            self.targetVisitToCobraVisitArcs[arc.id] = arc
        elif type(arc) == CobraVisitToCobraArc:
            self.cobraVisitToCobraArcs[arc.id] = arc
        elif type(arc) == CobraToSinkArc:
            self.cobraToSinkArcs[arc.id] = arc
        elif type(arc) == CobraVisitToCobraArc:
            self.cobraVisitToCobraArcs[arc.id] = arc
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


class TargetToTargetVisitArc(Arc):
    def __init__(self, startnode, endnode):
        super(TargetToTargetVisitArc, self).__init__(startnode, endnode)


class TargetVisitToCobraVisitArc(Arc):
    def __init__(self, startnode, endnode):
        super(TargetVisitToCobraVisitArc, self).__init__(startnode, endnode)
        self.d = None  # Distance from target to cobra in mm.
        self.visit = None


class CobraVisitToCobraArc(Arc):
    def __init__(self, startnode, endnode):
        super(CobraVisitToCobraArc, self).__init__(startnode, endnode)


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
    def __init__(self, id):
        super(CalTargetClass, self).__init__(id)
        self.visit = None


class SciTargetClass(TargetClass):
    def __init__(self, id):
        super(SciTargetClass, self).__init__(id)
        # Describes cost of partial completion.
        # Typically this should be a higher cost
        # than not observing the target at all.
        self.cost_partial_compl = None


class Target(Node):
    def __init__(self, id, x, y, gain=3):
        super(Target, self).__init__(id)
        self.x = x
        self.y = y
        self.gain = gain  # number of required exposures
        self.collision_group = None
        self.targetVisits = {}


class SciTarget(Target):
    def __init__(self, id, x, y, gain=3):
        super(SciTarget, self).__init__(id, x, y, gain)


class CalTarget(Target):
    """
    Calibration targets are different in the sense that they can be reobserved
    an arbitrary number of times while their required number of visits
    is always one. A reobservation will therefore never be enforced.
    """
    def __init__(self, id, x, y, visit, gain=3):
        super(CalTarget, self).__init__(id, x, y, gain)
        self.visit = visit


class TargetVisit(Node):
    def __init__(self, id, target, visit):
        super(TargetVisit, self).__init__(id)
        self.visit = visit
        self.target = target
