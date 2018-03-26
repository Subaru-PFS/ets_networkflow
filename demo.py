import pyETS
import libETS.data_model as dm
import numpy as np


class SkyCalibTarget(dm.CalibTarget):
    """Class for sky calibration targets"""
    nonObservationCost = 10000.

    @staticmethod
    def numRequired():
        return 1
    @staticmethod
    def classname():
        return "sky_P1"

class StarCalibTarget(dm.CalibTarget):
    """Class for star calibration targets"""
    nonObservationCost = 10000.

    @staticmethod
    def numRequired():
        return 1
    @staticmethod
    def classname():
        return "cal_P1"

# ... maybe other calibration classes

# define locations of the input files
catalog_path = "/home/martin/codes/ets_fiber_assigner/pfs_target_list"
fscience_targets = catalog_path+"/pfs_preliminary_target_cosmology.dat"
fcal_stars       = catalog_path+"/pfs_preliminary_target_cosmology_fcstars.dat"
fsky_pos         = catalog_path+"/pfs_preliminary_target_cosmology_sky.dat"

# read all targets into a sigle list, giving them their proper types
tgt = dm.readScientificFromFile(fscience_targets)
tgt += dm.readCalibrationFromFile(fcal_stars, StarCalibTarget)
tgt += dm.readCalibrationFromFile(fsky_pos, SkyCalibTarget)

# get a complete, idealized focal plane configuration
cobras = dm.getFullFocalPlane()

# build reduced Cobra list to speed up calculation
cobras = [c for c in cobras if abs(c.center) < 200]

# point the telescope at the center of all science targets
raTel, decTel = dm.telescopeRaDecFromFile(fscience_targets)
posang = 0.
otime = "2016-04-03T08:00:00Z"
telescope = dm.Telescope(cobras, 1., raTel, decTel, posang, otime)

res=telescope.observeWithNetflow(tgt, 9, 300.)
#print(res)
