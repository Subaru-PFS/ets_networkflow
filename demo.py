import pyETS
import libETS.data_model as dm
import numpy as np

# define locations of the input files
catalog_path = "/home/martin/codes/ets_fiber_assigner/pfs_target_list"
fscience_targets = catalog_path+"/pfs_preliminary_target_cosmology.dat"
#fscience_targets = catalog_path+"/pfs_preliminary_target_galaxy.dat"
fcal_stars       = catalog_path+"/pfs_preliminary_target_cosmology_fcstars.dat"
fsky_pos         = catalog_path+"/pfs_preliminary_target_cosmology_sky.dat"

# read all targets into a sigle list, giving them their proper types
tgt = dm.readScientificFromFile(fscience_targets, "sci")
tgt += dm.readCalibrationFromFile(fcal_stars, "cal")
tgt += dm.readCalibrationFromFile(fsky_pos, "sky")

# get a complete, idealized focal plane configuration
cobras = dm.getFullFocalPlane()

# build reduced Cobra list to speed up calculation
cobras = [c for c in cobras if abs(c.center) < 200]

# point the telescope at the center of all science targets
raTel, decTel = dm.telescopeRaDecFromFile(fscience_targets)
posang = 0.
otime = "2016-04-03T08:00:00Z"
telescope = dm.Telescope(cobras, 1., raTel, decTel, posang, otime)

classdict = {}
classdict["sci_P1"] = {"nonObservationCost": 100, "partialObservationCost": 1e9}
classdict["sci_P2"] = {"nonObservationCost": 90, "partialObservationCost": 1e9}
classdict["sci_P3"] = {"nonObservationCost": 80, "partialObservationCost": 1e9}
classdict["sci_P4"] = {"nonObservationCost": 70, "partialObservationCost": 1e9}
classdict["sci_P5"] = {"nonObservationCost": 60, "partialObservationCost": 1e9}
classdict["sci_P6"] = {"nonObservationCost": 50, "partialObservationCost": 1e9}
classdict["sci_P7"] = {"nonObservationCost": 40, "partialObservationCost": 1e9}
classdict["sky"] = {"numRequired": 1}
classdict["cal"] = {"numRequired": 1}

res = telescope.observeWithNetflow(tgt, classdict, 3, 300.)
print(res)
