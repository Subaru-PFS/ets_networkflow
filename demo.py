import pyETS
import libETS.data_model as dm
import numpy as np

class SkyCalibTarget(dm.CalibTarget):
    nonObservationCost = 10000.

    @staticmethod
    def numRequired():
        return 1
    @staticmethod
    def classname():
        return "sky_P1"

class StarCalibTarget(dm.CalibTarget):
    nonObservationCost = 10000.

    @staticmethod
    def numRequired():
        return 1
    @staticmethod
    def classname():
        return "cal_P1"

# ... maybe other calibration classes

def telescopeRaDecFromFile(file):
    with open(file) as f:
        ras=[]
        decs=[]
        ll=f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                ra,dec = (float(tt[1]), float(tt[2]))
                ras.append(ra)
                decs.append(dec)
    return float(np.average(ras)), float(np.average(decs))


def readScientificFromFile(file):
    with open(file) as f:
        res = []
        ll=f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_,ra,dec,tm,pri = (str(tt[0]), float(tt[1]), float(tt[2]),
                                       float(tt[3]), int(tt[4]))
                res.append(dm.ScienceTarget(id_, ra, dec, tm, pri))
    return res


def readCalibrationFromFile(file, cls):
    with open(file) as f:
        res = []
        ll=f.readlines()
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_,ra,dec = (str(tt[0]), float(tt[1]), float(tt[2]))
                res.append(cls(id_, ra, dec))
    return res


def getFocalPlane():
    cobras = pyETS.getAllCobras()
    res = []
    for i in range(len(cobras)):
        ID = "{}".format(i)
        res.append(dm.Cobra(ID, cobras[i][0], cobras[i][3], cobras[i][4],
                            cobras[i][1], cobras[i][2]))
    return res


catalog_path = "/home/martin/codes/ets_fiber_assigner/pfs_target_list"
fscience_targets = catalog_path+"/pfs_preliminary_target_cosmology.dat"
fcal_stars       = catalog_path+"/pfs_preliminary_target_cosmology_fcstars.dat"
fsky_pos         = catalog_path+"/pfs_preliminary_target_cosmology_sky.dat"

tgt = readScientificFromFile(fscience_targets)
tgt += readCalibrationFromFile(fcal_stars, StarCalibTarget)
tgt += readCalibrationFromFile(fsky_pos, SkyCalibTarget)

cobras = getFocalPlane()

# build reduced Cobra list to speed up calculation
cobras = [c for c in cobras if abs(c.center) < 20]

# point the telescope at the center of all science targets
raTel, decTel = telescopeRaDecFromFile(fscience_targets)
posang = 0.
otime = "2016-04-03T08:00:00Z"
telescope = dm.Telescope(cobras, 1., raTel, decTel, posang, otime)

#res=telescope.observeWithETS(tgt, 20, 300, "draining")
res=telescope.observeWithNetflow(tgt, 21, 300.)
print(res)
