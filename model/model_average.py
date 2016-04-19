#!/usr/bin/env python
import sys
import argparse
import numpy as np
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
# Prevent propagating to higher loggers.
logger.propagate = 0
# Console log handler.
ch = logging.StreamHandler()
# Add formatter
FORMAT = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
formatter = logging.Formatter(FORMAT)
ch.setFormatter(formatter)
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)

def _read_model_txt(model_file):
    """
    type model_file : char (model name)

    mode file layout:
        n_model
        lon lat
        nlayer
        thick density vp vs
        ...
        lon lat
        nlayer
        thick density vp vs
        ...
    """
    model_all = []
    n_line = 0
    txt_file = open(model_file, "r")
    mlines = txt_file.readlines()

    n_model= int(mlines[0])
    n_line += 1
    print("%d models" % n_model)

    for i in range(1,n_model):
        model = {}
        lon, lat = mlines[n_line].strip().split()
        n_model_line = int( mlines[n_line + 1].strip() )

        model["lon"] = float(lon)
        model["lat"] = float(lat)
        model["layer"] = int(n_model_line)
        model["modelpara"] = []

        n_line += 2
        for j in range(n_model_line):
            thick, density, vp, vs = mlines[n_line + j].strip().split()
            model["modelpara"].append([float(thick), float(density),
                                       float(vp), float(vs)])

        n_line += n_model_line
        model_all.append(model)
    logger.debug("%d layers in model" % n_model_line)

    return model_all


def _average_ocean_model(models):
    """
    type model : model list
    """

    lon_min = -132.0
    lon_max = -120.0

    lat_min = 38.0
    lat_max = 50.0

    vs_upper_crust_min = 2.63 # km/s
    thick_sediment_max = 2.0 # km

    layer_max = 25

    model_mean = []
    model_err = []

    # calculate the average thickness of water layer (thickness of water must
    # larger than 0.5 km)
    # model contains layers [thick, density, vp, vs]
    water_thick = []
    water_density = []
    water_vp = []
    water_vs = []

    sediment_thick = []
    sediment_density = []
    sediment_vp = []
    sediment_vs = []

    model_thick = [[] for x in xrange(layer_max)]
    model_density = [[] for x in xrange(layer_max)]
    model_vp = [[] for x in xrange(layer_max)]
    model_vs = [[] for x in xrange(layer_max)]

    ocean_count = 0

    for model in models:
        model_layers = model["modelpara"]

        # average model for Juan de Fuca, Gorda and pacific plate
        if  model["lat"] < lat_max and model["lat"] > lat_min and\
            model["lon"] > lon_min and model["lon"] < lon_max:

            # for oceanic model only
            # shear wavespeed in water is 0.0
            # vs = 2.63 is a typical shear wavespeed of upper crust
            # thickness of sediments is around 2.0 km in study area
            # oceanic model has 22 layers in general if not in all cases
            if model_layers[0][3] <= 0.0 and\
                    model_layers[1][0] < thick_sediment_max and\
                    model_layers[1][3] < vs_upper_crust_min and\
                    model["layer"] == 22:

                water_thick.append(model_layers[0][0])
                water_density.append(model_layers[0][1])
                water_vp.append(model_layers[0][2])

                sediment_thick.append(model_layers[1][0])
                sediment_density.append(model_layers[1][1])
                sediment_vp.append(model_layers[1][2])
                sediment_vs.append(model_layers[1][3])

                for i in range(2, model["layer"]):
                    model_thick[i].append(model_layers[i][0])
                    model_density[i].append(model_layers[i][1])
                    model_vp[i].append(model_layers[i][2])
                    model_vs[i].append(model_layers[i][3])

    water_thick_mean = np.mean(water_thick)
    water_thick_err = np.std(water_thick)
    water_density_mean = np.mean(water_density)
    water_vp_mean = np.mean(water_vp)

    sediment_thick_mean = np.mean(sediment_thick)
    sediment_thick_err = np.std(sediment_thick)
    sediment_density_mean = np.mean(sediment_density)
    sediment_density_err = np.std(sediment_density)
    sediment_vp_mean = np.mean(sediment_vp)
    sediment_vp_err = np.std(sediment_vp)
    sediment_vs_mean = np.mean(sediment_vs)
    sediment_vs_err = np.std(sediment_vs)

    layer_water = [water_thick_mean, water_density_mean, water_vp_mean, 0.0]
    layer_water_err = [water_thick_err, 0.0, 0.0, 0.0]

    layer_sediment = [sediment_thick_mean, sediment_density_mean,
                      sediment_vp_mean, sediment_vs_mean]
    layer_sediment_err = [sediment_thick_err, sediment_density_err,
                      sediment_vp_err, sediment_vs_err]

    model_mean.append(layer_water)
    model_err.append(layer_water_err)

    model_mean.append(layer_sediment)
    model_err.append(layer_sediment_err)

    for i in range(2, layer_max):
        # make sure is not empty
        if model_thick[i]:
            model_mean.append( [np.mean(model_thick[i]),
                                   np.mean(model_density[i]),
                                   np.mean(model_vp[i]),
                                   np.mean(model_vs[i])] )

            model_err.append( [np.std(model_thick[i]),
                                  np.std(model_density[i]),
                                  np.std(model_vp[i]),
                                  np.std(model_vs[i])] )

    logger.debug("%d layers in mean model" % len(model_mean))

    return model_mean, model_err

def _write_model(model_list, model_fname):
    """
    type model_list :: [[layer 1], ..., [layer N], float
    type model_fname :: model name, char
    """
    # starting value of Qmu
    Qmu_0 = 100.0

    fh = open(model_fname, "w")

    for i in range(len(model_list)):
        fh.write("%8.3f %8.3f %8.3f %8.3f %8.3f\n" % (model_list[i][0],
                                                      model_list[i][1],
                                                      model_list[i][2],
                                                      model_list[i][3],
                                                      Qmu_0
                                                      ))
    fh.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', action='store', dest='modelfile',
                        required=True, help="input model file")
    parser.add_argument('-o', action='store', dest='model_mean',
                        required=True, help="output model file")
    parser.add_argument('-e', action='store', dest='model_err',
                        required=True, help="output model err")
    args = parser.parse_args()

    models = _read_model_txt(args.modelfile)
    model_mean, model_err = _average_ocean_model(models)

    _write_model(model_mean, args.model_mean)
    _write_model(model_err, args.model_err)

