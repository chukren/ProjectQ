#!/usr/bin/env python

import sys
import argparse
import numpy as np
import logging


def _update_prem_model(filename, model, outputmodel):
    """
    type filename :: char (name of input prem files)
    type model :: dictionary of model including
        thick
        density
        vp
        vs
        Qkappa
        Qmu

    read prem model
        ISOTROPIC PREM
        0 1.0 1
            150 11 32 6371
        radius density vpv vsv Qk Qmu vph vsh eta
        ...
    """
    output = open(outputmodel, "w")

    txt_file = open(filename, "r")
    mlines = txt_file.readlines()
    nlayer = len(mlines)
    nheader = 3

    # sewing layer should be less than 10 km
    gap = 10000.0

    model_thick = model["thick"]
    model_density = model["density"]
    model_vp = model["vp"]
    model_vs = model["vs"]
    model_Qkappa = model["Qkappa"]
    model_Qmu = model["Qmu"]

    new_model = []

    # PREM format
    formatstr = \
        "%10.1f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f"

    # count layers of new model
    count_layer = 0

    for i in range(nlayer):
        # model starts from fourth line
        if i > nheader - 1:
            radius, density, vpv, vsv, Qkappa, Qmu, vph, vsh, eta =\
                mlines[i].rstrip().split()

            # copy cores and lower mantle
            if float(radius) <= (model_thick[-1] - gap):
                new_model.append(mlines[i].rstrip())
                count_layer += 1

            # replace upper mantle with inversion model
            else:
                new_model.append(formatstr % (model_thick[-1], float(density),
                                             float(vpv), float(vsv),
                                             float(Qkappa), float(Qmu),
                                             float(vph), float(vsh),
                                             float(eta)))
                count_layer += 1
                # dump new model in PREM format
                for x in xrange(len(model_thick)):
                    thick = model_thick.pop()
                    density = model_density.pop()
                    vp = model_vp.pop()
                    vs = model_vs.pop()
                    Qkappa = model_Qkappa.pop()
                    Qmu = model_Qmu.pop()
                    new_model.append(formatstr % (thick, density, vp, vs,
                                                  Qkappa, Qmu, vp, vs, 1.0))
                    count_layer += 1
                break

    # read and update header -- 3 lines
    for i in range(3):
        # update layers
        if i == 2:
            tmp, i_core, o_core, radius_km = mlines[i].rstrip().split()
            output.write("%12d %12s %12s %7s\n" % (count_layer, i_core, o_core,
                                                   radius_km))
        else:
            output.write("%s" % mlines[i])

    # dump model
    for line in new_model:
        output.write("%s\n" % line)

    output.close()


def _read_model_text(filename):
    """
    read text model
    """
    r_earth_m = 6371000.0
    eta = 1.0

    txt_file = open(filename, "r")
    mlines = txt_file.readlines()
    nlayer = len(mlines)

    model_thick = []
    model_density = []
    model_vp = []
    model_vs = []
    model_Qkappa = []
    model_Qmu = []

    thick_m = 0.0

    for i in range(nlayer):
        thick_str, density_str, vp_str, vs_str, Qmu_str = \
            mlines[i].rstrip().split()
        thick_m = thick_m + float(thick_str) * 1000.0
        thick = r_earth_m - thick_m
        density = float(density_str) * 1000.0
        vp = float(vp_str) * 1000.0
        vs = float(vs_str) * 1000.0
        Qmu = float(Qmu_str)
        Qkappa = 57827.0


        # values at top of a layer
        if i == 0:
            model_thick.append(r_earth_m)
        else:
            model_thick.append(model_thick[-1])
        model_density.append(density)
        model_vp.append(vp)
        model_vs.append(vs)
        model_Qkappa.append(Qkappa)
        model_Qmu.append(Qmu)

        # values at bottom of a layer
        model_thick.append(thick)
        model_density.append(density)
        model_vp.append(vp)
        model_vs.append(vs)
        model_Qkappa.append(Qkappa)
        model_Qmu.append(Qmu)

    print("Total thickness of this model is %9.1f km"
          % ((r_earth_m - model_thick[-1])/1000.0) )


    model = {}
    model["thick"] = model_thick
    model["density"] = model_density
    model["vp"] = model_vp
    model["vs"] = model_vs
    model["Qkappa"] = model_Qkappa
    model["Qmu"] = model_Qmu

    return model


def _write_model_prem(filename):
    """
    write output model
    """
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='inputmodel',
                        required=True, help="input model file")
    parser.add_argument('-p', action='store', dest='prem',
                        required=True, help="input PREM model file")
    parser.add_argument('-o', action='store', dest='outputmodel',
                        required=True, help="output model file")

    args = parser.parse_args()
    model = _read_model_text(args.inputmodel)
    _update_prem_model(args.prem, model, args.outputmodel)
