#!/usr/bin/env python
import sys
import argparse
import json
import numpy as np

ftemplate = "PAR_RAYLEIGH"
model_list = "model.list"

def _read_text(filelist):
    with open(filelist, "r") as fh:
        models = fh.readlines()
    return models

def _read_json_period(period_list_json):
    with open(period_list_json, "r") as fh:
        period_band = json.load(fh)
    return period_band

def read_modeoutput(fmodeoutput, period_band_list):
    mode_output = fmodeoutput.strip()
    with open(mode_output) as fh:
        lines = fh.readlines()

    atten_coef_period = {}

    for period in period_band_list:
        fmin = period_band_list[period]["min"]
        fmax = period_band_list[period]["max"]

        atten_coef = []

        for line in lines:
            if "0 R" in line and "skip" not in line:
                (mode_n, mode_type, mode_l, mode_freq, mode_period,
                 phsv, grpv, q_rayl, qaulity, dummy) = line.strip().split()

                atten_coef_temp = \
                    np.pi * float(mode_freq) / (float(grpv) * float(q_rayl))

                if fmin <= float(mode_freq) <= fmax:
                    atten_coef.append(atten_coef_temp)
                    # print period, mode_freq, q_rayl, fmin, fmax

        atten_coef_mean = np.mean(atten_coef)
        atten_coef_stdv = np.std(atten_coef)

        print("%s %f %f" % (period, atten_coef_mean, atten_coef_stdv))

        atten_coef_period[period] = atten_coef_mean

    return atten_coef_period

def read_param(ftemplate):
    with open(ftemplate, "r") as fh:
        lines = fh.readlines()
    return lines

def read_modellist(model_list):
    with open(model_list, "r") as fh:
        lines = fh.readlines()
    return lines

def generate_param(param_template, model, mode_output, fparam):
    foutput = open(fparam, "w")
    foutput.write(model)
    foutput.write(mode_output)
    for i in range(2, len(param_template)):
        foutput.write(param_template[i])
    foutput.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', action='store', dest='period_band',
                        required=True, help="input model file")
    parser.add_argument('-f', action='store', dest='mode_output_list',
                        required=True, help="input mode output files")
    parser.add_argument('-r', action='store', dest='reference_mode_output',
                        required=True, help="input reference mode output file")
    parser.add_argument('-g', action='store', dest='gradient_output',
                        required=True, help="input reference mode output file")
    args = parser.parse_args()

    period_list = _read_json_period(args.period_band)
    mode_outputs = _read_text(args.mode_output_list)

    # calculate atten_coef at each period for perturbed models
    atten_coef_matrix = {}
    for model in mode_outputs:
        (dummy, layer_id) =  model.rstrip().split("_")
        print ("%s" % layer_id)
        atten_coef_period = read_modeoutput(model, period_list)
        atten_coef_matrix[layer_id] = atten_coef_period

    # calculate atten_coef at each period for reference model
    model = args.reference_mode_output
    (dummy, layer_id) =  model.rstrip().split("_")
    print ("%s" % layer_id)
    atten_coef_period_ref = read_modeoutput(model, period_list)

    # fixed perturbation
    delta_q = 50.

    # calculate gradient
    gradient = {}
    for layer_id in atten_coef_matrix:
        synt = atten_coef_matrix[layer_id]
        gradient_period = {}
        for period in atten_coef_period_ref:
            atten_p = synt[period]
            atten_r = atten_coef_period_ref[period]
            g = (atten_p - atten_r) / delta_q
            gradient_period[period] = g
            print ("%s %s %e %e %e" % (layer_id, period, atten_p, atten_r, g))
        gradient[layer_id] = gradient_period

    print json.dumps(gradient, indent=2, sort_keys=True)

