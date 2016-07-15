#!/usr/bin/env python
import sys
import argparse
import json
import numpy as np

ftemplate = "PAR_RAYLEIGH"
model_list = "model.list"
model_index = "/home/youyir/Cascadia/ProjectQ/model/model.index"
data_index = "/home/youyir/Cascadia/ProjectQ/data/data.index"

def _read_text(filelist):
    with open(filelist, "r") as fh:
        lines = fh.readlines()
    return lines

def _read_json_period(period_list_json):
    with open(period_list_json, "r") as fh:
        period_band = json.load(fh)
    return period_band

def read_modeoutput(fmodeoutput, period_band_list):
    mode_output = fmodeoutput.strip()
    with open(mode_output) as fh:
        lines = fh.readlines()

    atten_coef_period = {}

    for period in sorted(period_band_list):
        fmin = period_band_list[period]["min"]
        fmax = period_band_list[period]["max"]

        atten_coef = []
        model_grpv = []

        for line in lines:
            if "0 R" in line and "skip" not in line:
                (mode_n, mode_type, mode_l, mode_freq, mode_period,
                 phsv, grpv, q_rayl, qaulity, dummy) = line.strip().split()

                #print("%s %f" % (grpv, float(grpv)))
                atten_coef_temp = \
                    np.pi * float(mode_freq) / (float(grpv) * float(q_rayl))

                if fmin <= float(mode_freq) <= fmax:
                    atten_coef.append(atten_coef_temp)
                    model_grpv.append(float(grpv))
                    #print period, mode_freq, grpv, q_rayl, fmin, fmax
                    #print len(atten_coef), len(model_grpv)

        atten_coef_mean = np.mean(atten_coef)
        atten_coef_stdv = np.std(atten_coef)

        model_grpv_mean = np.mean(model_grpv)
        model_grpv_stdv = np.std(model_grpv)

        #print("%s %f %f" % (period, atten_coef_mean, model_grpv_mean))

        atten_coef_period[period] = atten_coef_mean

    return atten_coef_period

def read_param(ftemplate):
    with open(ftemplate, "r") as fh:
        lines = fh.readlines()
    return lines

def read_index(f_index):
    with open(f_index, "r") as fh:
        lines = fh.readlines()

    index_dict = {}
    for line in lines:
        (index, key) = line.strip().split()
        index_dict[key] = int(index)
    return index_dict

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

    # calculate data matrix
    # can be move out to config later
    fraw_data = "/home/youyir/Cascadia/ProjectQ/data/atten_raw.dat"
    fraw_data_json = "/home/youyir/Cascadia/ProjectQ/data/atten.json"
    raw_data = _read_text(fraw_data)
    data = {}
    for line in raw_data:
        (period_str, atten_coef, atten_stdv) = line.strip().split()
        period_key = period_str.split(".")[0]
        #print period_key
        item = {}
        item["data"] = float(atten_coef)
        item["stdv"] = float(atten_stdv)
        data[period_key] = item


    # calculate atten_coef at each period for perturbed models
    atten_coef_matrix = {}
    for model in mode_outputs:
        (dummy, layer_id) =  model.rstrip().split("_")
        #print ("%s" % layer_id)
        atten_coef_period = read_modeoutput(model, period_list)
        atten_coef_matrix[layer_id] = atten_coef_period

    # calculate atten_coef at each period for reference model
    model = args.reference_mode_output
    (dummy, layer_id) =  model.rstrip().split("_")
    #print ("%s" % layer_id)
    atten_coef_period_ref = read_modeoutput(model, period_list)

    # read index for formated output
    m_index = read_index(model_index)
    d_index = read_index(data_index)

    # fixed perturbation
    delta_q = 50.

    # calculate gradient
    f_kernel = open("kernel_matrix.inp","w")
    gradient = {}
    for layer_id in sorted(atten_coef_matrix):
        synt = atten_coef_matrix[layer_id]
        gradient_period = {}
        for period in sorted(atten_coef_period_ref):
            atten_p = synt[period]
            atten_r = atten_coef_period_ref[period]
            atten_stdv = data[period]["stdv"]
            g = (atten_p - atten_r) / delta_q / atten_stdv
            gradient_period[period] = g
            print ("%s %s %e %e %e" % (layer_id, period, atten_p, atten_r, g))
            f_kernel.write("%3d %3d %e\n" %
                           (d_index[period], m_index[layer_id], g))
        gradient[layer_id] = gradient_period

    f_kernel.close()

    with open(args.gradient_output, 'w') as f:
        json.dump(gradient, f, indent=2, sort_keys=True)

    f_synt = open("atten_coef_syn.dat", "w")
    f_data = open("data_matrix.inp", "w")

    d_matrix = {}

    for period in sorted(atten_coef_period_ref):
        atten_synt = atten_coef_period_ref[period]
        atten_data = data[period]["data"]
        atten_stdv = data[period]["stdv"]
        delta_d = (atten_data - atten_synt) / atten_stdv
        d_matrix[period] = delta_d
        f_data.write("%3d %e\n" % (d_index[period], delta_d))
        f_synt.write("%s %13e %13e %13e %13e\n" % (period, atten_synt,
                                                   atten_data, atten_stdv,
                                                   delta_d))

    f_data.close()

    with open(fraw_data_json, 'w') as f:
        json.dump(data, f, indent=2, sort_keys=True)

    with open("atten_coef_syn.json", "w") as f:
        json.dump(atten_coef_period_ref, f, indent=2, sort_keys=True)

    with open("data_matrix.json", "w") as f:
        json.dump(d_matrix, f, indent=2, sort_keys=True)



