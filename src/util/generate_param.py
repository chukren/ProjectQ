#!/usr/bin/env python

ftemplate = "PAR_RAYLEIGH"
model_list = "model.list"

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

param = read_param(ftemplate)

models = read_modellist(model_list)

for model in models:
    (modelname, appendix) = model.split("_")
    mode_output = "%s_%s" % ("rayleigh.out", appendix)
    fparam = "%s_%s" % ("param", appendix.rstrip())
    generate_param(param, model, mode_output, fparam)
