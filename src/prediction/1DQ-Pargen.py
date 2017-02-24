#!/usr/bin/env python
import os, sys
import numpy as np
import argparse

agelist = {"1MA": 1.0,
           "3MA": 3.0,
           "5MA": 5.0,
           "7MA": 7.0,
           "9MA": 9.0,
           "11MA": 11.0}

periodlist = {"50s": 0.02,
              "20s": 0.05}

def make_input_file(agelist, periodlist, solidus):
    for age_k, age_v in agelist.iteritems():
        for period_k, period_v in periodlist.iteritems():
            filelabel = "%s.%s.%s" % (age_k, period_k, solidus)

            predir = os.path.join(filelabel)
            if not os.path.exists(predir):
                os.makedirs(predir)

            filename = os.path.join(predir, "PAR.%s.YT2016" % filelabel)

            with open(filename, "w") as f:
                f.write("%-6.2f\n" % age_v)
                f.write("%-6.2f\n" % period_v)
                f.write("yamauchi.%s.dat\n" % filelabel)
                f.write("solidus.yamauchi.%s.dat\n" % filelabel)
                f.write("modulus.yamauchi.%s.dat\n" % filelabel)
                f.write("background.yamauchi.%s.dat\n" % filelabel)
            print("%s" % predir)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action='store', dest='solidus',
                        required=True)
    args = parser.parse_args()
    solidus = args.solidus
    make_input_file(agelist, periodlist, solidus)
