# TODO:
# - add instruction for printing on .csv

import random
import numpy as np
import sys
import getopt
import argparse

from producecfg import *
from libSOGA import *
from sogaPreprocessor import transform2gm

from time import time

random.seed(0)
np.random.seed(0)

def getCliCmd():
    # Create ArgumentParser object
    parser = argparse.ArgumentParser(description="SOGA CLI")

    parser.add_argument("-f","--modelfile", help="SOGA model path",required=True)
    parser.add_argument("-o","--outputfile", help="Output file path",required=False)

    # Add optional flag
    parser.add_argument("-c", "--covariance", action="store_true", help="Output covariance",required=False)

    # Add list of strings
    parser.add_argument("-v","--vars", nargs="*", default=[],help="List of output variables",required=False)

    # Parse the command-line arguments
    args = parser.parse_args()

    return args

def printOutput(output_dist,preprocTime,cfgTime,sogaTime,args):
    var_list = []
    var_idx = []

    if args.outputfile is not None:
        sys.stdout=open(args.outputfile, 'w')

    mwidth=np.max([len(preprocTime),len(cfgTime),len(sogaTime)])

    print(f'SOGA preprocessing in: {preprocTime.rjust(mwidth)} s')
    print(f'      CFG produced in: {cfgTime.rjust(mwidth)} s')
    print(f'              Runtime: {sogaTime.rjust(mwidth)} s')
    print(f"c: {(output_dist.gm.n_comp()):d}")
    print(f"d: {(len(output_dist.var_list)):d}")

    if len(args.vars) == 0:
        for var, val in zip(output_dist.var_list, output_dist.gm.mean()):
            print('E['+var+']:', round(val,5))
    else:
        for var in args.vars:
            i = output_dist.var_list.index(var)
            print('E['+var+']:', round(output_dist.gm.mean()[i],5))
            var_idx.append(i)
    print('\n')

    if args.covariance:
        if len(args.vars) == 0:
            print('Covariance:\n', np.around(output_dist.gm.cov(),5))
        else:
            print('Covariance:\n', np.around(output_dist.gm.cov()[var_idx][:,var_idx], 5))

def printBanner():
    print('/ ___| / _ \ / ___|  / \\\n\___ \| | | | |  _  / _ \\\n ___) | |_| | |_| |/ ___ \\\n|____/ \___/ \____/_/   \_\\\n')


def SOGA():
    printBanner()

    args=getCliCmd()
    preproc_strt=time()
    compiledFile=transform2gm(args.modelfile)
    preproc_end=time()
    
    cfg_start = time()
    cfg = produce_cfg(compiledFile)
    cfg_end = time()

    comp_start = time()
    output_dist = soga_analyze(cfg)
    comp_end = time()

    preprocTime=f"{preproc_end-preproc_strt:<.3f}"
    cfgTime=f"{cfg_end-cfg_start:<.3f}"
    sogaTime=f"{comp_end-comp_start:<.3f}"

    printOutput(output_dist=output_dist,preprocTime=preprocTime
        ,cfgTime=cfgTime,sogaTime=sogaTime,args=args)


if __name__ == '__main__':
    SOGA()    