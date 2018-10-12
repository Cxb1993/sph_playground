#!/usr/bin/env python

# ./execute
#x --dim 2
#x --equations diffusion
#x --ics pulse
# --resolution 64
# --ddw 2nw_sd
# --tfinish .05
# --hfac 1.
# --nsnapshots 10
# --coordsys cartesian
# --process backcompatibility
# --usedumps yes
# --influencefile influence.info
# --resultfile result.info
# --adden yes
# --artts yes
# --au -1

#!/usr/bin/env python

# in rev1 comment we asked to see how anisotropic diffusion goes on glass like lattice
# first I run hydro code with periodic conditions with uniform 'T'
# + then I
# check if all the particles are inside of the domain
# + and make
# x-borders to be fixed
# T = 1 for x < 0, T = 2 for x >= 0

from runnerbase import Context, Setup
import os
import math

def defaultSetup():
    setup   = Setup()
    setup.xmin = -1.0
    setup.xmax =  1.0
    setup.dim  =  2.0
    setup.ics  = "pulse"
    setup.process   = "backcompatibility"
    setup.tfinish   = 0.05
    setup.equations = "diffusion"
    setup.kernel    = "m6"
    setup.hfac      = 1.0
    setup.coordsys  = "cartesian"
    setup.silent    = "no"
    setup.usedumps  = "yes"
    setup.adden     = "yes"
    setup.nsnapshots = 10.0
    setup.resultfile = "result.info"
    setup.influencefile = "influence.info"
    setup.ddw        = "2nw_ds"
    return setup

def main():
    # Read hcready dumpfiles after 002_anisopy and change derivatives type
    # Same particle placement but different diffusion
    tmp = Context()
    resfile = tmp.GetDateTimeString()
    tmp.setup = defaultSetup()
    tmp.SetThreadsOMP(6)
    tmp.SimpleMake()
    # for rest in [16, 32, 64, 128, 256, 512, 1024]:
    for s_rest in [16, 32, 64]:
        for s_artts in ['yes', 'no']:
            pulse = Context()
            pulse.setup = defaultSetup()
            pulse.setup.resolution = s_rest
            if (s_artts == "no"):
                pulse.setup.artts = 'no'
                pulse.setup.au    = -1
            elif (s_artts == "yes"):
                pulse.setup.artts = 'yes'
                pulse.setup.au    = -1

            pulse.setup.resolution = s_rest
            rs = '{0:0>4}'.format(s_rest)
            ss = str(s_artts)
            pulse.setup.resultfile = resfile+'-'+ss+".info"
            pulse.CleanRun()
            pulse.BackupOutput("output-"+rs+"-"+ss)
            print(rs, '|', ss,'|', pulse.PrintCell(pulse.setup.resultfile, -1, 3))
        print('')

main()
