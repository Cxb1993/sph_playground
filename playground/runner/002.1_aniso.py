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
    setup.ics  = "shock12"
    setup.process   = "relaxation"
    setup.tfinish   = 20.0
    setup.equations = "hydro"
    setup.kernel    = "m6"
    setup.hfac      = 1.0
    setup.coordsys  = "cartesian"
    setup.silent    = "no"
    setup.usedumps  = "yes"
    setup.adden     = "yes"
    setup.artts     = "yes"
    setup.au        = 1.0
    setup.nsnapshots = 10.0
    setup.resultfile = "result.info"
    setup.influencefile = "influence.info"
    return setup

def main():
    # Read hcready dumpfiles after 002_anisopy and change derivatives type
    # Same particle placement but different diffusion
    ddwt = "2nw_sd"
    hc12 = Context()
    hc12.ReadDumpFile("output/dm.hcready", "output/fd.hcready")
    hc12.setup.ddw = hc12.esd[ddwt]
    hc12.Apply()
    hc12.ContinueRun()
    hc12.BackupOutput("output-002-1")
    # print("Done: | ", ddwt, " | ", outrest , "|",
    #     hc12.PrintCell(hc12.setup.resultfile, -1, 3))
    print("Error: ", hc12.PrintCell(hc12.setup.resultfile, -1, 3))

main()
