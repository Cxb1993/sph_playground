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

def defaultSetup():
    setup   = Setup()
    setup.xmin = -1.0
    setup.xmax =  1.0
    setup.dim = 2.0
    setup.ics = "shock12"
    setup.process   = "relaxation"
    setup.tfinish   = 20.0
    setup.equations = "hydro"
    setup.kernel    = "m6"
    setup.hfac      = 1.0
    setup.nsnapshots = 10.0
    setup.resultfile = "result.info"
    setup.influencefile = "influence.info"
    setup.coordsys  = "cartesian"
    setup.silent    = "no"
    setup.usedumps  = "yes"
    setup.adden     = "yes"
    setup.artts     = "yes"
    setup.au        = 1.0
    return setup

def main():
    tmp = Context()
    resfile = tmp.GetDateTimeString()
    tmp.setup = defaultSetup()
    tmp.SetThreadsOMP(4)
    tmp.SimpleMake()
    ddwt = "2nw_sd"
    relax = Context()
    relax.setup = defaultSetup()
    relax.setup.resolution = 64
    relax.setup.ddw        = ddwt
    relax.setup.resultfile = resfile + "-" + ddwt + ".info"
    relax.CleanRun()
    relax.ReadFinalDump()
    fd = "output/"+"fd-" + str(ddwt)+".relaxed"
    dm = "output/"+"dm-" + str(ddwt)+".relaxed"
    relax.BackupDumps(fd, dm)
    hc12                    = relax.CopyContext()
    hc12.setup.eqs          = hc12.eeq['diffusion']
    hc12.setup.tfinish      = 0.05
    hc12.setup.dtprint      = hc12.setup.tfinish/10.
    hc12.setup.process      = hc12.epc['backcompatibility']
    hc12.setup.time         = 0.0
    hc12.setup.au           = -1
    addedN = hc12.AddParticles(
        dim='x',
        type='fixed',
        method='periodic'
    )
    hc12.setup.fixedpn = addedN
    hc12.ModifyParticles(
        condition   = lambda rx: True,
        condarg     = 'rx',
        properties  = ['u', 't'],
        value       = lambda rx: 1.0 if (rx <= 0.0) else 2.0,
        valuearg    = 'rx'
    )
    hc12.ModifyParticles(
        condition   = lambda rx: True,
        condarg     = 'rx',
        properties  = [
            'dtdx', 'dtdy', 'dtdz',
            'vx', 'vy', 'vz',
            'ax', 'ay', 'az'],
        value       = lambda rx: 0.0,
        valuearg    = 'rx'
    )
    # hc12.PrintState()
    hc12.Apply()
    hc12.BackupDumps(
        "output/"+"fd.hcready",
        "output/"+"dm.hcready"
    )
    hc12.ContinueRun()
    relax.BackupOutput("output-"+ '{0:0^4}'.format(relax.setup.resolution))
    # print("Done: | ", ddwt, " | ", outrest , "|",
    #     hc12.PrintCell(hc12.setup.resultfile, -1, 3))
    print(hc12.PrintCell(hc12.setup.resultfile, -1, 3))

main()
