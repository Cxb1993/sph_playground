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
    tmp.SetThreadsOMP(8)
    tmp.SimpleMake()

    # for ddwt in ["2nw_sd"]:
        # for rest in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0]:
        #     for artts in ["yes"]:
    for rest in [16, 32, 64, 128, 256, 512, 1024]:
        for ddwt in ["2nw_sd", "2nw_ds", "n2w", "fab", "fw"]:
            relax = Context()
            relax.setup = defaultSetup()
            relax.setup.resolution = rest
            relax.setup.ddw        = ddwt
            relax.CleanRun()
            relax.ReadFinalDump()
            fd = "output/"+"fd-" + str(ddwt) + "-" + str(rest)+".relaxed"
            dm = "output/"+"dm-" + str(ddwt) + "-" + str(rest)+".relaxed"
            relax.BackupDumps(fd, dm)
            hc12                    = relax.CopyContext()
            hc12.setup.ddw          = hc12.esd[ddwt]
            hc12.setup.eqs          = hc12.eeq['diffusion']
            hc12.setup.tfinish      = 0.05
            hc12.setup.dtprint      = hc12.setup.tfinish/10.
            hc12.setup.resultfile   = resfile + "-" + ddwt + ".info"
            hc12.setup.process      = hc12.epc['backcompatibility']
            hc12.setup.time         = 0.0
            hc12.setup.au           = 10.0
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
                "output/"+"fd-" + str(ddwt) + "-" + str(rest)+".hcready",
                "output/"+"dm-" + str(ddwt) + "-" + str(rest)+".hcready"
            )
            hc12.ContinueRun()
            relax.BackupOutput("output-"+str(ddwt) + "-" + str(rest))
            print("Done: | ", ddwt, " | ", rest, "|")

main()
