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
    tmp.SetThreadsOMP(6)
    tmp.SimpleMake()

    # for rest in [16, 32]:
    for rest in [16, 32, 64, 128, 256, 512, 1024]:
        relax = Context()
        relax.setup = defaultSetup()
        relax.setup.resolution = rest
        relax.setup.ddw        = "fab"
        relax.CleanRun()
        relax.ReadFinalDump()
        fd = "output/"+"fd-" + str(rest)+".relaxed"
        dm = "output/"+"dm-" + str(rest)+".relaxed"
        relax.BackupDumps(fd, dm)
        for ddwt in ["2nw_sd", "2nw_ds"]:
            for s_artts in ["yes", "no"]:
                hc12                    = relax.CopyContext()
                hc12.setup.ddw          = hc12.esd[ddwt]
                hc12.setup.eqs          = hc12.eeq['diffusion']
                hc12.setup.tfinish      = 0.05
                hc12.setup.dtprint      = hc12.setup.tfinish/10.
                hc12.setup.resultfile   = resfile+"-"+ddwt+'-'+s_artts+".info"
                hc12.setup.process      = hc12.epc['backcompatibility']
                hc12.setup.time         = 0.0
                hc12.setup.artts        = hc12.eif[s_artts]
                if (ddwt == "2nw_ds") or (ddwt == "2nw_sd"):
                    hc12.setup.au = -1.0
                else:
                    hc12.setup.au = 0.0

                addedN = hc12.AddParticles(
                    dim='x',
                    type='fixed',
                    method='periodic'
                )
                hc12.setup.fixedpn = addedN

                initialTemp = lambda rx: 1.0 if (rx <= 0.0) else 2.0
                if (hc12.setup.artts == hc12.eif["no"]):
                    h = 2./relax.setup.resolution
                    hc12.setup.time = h*h/4.
                    initialTemp = lambda rx: 1.5 + 0.5 * math.erf(rx/h)
                else:
                    initialTemp = lambda rx: 1.0 if (rx <= 0.0) else 2.0

                hc12.ModifyParticles(
                    condition   = lambda rx: True,
                    condarg     = 'rx',
                    properties  = ['u', 't'],
                    value       = initialTemp,
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
                ds = str(ddwt)
                rs = str(rest)
                ss = str(s_artts)
                relax.BackupOutput("output-"+ds+"-"+rs+"-"+ss)
                print("Done: | ", ds, " | ", rs, "|", ss, "|", hc12.PrintCell(hc12.setup.resultfile, -1, 3))

main()
