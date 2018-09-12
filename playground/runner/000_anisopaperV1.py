#!/usr/bin/env python

# in rev1 comment we asked to see how anisotropic diffusion goes on glass like lattice
# first I run hydro code with periodic conditions with uniform 'T'
# + then I
# check if all the particles are inside of the domain
# + and make
# x-borders to be fixed
# T = 1 for x < 0, T = 2 for x >= 0

from runnerbase import Context, Setup

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
    return setup

def main():
    tmp = Context()
    resfile = tmp.GetDateTimeString()
    tmp.setup = defaultSetup()
    tmp.SimpleMake()

    # for ddwt in ["2nw_sd", "2nw_ds", "n2w", "fab", "fw"]:
    for ddwt in ["2nw_sd", "2nw_ds"]:
        # for rest in [16, 32, 64, 128, 256, 512, 1024]:
        for rest in [16]:
            relax = Context()
            relax.SetThreadsOMP(8)
            relax.setup = defaultSetup()
            relax.setup.ddw        = ddwt
            relax.setup.resolution = rest
            relax.CleanRun()
            relax.ReadFinalDump()
            relax.BackupDumps(
                "output/"+"fd-" + str(ddwt) + "-" + str(rest)+".relaxed",
                "output/"+"dm-" + str(ddwt) + "-" + str(rest)+".relaxed"
            )
            hc12 = relax.CopyContext()
            hc12.setup.eqs          = hc12.eeq['diffusion']
            hc12.setup.tfinish      = 0.05
            hc12.setup.dtprint      = hc12.setup.tfinish/10.
            hc12.setup.resultfile   = resfile + "-" + ddwt + ".info"
            hc12.setup.process      = hc12.epc['backcompatibility']
            hc12.setup.time = 0.0
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
                properties  = ['dtdx', 'dtdy', 'dtdz'],
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
