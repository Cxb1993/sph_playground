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
    setup.ics = "pulse"
    setup.process   = "relaxation"
    setup.tfinish   = 100.0
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

    # for rest in [256]:
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
        # for ddwt in ["2nw_ds"]:
        for ddwt in ["2nw_ds", "2nw_ds", "n2w", "fab", "fw"]:
            if ddwt == "2nw_ds":
                smlist = ["no", "artterm"]
            else:
                smlist = ["no"]
            for s_artts in smlist:
                diff                    = relax.CopyContext()
                diff.setup.ddw          = diff.esd[ddwt]
                diff.setup.eqs          = diff.eeq['diffusion']
                diff.setup.tfinish      = 0.025
                diff.setup.dtprint      = diff.setup.tfinish/10.
                diff.setup.resultfile   = resfile+"-"+ddwt+'-'+s_artts+".info"
                diff.setup.process      = diff.epc['backcompatibility']
                diff.setup.time         = 0.0
                diff.setup.disotropic   = diff.eif['no']


                addedN = diff.AddParticles(
                    dim='x',
                    type='fixed',
                    method='periodic'
                )
                diff.setup.fixedpn = addedN

                initialTemp = lambda rx, ry, rz: ((2.*math.pi)**(-setup.dim/2.))/((0.1**2)**(setup.dim/2.))*math.exp(-0.5*(rx*rx + ry*ry + rz*rz)/(0.1**2))

                if (s_artts == "no"):
                    diff.setup.artts = diff.eif['no']
                elif (s_artts == "artterm"):
                    diff.setup.artts = diff.eif['yes']
                    diff.setup.au    = -1

                diff.ModifyParticles(
                    properties  = ['u', 't'],
                    value       = initialTemp,
                    valuearg    = ['rx', 'ry', 'rz']
                )
                diff.ModifyParticles(
                    properties  = [
                        'dtdx', 'dtdy', 'dtdz',
                        'vx', 'vy', 'vz',
                        'ax', 'ay', 'az',
                        'by, bz'],
                    value       = lambda rx: 0.0,
                )
                hc12.ModifyParticles(
                    properties  = ['bx'],
                    value       = lambda rx: 1.0,
                )
                # diff.PrintState()
                diff.Apply()
                ds = str(ddwt)
                rs = '{0:0>4}'.format(rest)
                ss = str(s_artts)
                diff.BackupDumps(
                    "output/"+"fd-"+ds+"-"+rs+"-"+ss+".hcready",
                    "output/"+"dm-"+ds+"-"+rs+"-"+ss+".hcready"
                )
                diff.ContinueRun()
                relax.BackupOutput("output-"+ds+"-"+rs+"-"+ss)
                print("Done: | ", ds, " | ", rs, "|", ss, "|", diff.PrintCell(diff.setup.resultfile, -1, 3))
        print('\n')
main()
