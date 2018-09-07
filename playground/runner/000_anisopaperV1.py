#!/usr/bin/env python

# in rev1 comment we asked to see how anisotropic diffusion goes on glass like lattice
# first I run hydro code with periodic conditions with uniform 'T'
# + then I
# check if all the particles are inside of the domain
# + and make
# x-borders to be fixed
# T = 1 for x < 0, T = 2 for x >= 0

from runnerbase import Context

def main():
    # for ddwt in ["2nw-+", "2nw+-"]
    problem = Context()
    resfile = problem.GetDateTimeString()
    # rightPlacing = False
    # it = 0
    # while(not rightPlacing):
    #     problem = Context()
    #     problem.setup.xmin = -1
    #     problem.setup.xmax =  1
    #     problem.setup.resolution = 16
    #     problem.setup.dimentions = 2
    #     problem.setup.initialConditions = "shock12"
    #     problem.setup.finishTime = 5
    #     problem.setup.equationsSet = "hydro"
    #     problem.setup.kernel = "m6"
    #     problem.setup.hfac = 1.0
    #     problem.setup.snapshotsNumber = 100
    #     problem.setup.resultFile = "result.info"
    #     # # problem.setup.particlesPlacement  = "allPeriodicXReal"
    #     # # problem.setup.particlesArrangment = "random"
    #     # # problem.setup.initialConditions   = "relaxation"
    #     # # problem.CleanSTDLogs()
    #     problem.SetThreadsOMP(4)
    #     if (it == 0):
    #         problem.SimpleMake()
    #     problem.CleanRun()
    #     problem.ReadFinalDump()
    #     problem.AfterRunUpdate()
    #     rightPlacing = problem.IsParticlesPositionsCorrect()
    #     it += 1
    #     rightPlacing = True
    # print("Spent iterations befor right placing: " + str(it))
    problem.ReadFinalDump()
    problem.AfterRunUpdate()
    problem.PrintState()
    # problem.BackupOutput()
    problem.setup.eqs = problem.eeq['diffusion']
    problem.setup.tfinish = 0.05
    problem.setup.ddw = problem.esd['2nw_sd']
    problem.setup.resultfile = resfile + "2nw-+.info"
    problem.PrintState()
main()
