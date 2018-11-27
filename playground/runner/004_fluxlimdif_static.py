#!/usr/bin/env python

from runnerbase import Context, Setup
import os
import math

def defaultSetup():
    setup = Setup()
    setup.xmin          = -1.0
    setup.xmax          =  1.0
    setup.dim           = 1.0
    setup.placement     = "uniform"
    setup.tfinish       = 0.0
    setup.kernel        = "m4"
    setup.hfac          = 1.0
    setup.nsnapshots    = 10.0
    setup.resultfile    = "result.info"
    setup.influencefile = "influence.info"
    setup.coordsys      = "cartesian"
    setup.silent        = "no"
    setup.usedumps      = "yes"
    setup.adden         = "yes"
    setup.artts         = "yes"
    setup.au            = 1.0
    setup.resolution    = 128
    return setup

def main():
    place = Context()
    resfile = place.GetDateTimeString()
    place.setup = defaultSetup()
    place.CleanAllLogs()
    place.SetThreadsOMP(4)
    place.SimpleMake()
    place.CleanRun()
    place.ReadFinalDump()
    place.PrintState()
    place.setup.tfinish         = 100.0
    place.setup.dtprint         = place.setup.tfinish/100.
    place.setup.resultfile      = resfile+".info"
    place.setup.process         = place.epc['fullyperiodic']
    place.setup.time            = 0.0
    place.setup.gamma           = 5./3.
    place.setup.disotropic      = 1.0
    place.setup.dcondconst      = 0.4
    place.setup.eqs             = place.eeq['manual']
    place.setup.eqonhydro       = place.eif['yes']
    place.setup.eqondiff        = place.eif['yes']
    place.setup.ddw             = place.esd['2nw_ds']
    place.setup.eqonfluxlim     = place.eif['yes']
    place.PrintState()

    def c(i):
        return lambda x: i

    rho = 1e-7
    mass = place.CalcUniformMass(rho)
    h = place.setup.hfac * place.setup.spacing
    e = 1e-2
    # e = u * rho => u = e / rho
    u = e/rho
    # u = prs/(place.setup.gamma -1)/rho
    prs = u*rho*(place.setup.gamma -1)
    # E = ksi * rho => ksi = E / rho
    E = 1e+12

    place.ModifyParticlesProperties(
        properties  = ['kappa', 't',    'c',        'm',   'den',  'h',  'p',    'u'],
        value       = [c(0.4), c(E), c(1.235e+8), c(mass), c(rho), c(h), c(prs), c(u)]
    )

    place.Apply()
    place.ContinueRun()
    
main()
