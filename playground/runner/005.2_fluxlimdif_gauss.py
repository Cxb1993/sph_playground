#!/usr/bin/env python

from fortdriver.core import Context, Setup
from fortdriver.units import PhysicalConstants
import os
import math
import time

def defaultSetup():
    setup = Setup()
    setup.xmin          = -1
    setup.xmax          =  1
    setup.dim           = 1.0
    setup.placement     = "uniform"
    setup.tfinish       = 0.0
    setup.kernel        = "m6"
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
    setup.resolution    = 1024
    return setup

def main():
    pc = PhysicalConstants()
    # unit = pc.Scale()
    stss = [8,16,32,64,128,256,512,1024]
    for steps in stss:
        # steps = stss[0]
        place = Context()
        resfile = place.GetDateTimeString()
        place.setup = defaultSetup()
        place.CleanAllLogs()
        # place.SetThreadsOMP(2)
        # place.SimpleMake(debug=True)
        place.SimpleMake()
        place.CleanRun()
        place.BackupDump("output/dump_full.h5")
        place.ReadDumpHDF("output/dump_full.h5")
        # place.PrintState()
        place.setup.tfinish         = 2e-12
        place.setup.dtprint         = place.setup.tfinish/1e2
        place.setup.resultfile      = resfile+".info"
        place.setup.process         = place.epc['borderless']
        place.setup.time            = 0.0
        place.setup.disotropic      = 1.0
        place.setup.eqs             = place.eeq['manual']
        place.setup.ddw             = place.esd['2nw_ds']
        place.setup.eqonhydro       = place.eif['no']
        place.setup.eqondiff        = place.eif['yes']
        place.setup.eqonfluxlim     = place.eif['yes']
        place.setup.eqonradexch     = place.eif['no']
        place.setup.eqonsts         = place.eif['yes']
        place.setup.ststype         = place.ests['fixeds']
        place.setup.stsfixeds       = steps

        def c(i):
            return lambda x: i

        rho = 1
        mass = place.CalcUniformMass(rho)
        h = place.setup.hfac * place.setup.spacing
        csound = 1
        kappa = 1
        mu = 2.0
        gamma = 5./3.

        # 1./sqrt(4.*pi*(29979245800./3./1.*t+0.00125))*exp(-x*x/4./(29979245800./3./1.*t+0.00125))

        place.setup.dcondconst    = kappa
        place.setup.molecularmass = mu
        place.setup.gamma         = gamma

        sigma = 1./20.
        s2 = sigma*sigma
        def ksi(x):
            return 1./math.sqrt(2.*math.pi*s2)*math.exp(-x*x/2./s2)
        def u(x):
            return (ksi(x)*pc.C*rho/4./pc.STEBOLTZ)**(1./4.)*pc.RG/(gamma - 1.)/mu
        def prs(x):
            return u(x)*rho*(gamma - 1.0)

        place.AddParticles(
            dim='x',
            type='fixedreal',
            method='periodic'
        )

        place.ModifyParticlesProperties(
            properties  = ['kappa', 't',    'c',      'm',   'den',  'h',   'p',   'u',  'om'],
            value       = [c(kappa), ksi, c(csound), c(mass), c(rho), c(h),  prs,   u,  c(1.0)]
        )
        place.Apply()
        place.ContinueRun()
        place.BackupOutput("output_gauss_1024_"+'{0:0>4}'.format(steps))
        print("Done: " + str(steps))
main()
