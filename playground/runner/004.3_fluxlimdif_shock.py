#!/usr/bin/env python

from fortdriver.core import Context, Setup
from fortdriver.units import PhysicalConstants
import os
import math

def defaultSetup():
    setup = Setup()
    setup.xmin          = -1e15
    setup.xmax          =  1e15
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
    setup.resolution    = 512
    return setup

def main():
    pc = PhysicalConstants()
    # unit = pc.Scale()
    # stss = [8,16,32,64,128,256,512,1024]
    # stss = [16,32,64,128,256,512,1024]
    # for steps in stss:
    # steps = 512
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
    place.setup.tfinish         = 1e9
    place.setup.dtprint         = place.setup.tfinish/1e2
    place.setup.resultfile      = resfile+".info"
    place.setup.process         = place.epc['borderless']
    place.setup.time            = 0.0
    place.setup.disotropic      = 1.0
    place.setup.eqs             = place.eeq['manual']
    place.setup.eqonhydro       = place.eif['yes']
    # place.setup.ddw             = place.esd['2nw_ds']
    place.setup.ddw             = place.esd['fab']
    place.setup.eqondiff        = place.eif['yes']
    place.setup.eqonfluxlim     = place.eif['yes']
    # place.setup.eqonfluxlim     = place.eif['no']
    # place.setup.eqondiff        = place.eif['no']
    place.setup.eqonradexch     = place.eif['yes']
    place.setup.eqonsts         = place.eif['yes']
    place.setup.ststype         = place.ests['auto']
    # place.setup.ststype         = place.ests['fixeds']
    # place.setup.stsfixeds       = steps

    def c(i):
        return lambda x: i

    rho = 1e-10
    mass = place.CalcUniformMass(rho)
    h = place.setup.hfac * place.setup.spacing
    csound = 3.2e5
    v0 = csound
    kappa = 4e-3
    mu = 2.0
    gamma = 5./3.

    place.setup.dcondconst    = kappa
    place.setup.molecularmass = mu
    place.setup.gamma         = gamma
    # E = ksi * rho
    # E = 4*sigmab*T_rad^4/lightspeed
    # ksi = 4*sigmab*T_rad^4/lightspeed/rho
    # ksi = 4.0*pc.STEBOLTZ*((gamma-1)*mu*u/R_g)**4.0/pc.C/rho
    T = 1500.0
    # E = 4*pc.STEBOLTZ*T**4/pc.C
    ksi = 4.0*pc.STEBOLTZ*T**4.0/pc.C/rho
    # T_gas = (gamma - 1) * mu * u / R_g
    # T_rad = (E*lightspeed/4/sigmab)^(1/4)
    # T_rad = (ksi*rho*lightspeed/4/sigmab)^(1/4)
    u = T*pc.RG/(gamma - 1.0)/mu
    prs = u*rho*(gamma - 1.0)
    print('P = ', prs)
    print('T_gas = ', (gamma - 1.0) * mu * u / pc.RG)
    print('T_rad = ', (ksi*rho*pc.C/4.0/pc.STEBOLTZ)**(1./4.))

    place.ModifyParticlesProperties(
        properties  = ['kappa', 't',    'c',      'm',   'den',  'h',   'p',   'u',  'om'],
        value       = [c(kappa), c(ksi), c(csound), c(mass), c(rho), c(h), c(prs), c(u), c(1.0)]
    )

    place.AddParticles(
        dim='x',
        type='fixedreal',
        method='periodic'
    )

    speedFunc = lambda rx: v0 if (rx <= 0.0) else -v0
    place.ModifyParticlesProperties(
        properties  = ['vx'],
        value       = [speedFunc],
        valuearg    = ['rx']
    )

    # place.PrintState()
    place.Apply()
    # place.ContinueRun()
    # place.BackupOutput("output_shock_m3_fab_"+'{0:0>4}'.format(steps))
    # print("Done: " + str(steps))
main()
