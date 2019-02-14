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
    pc = PhysicalConstants()
    # unit = pc.Scale()
    place = Context()
    resfile = place.GetDateTimeString()
    place.setup = defaultSetup()
    place.CleanAllLogs()
    # place.SetThreadsOMP(2)
    # place.SimpleMake(debug=True)
    place.SimpleMake()
    place.CleanRun()
    place.ReadDumpHDF("./output/dump_full.h5")
    # place.PrintState()
    place.setup.tfinish         = 1e9
    place.setup.dtprint         = place.setup.tfinish/1e2
    place.setup.resultfile      = resfile+".info"
    place.setup.process         = place.epc['borderless']
    place.setup.time            = 0.0
    place.setup.gamma           = 5./3.
    place.setup.disotropic      = 1.0
    place.setup.eqs             = place.eeq['manual']
    place.setup.eqonhydro       = place.eif['yes']
    # place.setup.eqondiff        = place.eif['yes']
    place.setup.eqondiff        = place.eif['no']
    place.setup.ddw             = place.esd['2nw_ds']
    # place.setup.eqonfluxlim     = place.eif['yes']
    place.setup.eqonfluxlim     = place.eif['no']

    def c(i):
        return lambda x: i

    rho = 1e-10
    mass = place.CalcUniformMass(rho)
    h = place.setup.hfac * place.setup.spacing
    csound = 3.2e5
    v0 = csound
    kappa = 4e1
    place.setup.dcondconst = kappa

    # E = ksi * rho
    # E = 4*sigmab*T_rad^4/lightspeed
    T = 1500
    E = 4*pc.STEBOLTZ*T**4/pc.C
    # T_gass = (gamma - 1) * mu * u / R_g
    mu = 2.1
    u = T*pc.RG/(place.setup.gamma - 1)/mu
    prs = u*rho*(place.setup.gamma - 1)
    print('P =',prs)

    # addedN = place.AddParticles(
    #     dim='x',
    #     type='fixed',
    #     method='periodic'
    # )
    # place.setup.fixedpn = addedN

    place.ModifyParticlesProperties(
        properties  = ['kappa', 't',    'c',      'm',   'den',  'h',   'p',   'u',  'om'],
        value       = [c(kappa), c(E), c(csound), c(mass), c(rho), c(h), c(prs), c(u), c(1.0)]
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
main()
