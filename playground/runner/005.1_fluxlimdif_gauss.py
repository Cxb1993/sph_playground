#!/usr/bin/env python

from fortdriver.core import Context, Setup
from fortdriver.units import PhysicalConstants
import os
import math

def defaultSetup():
    setup = Setup()
    setup.xmin          = -5
    setup.xmax          =  5
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
    place.BackupDump("output/dump_full.h5")
    place.ReadDumpHDF("output/dump_full.h5")
    # place.PrintState()
    place.setup.tfinish         = 4e-10
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

    def c(i):
        return lambda x: i

    rho = 1
    mass = place.CalcUniformMass(rho)
    h = place.setup.hfac * place.setup.spacing
    csound = 1
    kappa = 0.1
    mu = 2.0
    gamma = 5./3.

    place.setup.dcondconst    = kappa
    place.setup.molecularmass = mu
    place.setup.gamma         = gamma

    sigma = 1./3.
    s2 = sigma*sigma
    def ksi(x):
        return 1./math.sqrt(2.*math.pi*s2)*math.exp(-x*x/2./s2)
    def u(x):
        return (ksi(x)*pc.C*rho/4./pc.STEBOLTZ)**(1./4.)*pc.RG/(gamma - 1.)/mu
    def prs(x):
        return u(x)*rho*(gamma - 1.0)

    print('P = ', prs(0))
    print('T_gas = ', (gamma - 1.0) * mu * u(0) / pc.RG)
    print('T_rad = ', (ksi(0)*rho*pc.C/4.0/pc.STEBOLTZ)**(1./4.))

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
main()
