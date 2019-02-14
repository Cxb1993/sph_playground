## /usr/bin/env python3

import math

class PhysicalConstants(object):
    PI =  3.1415926536e0
    C = 2.997924e10                     # Speed of light            cm/s
    GG = 6.672041e-8                    # Gravitational constant    dyn cm^2 g^-2
    RG = 8.314e7                        # Gas constant              ergs mole^-1 K^-1
    CGSMU0 = 4.*PI
    ELECTRON_MASS_CGS = 9.10938291e-28  # Electron mass             g
    PROTON_MASS_CGS = 1.67262158e-24    # Proton mass               g
    ATOMIC_MASS_UBIT = 1.660538921e-24  # Atomic mass unit          g
    CROSS_SECTION_H2_CGS = 2.367e-15    # Hydrogen molecule cs      cm^-2
    RADCONST = 7.5646e-15               # Radiation constant        erg cm^-3 K^-4
    KBOLTZ = 1.38066e-16                # Boltzmann constant        erg/K
    KB_ON_MH = KBOLTZ/PROTON_MASS_CGS   # kB/m_H                    erg/K/g
    EV     = 1.60219e-12                # electron volt             erg
    QE     = 4.8032068e-10              # charge on electron        esu
    PLANCKH  =   6.6260755e-27          # Planck's Constant         erg.s
    PLANCKHBAR = 1.05457266e-27         # Planck's Constant/(2pi)   erg.s
    THOMCS     = 6.6525e-25             # Thomson cross section     cm^2
    FINESTR    = 7.2974e-3              # Fine structure constant   unitless
    STEBOLTZ   = 5.67051e-5             # Stefan-Boltzmann constant erg cm^-2K^-4 s^-1
    AVOGADRO   = 6.0221408577e23        # Avogadro's number         mole^-1
    RO         = 3.00000000             # Rossby number without dimension
    # Astronomical constants (cgs units)

    # Solar mass and radius
    SOLARM = 1.9891e33                  # Mass of the Sun           g
    SOLARR = 6.959500e10                # Radius of the Sun         cm
    SOLARL = 3.9e33                     # Luminosity of the Sun     erg/s

    # Earth mass and radius
    EARTHM = 5.979e27                   # Mass of the Earth         g
    EARTHR = 6.371315e8                 # Radius of the Earth       cm
    JUPITERM = 1.89813e30               # Mass of Jupiter           g
    CERESM = 8.958e23                   # Mass of Ceres             g

    # Distance scale
    AU = 1.496e13                       # Astronomical unit         cm
    LY = 9.4605e17                      # Light year                cm
    PC = 3.086e18                       # Parsec                    cm
    KPC = 3.086e21                      # Kiloparsec                cm
    MPC = 3.086e24                      # Megaparsec                cm
    KM = 1.0e5                          # Kilometer                 cm

    # Time scale
    SECONDS = 1.0e0
    MINUTES = 6.0e1
    HOURS = 3.6e3
    DAYS = 8.64e4
    YEARS = 3.1556926e7

    # Energy conversion
    EVTOK = 1.1604519e4                 # Degrees kelvin per eV     K/eV

    def Scale(self, idist=None, imass=None, itime=None, ig=None, ic=None):
        dist = idist if (idist != None) else 1.0
        mass = imass if (imass != None) else 1.0
        time = itime if (itime != None) else 1.0

        if (ic != None):
            if (imass != None):
                dist = self.GG*mass/self.C/self.C
                time = dist/self.C
            elif (idist != None):
                time = dist/self.C
                mass = self.C*self.C*dist/self.GG
            elif (itime != None):
                dist = time*self.C
                mass = self.C*self.C*dist/self.GG
            else:
                dist = self.GG*mass/self.C**2
                time = dist/self.C
        elif (ig != None):
            if (imass != None) and (idist != None):
                time = math.sqrt(dist**3/(self.GG*mass))
            elif (idist != None) and (itime != None):
                mass = dist**2/(self.GG*time**2)
            elif (imass != None) and (itime != None):
                dist = (time**2*(self.GG*mass))**(1.0/3.0)
            elif (itime != None):
                mass = dist**2/(self.GG*time**2)     # dist is 1
            else:
                time = math.sqrt(dist**3/(self.GG*mass))  # dist and mass are 1

        return Unit(dist, mass, time)

class Unit:
    def __init__(self, dist, mass, time):
        self.dist       = dist
        self.mass       = mass
        self.time       = time

        self.charge     = math.sqrt(mass*dist/PhysicalConstants.CGSMU0)
        self.Bfield     = mass/(time*self.charge)

        self.velocity   = dist/time
        self.accel      = dist/time/time
        self.density    = mass/dist**3
        self.pressure   = mass/(dist*time**2)
        self.ergg       = self.velocity**2
        self.energ      = mass*self.ergg
    def print(self):
        for key in self.__dict__.keys():
            print(key, ' = ', self.__dict__[key])
