#!/usr/bin/env python

import argparse
from scipy.io import FortranFile
import numpy
import time
import datetime
import math
import os
import subprocess as sp
import time
import copy
import h5py

class RunnerArgument:
    parser = argparse.ArgumentParser(description='I\'m the law.')

    parser.add_argument('-c', '--constf90', type=str, nargs='?',
                        help='fortran file with constant descriptions',
                        default='./src/const.f90')
    parser.add_argument('-m', '--dumpmap', type=str, nargs='?',
                        help='dumpmap file',
                        default='./output/dumpmap')
    parser.add_argument('-o', '--dumpfile', type=str, nargs='?',
                        help='dump file',
                        default='./output/finaldump')
    parser.add_argument('-n', '--modifieddump', type=str, nargs='?',
                        help='new dump file',
                        default='./output/finaldump')

    args = parser.parse_args()

class Setup:
    stderrfile = "__stderr.info"
    stdoutfile = "__stdout.info"
    logfile = "__runslog.info"
    kernel = "m4"

    def toArgs(self):
        args = []
        for k in self.__dict__.keys():
            args.append("--"+k)
            args.append(str(self.__dict__[k]))
        return args

class Context:
    pwd = os.getcwd()
    env = os.environ.copy()
    rargs = RunnerArgument().args

    setup = Setup()
    store = []
    state = []

    def __init__(self):
        self.ReadConstFile()
        for key in self.ec.keys():
            self.setup.__dict__[key] = 0.0

    def __tolog__(self, s):
        if self.setup.logfile != "none":
            f = open(self.pwd + "/" +self.setup.logfile, 'a+')
            f.write(s + "\n")
            f.close()
    def __tostdout__(self, s):
        if self.setup.stdoutfile != "none":
            f = open(self.pwd + "/" +self.setup.stdoutfile, 'a+')
            f.write(s + "\n")
            f.close()
    def __tostderr__(self, s):
        if self.setup.stderrfile != "none":
            f = open(self.pwd + "/" + self.setup.stderrfile, 'a+')
            f.write(s + "\n")
            f.close()

    def __waitwithtimeout__(self, process, timeout):
        t0 = time.time()
        tn = t0
        tc = t0
        dtn = 1.0 if timeout > 0 else 0.0
        dtc = 1.0
        s = "CMD [" + " ".join(str(x) for x in process.args) + "] SENT TO RUN"
        self.__tolog__(s)

        while ((tn <= t0 + timeout) and (process.poll() == None)):
            time.sleep(dtc)
            tn += dtn
            tc += dtc

        s = ""
        if (process.poll() == None):
            process.kill()
            s = "CMD [" + " ".join(str(x) for x in process.args) + "] KILLED AFTER " + str(tc - t0) + " SECONDS"
        else:
            s = "CMD [" + " ".join(str(x) for x in process.args) + "] EXECUTED IN " + str(tc - t0) + " SECONDS"
        so = process.stdout.read().decode('ascii')
        # se = process.stderr.read().decode('ascii')
        self.__tolog__(s)
        if (so != ""):
            self.__tostdout__(s)
            self.__tostdout__(so)
        # if (se != ""):
        #     self.__tostderr__(s)
        #     self.__tostderr__(se)
        return so

    def CleanSTDLogs(self):
        if self.setup.stdoutfile != "none":
            f = open(self.pwd + "/" +self.setup.stdoutfile, 'w+')
            f.write("")
            f.close()
        if self.setup.stderrfile != "none":
            f = open(self.pwd + "/" + self.setup.stderrfile, 'w+')
            f.write("")
            f.close()
    def CleanAllLogs(self):
        if self.setup.logfile != "none":
            f = open(self.pwd + "/" +self.setup.logfile, 'w+')
            f.write("")
            f.close()
        if self.setup.stdoutfile != "none":
            f = open(self.pwd + "/" +self.setup.stdoutfile, 'w+')
            f.write("")
            f.close()
        if self.setup.stderrfile != "none":
            f = open(self.pwd + "/" + self.setup.stderrfile, 'w+')
            f.write("")
            f.close()

    def SetThreadsOMP(self, i):
        if (i != 1):
            self.env["OMP_NUM_THREADS"] = str(i)

    def GetDateTimeString(self):
        strres = str(datetime.datetime.today().year)
        strres += "-" + str(datetime.datetime.today().month)
        strres += "-" + str(datetime.datetime.today().day)
        strres += "-" + str(datetime.datetime.today().hour)
        strres += "-" + str(datetime.datetime.today().minute)
        strres += "-" + str(datetime.datetime.today().second)
        return strres

    def SimpleMake(self, debug=False):
        make = sp.Popen([
            "make",
            "debug=t" if debug else "debug=f",
            "kernel=" + self.setup.kernel,
            "useomp=t" if ("OMP_NUM_THREADS" in self.env) else "useomp=f"
            ], stdout=sp.PIPE, cwd=self.pwd, env=self.env)
        self.__waitwithtimeout__(make, 0.0)

    def ContinueRun(self):
        run = sp.Popen([
            "./execute"
        ], stdout=sp.PIPE, stderr=sp.PIPE, cwd=self.pwd, env=self.env)
        self.__waitwithtimeout__(run, 0.0)

    def SimpleRun(self):
        cmd = ["./execute"]
        cmd += self.setup.toArgs()
        run = sp.Popen(cmd, stdout=sp.PIPE, cwd=self.pwd, env=self.env)
        self.__waitwithtimeout__(run, 0.0)

    def CleanRun(self):
        if not os.path.exists(self.pwd + "/output"):
            os.mkdir(self.pwd + "/output")
        else:
            os.system("rm -rf " + self.pwd + "/output/*")
        self.SimpleRun()

    def CleanFinalDumps(self):
        if not os.path.exists(self.pwd + "/output"):
            os.mkdir(self.pwd + "/output")
        else:
            os.system("rm -rf " + self.pwd + "/output/finaldump")
            os.system("rm -rf " + self.pwd + "/output/fulldump")

    def MakeAndRun(self, debug=False):
        self.SimpleMake(debug)
        self.CleanRun()

    def PrintCell(self, file, row, column):
        with open(file.split()[0], 'r') as f:
            curline = f.readline()
            lastline = curline
            r = 1
            while(curline and (r != row)):
                lastline = curline
                curline = f.readline()
                r += 1
            splitted = lastline.split()
            return splitted[column-1]

###############################################################
###############################################################

    # All the enums from the const.f90
    # fields in state array
    ec      = {}
    ecRev   = {}
    # fields in store arrat
    es      = {}
    esRev   = {}
    # particle type
    ept     = {}
    eptRev  = {}
    # initial setups
    ett     = {}
    ettRev  = {}
    # border enumerators
    ebc     = {}
    ebcRev  = {}
    # equation types
    eeq     = {}
    eeqRev  = {}
    # second derivatives
    esd     = {}
    esdRev  = {}
    # an/inosotropy type
    edi     = {}
    ediRev  = {}
    # components of equation array
    eqs     = {}
    eqsRev  = {}
    # processess described in iterate
    epc     = {}
    epcRev  = {}
    # yes / no constants
    eif     = {}
    eifRev  = {}

    def ReadConstFile(self):
        with open(self.rargs.constf90, 'r') as cf90:
            curline = cf90.readline()
            while(curline):
                splitted = curline.split()
                if (len(splitted)==3):
                    if((splitted[0][0:3] == "ec_") and (splitted[2][-2:] == ",&")):
                        self.ec[splitted[0][3:]] = int(splitted[2][:-2])
                        self.ecRev[int(splitted[2][:-2])] = splitted[0][3:]
                    if((splitted[0][0:3] == "es_") and (splitted[2][-2:] == ",&")):
                        if (splitted[0][3:] != "total"):
                            self.es[splitted[0][3:]] = int(splitted[2][:-2])
                            self.esRev[int(splitted[2][:-2])] = splitted[0][3:]
                    if (splitted[0][0:4] == "ept_"):
                        if (splitted[2][-2:] == ",&"):
                            self.ept[splitted[0][4:]] = int(splitted[2][:-2])
                            self.eptRev[int(splitted[2][:-2])] = splitted[0][4:]
                        else:
                            self.ept[splitted[0][4:]] = int(splitted[2])
                            self.eptRev[int(splitted[2])] = splitted[0][4:]
                    if (splitted[0][0:4] == "ett_"):
                        if (splitted[2][-2:] == ",&"):
                            self.ett[splitted[0][4:]] = int(splitted[2][:-2])
                            self.ettRev[int(splitted[2][:-2])] = splitted[0][4:]
                        else:
                            self.ett[splitted[0][4:]] = int(splitted[2])
                            self.ettRev[int(splitted[2])] = splitted[0][4:]
                    if (splitted[0][0:4] == "ebc_"):
                        if (splitted[2][-2:] == ",&"):
                            self.ebc[splitted[0][4:]] = int(splitted[2][:-2])
                            self.ebcRev[int(splitted[2][:-2])] = splitted[0][4:]
                        else:
                            self.ebc[splitted[0][4:]] = int(splitted[2])
                            self.ebcRev[int(splitted[2])] = splitted[0][4:]
                    if (splitted[0][0:4] == "eeq_"):
                        if (splitted[2][-2:] == ",&"):
                            self.eeq[splitted[0][4:]] = int(splitted[2][:-2])
                            self.eeqRev[int(splitted[2][:-2])] = splitted[0][4:]
                        else:
                            self.eeq[splitted[0][4:]] = int(splitted[2])
                            self.eeqRev[int(splitted[2])] = splitted[0][4:]
                    if (splitted[0][0:4] == "esd_"):
                        if (splitted[2][-2:] == ",&"):
                            self.esd[splitted[0][4:]] = int(splitted[2][:-2])
                            self.esdRev[int(splitted[2][:-2])] = splitted[0][4:]
                        else:
                            self.esd[splitted[0][4:]] = int(splitted[2])
                            self.esdRev[int(splitted[2])] = splitted[0][4:]
                    if (splitted[0][0:4] == "edi_"):
                        if (splitted[2][-2:] == ",&"):
                            self.edi[splitted[0][4:]] = int(splitted[2][:-2])
                            self.ediRev[int(splitted[2][:-2])] = splitted[0][4:]
                        else:
                            self.edi[splitted[0][4:]] = int(splitted[2])
                            self.ediRev[int(splitted[2])] = splitted[0][4:]
                    if (splitted[0][0:4] == "eqs_"):
                        if (splitted[2][-2:] == ",&"):
                            self.eqs[splitted[0][4:]] = int(splitted[2][:-2])
                            self.eqsRev[int(splitted[2][:-2])] = splitted[0][4:]
                        else:
                            self.eqs[splitted[0][4:]] = int(splitted[2])
                            self.eqsRev[int(splitted[2])] = splitted[0][4:]
                    if (splitted[0][0:4] == "epc_"):
                        if (splitted[2][-2:] == ",&"):
                            self.epc[splitted[0][4:]] = int(splitted[2][:-2])
                            self.epcRev[int(splitted[2][:-2])] = splitted[0][4:]
                        else:
                            self.epc[splitted[0][4:]] = int(splitted[2])
                            self.epcRev[int(splitted[2])] = splitted[0][4:]
                    if (splitted[0][0:4] == "eif_"):
                        if (splitted[2][-2:] == ",&"):
                            self.eif[splitted[0][4:]] = float(splitted[2][:-2])
                            self.eifRev[float(splitted[2][:-2])] = splitted[0][4:]
                        else:
                            self.eif[splitted[0][4:]] = float(splitted[2])
                            self.eifRev[float(splitted[2])] = splitted[0][4:]

                curline = cf90.readline()

    def ReadDumpHDF(self, filename):
        self.file = h5py.File(filename, 'r+')
        self.store = self.file['store']
        self.state = self.file['state']
        for i in range(1,len(list(self.state))+1):
            self.setup.__dict__[self.ecRev[i]] = self.state[i-1]

    def PrintSetup(self):
        print("\t +", '{0:-^15}'.format("field name"),"+", '{0:-^40}'.format("setup"),"+")
        for i in self.ecRev:
            e = self.ecRev[i]
            if e in self.setup.__dict__.keys():
                print("\t |", '{0:^15}'.format(e), "|",
                '{0:^40}'.format(self.setup.__dict__[e]), "|")
        print("\t +", '{0:-^15}'.format(""), "+", '{0:-^40}'.format(""),"+")

    def PrintState(self):
        print("\t +", '{0:-^15}'.format("field name"),"+", '{0:-^40}'.format("state"),"+")
        for i in self.ecRev:
            print("\t |", '{0:^15}'.format(self.ecRev[i]), "|",
                '{0:^40}'.format(self.state[i-1]), "|")
        print("\t +", '{0:-^15}'.format(""), "+", '{0:-^40}'.format(""),"+")

    def ParticlesWithWrongPositions(self):
        realPartNumb = int(self.setup.realpn)
        idx = []
        for i in range(0,realPartNumb):
            irx = self.store[self.idx('rx',i)]
            iry = self.store[self.idx('ry',i)]
            irz = self.store[self.idx('rz',i)]
            if not ((irx >= self.setup.xmin) and (irx <= self.setup.xmax) and
                (iry >= self.setup.ymin) and (iry <= self.setup.ymax) and
                (irz >= self.setup.zmin) and (irz <= self.setup.zmax)):
                print("Particle " + str(i) + "r = [" +
                    str(irx) + ";" +str(iry) + ";" + str(irz) +"] is out of the box [" +
                    str(self.setup.xmin) + ":" + str(self.setup.xmax) +"X"+
                    str(self.setup.ymin) + ":" + str(self.setup.ymax) +"X"+
                    str(self.setup.zmin) + ":" + str(self.setup.zmax) +"]")
                idx.append(i)
        return idx

    def PrintPropertyKeys(self):
        print("Properties: " + str(self.es.keys()))

    def GetState(self, key):
        return self.state[self.ec[key]-1]
    def SetState(self, key, val):
        self.state[self.ec[key]-1] = val

    def BackupOutput(self, str):
        if (len(str.split()) != 1):
            print("Cannot copy with " + str)
        else:
            # print("cp -r " + self.pwd + "/output " + self.pwd + "/" + str)
            os.system("cp -r " + self.pwd + "/output " + self.pwd + "/" + str)

    def BackupDump(self, file):
        if (len(file.split()) != 1):
            print("Cannot copy with [" + file + "]")
        else:
            os.system("cp " + self.pwd + "/" + file + " " + self.pwd + "/" + file + ".bck")

    def AddParticles(self, dim:str, type:str, method:str):
        prop  = 'r'+dim
        iprop = self.es[prop]-1
        ftype  = self.ept[type]
        itype = self.es['type']-1

        bs  = self.setup.bordsize
        min = self.setup.__dict__[dim + 'min']
        max = self.setup.__dict__[dim + 'max']
        res = self.setup.resolution
        dr  = (max-min)/res
        bmin = min + bs
        bmax = max - bs

        newstore = self.store[...]
        nreal = int(self.setup.realpn)
        addedNumber = 0
        for i in range(0,nreal):
            oldr = newstore[i][iprop]
            newr = oldr
            if (oldr < bmin):
                newr = max + oldr - min
            if (oldr > bmax):
                newr = min + oldr - max
            if (abs(oldr-newr) > dr*1e-3):
                p = newstore[i].copy()
                p[itype] = ftype
                p[iprop] = newr
                newstore = numpy.append(newstore, [p], axis=0)
                addedNumber += 1
        self.store.resize(self.store.shape[0]+addedNumber, axis=0)
        self.store[...] = newstore
        if (type == 'fixed') or (type == 'fixedreal'):
            self.setup.fixedpn = addedNumber

    def CalcUniformMass(self, rho):
        mass = 1.
        lx = ly = lz = 1.0
        if (self.setup.dim > 0):
            lx = abs(self.setup.xmax - self.setup.xmin)
        if (self.setup.dim > 1):
            ly = abs(self.setup.ymax - self.setup.ymin)
        if (self.setup.dim > 2):
            lz = abs(self.setup.zmax - self.setup.zmin)
        mass = lx*ly*lz*rho/self.setup.realpn
        return mass

    def CopyContext(self):
        newContext = copy.deepcopy(self)
        newContext.setup = copy.deepcopy(self.setup)
        newContext.store = copy.deepcopy(self.store)
        return newContext

    def ModifyParticlesProperties(self,
        condition = [lambda rx: True],
        conditionarg = ['rx'],
        properties = ['rx'],
        value = [lambda rx: 1.0],
        valuearg = ['rx']):

        storenew = self.store[...]

        nreal = int(self.setup.realpn)
        nfixd = int(self.setup.fixedpn)
        nc = len(condition)
        ncx = len(conditionarg)
        np = len(properties)
        nv = len(value)
        nvx = len(valuearg)

        argval = lambda arg, i: self.store[i, self.es[arg]-1]
        iprop = lambda propname : self.es[propname]-1

        for i in range(0,nreal+nfixd):
            for j in range(np):
                fargj = valuearg[0] if nvx == 1 else valuearg[j]
                x1 = argval(fargj,i)
                storenew[i, iprop(properties[j])] = value[j](x1)

        self.store[...] = storenew
        self.file.flush()

    def Apply(self):
        for i in self.ecRev:
            e = self.ecRev[i]
            if e in self.setup.__dict__.keys():
                self.state[i-1] = self.setup.__dict__[e]
        self.file.flush()
