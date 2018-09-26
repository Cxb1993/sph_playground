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
    resultfile = ""
    influencefile = ""
    time = 0.0

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
        self.env["OMP_NUM_THREADS"] = str(i)

    def GetDateTimeString(self):
        strres = str(datetime.datetime.today().year)
        strres += "-" + str(datetime.datetime.today().month)
        strres += "-" + str(datetime.datetime.today().day)
        strres += "-" + str(datetime.datetime.today().hour)
        strres += "-" + str(datetime.datetime.today().minute)
        strres += "-" + str(datetime.datetime.today().second)
        return strres

    def SimpleMake(self):
        make = sp.Popen([
            "make",
            "debug=f",
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

    def MakeAndRun(self):
        self.SimpleMake()
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
########                                              #########
########                                              #########
########                                              #########
########                                              #########
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

    def ReadFinalDump(self):
        self.ReadDumpFile("./output/dumpmap", "./output/finaldump")

    def ReadDumpFile(self, dumpmap, dumpfile):
        with open(dumpmap, 'r') as fmap:
            fdump = FortranFile(dumpfile, 'r')
            curmapline = fmap.readline()
            while (curmapline):
                splitted = curmapline.split()
                if len(splitted) != 4:
                    return False
                if splitted[0] == 't':
                    t = fdump.read_reals(dtype=numpy.float64)[0]
                if splitted[0] == 'resultfile':
                    resultfile = "".join(map(chr, fdump.read_record(dtype=numpy.uint8)))
                    influencefile = resultfile[50:]
                    resultfile = resultfile[:50]
                if splitted[0] == 'state':
                    state = fdump.read_reals(dtype=numpy.float64)
                if splitted[0] == 'store':
                    store = fdump.read_reals(dtype=numpy.float64)
                curmapline = fmap.readline()

            for i in range(1,len(state)+1):
                self.setup.__dict__[self.ecRev[i]] = state[i-1]
            self.store = store
            self.setup.time = t
            self.setup.resultfile = resultfile
            self.setup.influencefile = influencefile
            return True
        return False

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

    def PrintState(self):
        print("\t +", '{0:-^15}'.format("field name"),"+", '{0:-^40}'.format("state"),"+")
        for i in range(1,len(self.setup.__dict__.keys())):
            if i in self.ecRev:
                print("\t |", '{0:^15}'.format(self.ecRev[i]), "|", '{0:^40}'.format(self.setup.__dict__[self.ecRev[i]]), "|")
        print('\t |  current time   |', '{0:^40}'.format(self.setup.time), "|")
        print('\t |  result file    |','{0:^40}'.format(self.setup.resultfile.split()[0]), "|")
        print('\t |  influence file |','{0:^40}'.format(self.setup.influencefile.split()[0]), "|")
        print("\t +", '{0:-^15}'.format(""), "+", '{0:-^40}'.format(""),"+")

    def PrintPropertyKeys(self):
        print("Properties: " + str(self.es.keys()))

    def StateVal(self, key):
        return self.State[self.ec[key]-1]

    def idx(self, propName, partIndx):
        # index in threaded array from squared data
        return (partIndx*len(self.es) + self.es[propName]-1)

    def BackupOutput(self, str):
        if (len(str.split()) != 1):
            print("Cannot copy with " + str)
        else:
            os.system("cp -r " + self.pwd + "/output " + self.pwd + "/" + str)

    def BackupDumps(self, strDump, strMap):
        if ((len(strDump.split()) != 1) or (len(strMap.split()) != 1)):
            print("Cannot copy with " + strDump + " | " + strMap)
        else:
            os.system("cp " + self.pwd + "/output/finaldump " + self.pwd + "/" + strDump)
            os.system("cp " + self.pwd + "/output/dumpmap " + self.pwd + "/" + strMap)

    def AddParticles(self, dim, type, method):
        locposidx  = self.es['r'+dim]-1
        loctypeidx = self.es['type']-1
        bd   = self.setup.bordsize
        min  = self.setup.__dict__[dim + 'min']
        max  = self.setup.__dict__[dim + 'max']
        bmin = min + bd
        bmax = max - bd

        storenew = copy.deepcopy(self.store)
        nreal = int(self.setup.realpn)
        idx = self.idx
        totprop = len(self.es)
        addedNumber = 0
        for i in range(0,nreal):
            oldr = self.store[idx('r'+dim,i)]
            newr = oldr
            if (oldr < bmin):
                newr = max + oldr - min
            if (oldr > bmax):
                newr = min + oldr - max
            if (abs(oldr-newr) > 10e-5):
                p1 = idx('rx',i)
                p2 = p1 + totprop
                newp = copy.deepcopy(self.store[p1:p2])
                newp[locposidx]  = newr
                newp[loctypeidx] = self.ept[type]
                storenew = numpy.append(storenew, newp)
                addedNumber += 1
        self.store = copy.deepcopy(storenew)
        return float(addedNumber)

    def CopyContext(self):
        newContext = copy.deepcopy(self)
        newContext.setup = copy.deepcopy(self.setup)
        newContext.store = copy.deepcopy(self.store)
        return newContext

    def ModifyParticles(self, condition, condarg, properties, value, valuearg):
        storenew = copy.deepcopy(self.store)
        nreal = int(self.setup.realpn)
        nfixd = int(self.setup.fixedpn)
        for i in range(0,nreal+nfixd):
            if condition(self.store[self.idx(condarg,i)]):
                for pa in properties:
                    storenew[self.idx(pa,i)] = value(self.store[self.idx(valuearg,i)])
        self.store = copy.deepcopy(storenew)

    def Apply(self):
        with open(self.rargs.dumpmap, 'r') as fmap:
            fdumpnew = FortranFile(self.rargs.dumpfile, 'w')
            curmapline = fmap.readline()
            while (curmapline):
                splitted = curmapline.split()
                if len(splitted) != 4:
                    return False
                if splitted[0] == 't':
                    fdumpnew.write_record(numpy.array([self.setup.time],dtype=numpy.float64))
                if splitted[0] == 'resultfile':
                    fdumpnew.write_record(
                        numpy.array((
                        '{0:<50}'.format(self.setup.resultfile) +
                        '{0:<50}'.format(self.setup.influencefile)
                        ).encode(encoding='utf_8')))
                if splitted[0] == 'state':
                    state = []
                    for i in self.ecRev:
                        state.append(self.setup.__dict__[self.ecRev[i]])
                    fdumpnew.write_record(numpy.array(state,dtype=numpy.float64))
                if splitted[0] == 'store':
                    fdumpnew.write_record(numpy.array(self.store,dtype=numpy.float64))
                curmapline = fmap.readline()
