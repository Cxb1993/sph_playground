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

    # xmin = 0.0
    # xmax = 0.0
    # ymin = 0.0
    # ymax = 0.0
    # zmin = 0.0
    # zmax = 0.0
    # resolution = -1
    # dimentions = -1
    # snapshotsNumber = -1
    # particlesPlacement  = "none"
    # particlesArrangment = "none"
    # initialConditions   = "none"
    # equationsSet        = "none"
    # secondDerivType     = "none"
    # resultFile          = "result.info"
    # finishTime = -1
    # hfac = -1
    kernel = "m4"
    stderrfile = "__stderr.info"
    stdoutfile = "__stdout.info"
    logfile = "__runslog.info"
    resultfile = ""
    influencefile = ""

    def toArgs(self):
        args = []
        args.append("--resolution")
        args.append(str(self.resolution))
        args.append("--dim")
        args.append(str(self.dimentions))
        args.append("--tfinish")
        args.append(str(self.finishTime))
        args.append("--hfac")
        args.append(str(self.hfac))
        args.append("--resultfile")
        args.append(str(self.resultFile))
        if (self.snapshotsNumber <= 0):
            args.append("--silent")
            args.append("yes")
        else:
            args.append("--nsnapshots")
            args.append(str(self.snapshotsNumber))
        if (self.particlesPlacement != "none"):
            args.append("--partplace")
            args.append(self.particlesPlacement)
            args.append("--box")
            args.append(str(self.xmin))
            args.append(str(self.xmax))
            args.append(str(self.ymin))
            args.append(str(self.ymax))
            args.append(str(self.zmin))
            args.append(str(self.zmax))
        if (self.particlesArrangment != "none"):
            args.append("--partarang")
            args.append(self.particlesArrangment)
        if (self.initialConditions != "none"):
            args.append("--initvar")
            args.append(self.initialConditions)
        if (self.equationsSet != "none"):
            args.append("--equations")
            args.append(self.equationsSet)
        if (self.secondDerivType != "none"):
            args.append("--ddw")
            args.append(self.secondDerivType)
        return args

class Context:
    pwd = os.getcwd()
    env = os.environ.copy()

    setup = Setup()
    rargs = RunnerArgument().args

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
        dtn = 0.5 if timeout > 0 else 0.0
        dtc = 0.5

        while ((tn <= t0 + timeout) and (process.poll() == None)):
            time.sleep(dtc)
            tn += dtn
            tc += dtc

        s = ""
        if (process.poll() == None):
            process.kill()
            s = "CMD " + " ".join(str(x) for x in process.args) + " KILLED AFTER " + str(tc - t0) + " SECONDS"
        else:
            s = "CMD " + " ".join(str(x) for x in process.args) + " EXECUTED IN " + str(tc - t0) + " SECONDS"
        so = process.stdout.read().decode('ascii')
        se = process.stderr.read().decode('ascii')
        self.__tolog__(s)
        if (so != ""):
            self.__tostdout__(s)
            self.__tostdout__(so)
        if (se != ""):
            self.__tostderr__(s)
            self.__tostderr__(se)

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
            "useomp=" + "t" if ("OMP_NUM_THREADS" in self.env) else "f"
            ], stdout=sp.PIPE, stderr=sp.PIPE, cwd=self.pwd, env=self.env)
        self.__waitwithtimeout__(make, 0.0)

    def ContinueRun(self):
        run = sp.Popen([
            "./execute"
        ], stdout=sp.PIPE, stderr=sp.PIPE, cwd=self.pwd, env=self.env)
        self.__waitwithtimeout__(run, 0.0)

    def SimpleRun(self):
        cmd = ["./execute"]
        cmd += self.setup.toArgs()
        run = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, cwd=self.pwd, env=self.env)
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

    CurrentTime = 0.0
    State = []
    Store = []

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

                curline = cf90.readline()

    def ReadFinalDump(self):
        with open(self.rargs.dumpmap, 'r') as fmap:
            fdump = FortranFile(self.rargs.dumpfile, 'r')
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
            self.CurrentTime = t
            self.State = state
            self.Store = store
            self.setup.resultfile = resultfile.split()[0]
            self.setup.influencefile = influencefile.split()[0]
            return True
        return False

    def IsParticlesPositionsCorrect(self):
        realPartNumb = int(self.State[self.ec['realpn']-1])
        for i in range(0,realPartNumb):
            rx = self.Store[self.partPropPos('rx',i)]
            ry = self.Store[self.partPropPos('ry',i)]
            rz = self.Store[self.partPropPos('rz',i)]
            if not ((rx >= self.setup.xmin) and (rx <= self.setup.xmax) and
                (ry >= self.setup.ymin) and (ry <= self.setup.ymax) and
                (rz >= self.setup.zmin) and (rz <= self.setup.zmax)):
                print("Particle " + str(i) + " is out of the box at [" + str(rx) + ";" +str(ry) + ";" + str(rz) +"]")
                return False
        return True

    def AfterRunUpdate(self):
        for i in range(1,len(self.State)+1):
            self.setup.__dict__[self.ecRev[i]] = self.State[i-1]

    def PrintState(self):
        print('\t +---field name ---+---run state---+---internal state---+')
        for i in range(1,len(self.State)+1):
            print("\t |", '{0:^15}'.format(self.ecRev[i]), "|", '{0:^13}'.format(self.State[i-1]),
                "|", '{0:^18}'.format(self.setup.__dict__[self.ecRev[i]]), "|")
        print('\t |   result file   |', '{0:^34}'.format(self.setup.resultfile), "|")
        print('\t +------------------------------------------------------+')

    def StateVal(self, key):
        return self.State[self.ec[key]-1]

    def partPropPos(self, propName, partIndx):
        return (partIndx*len(self.es) + self.es[propName]-1)

    def BackupOutput(self):
        os.system("cp -r " + self.pwd + "/output " + self.pwd + "/output.bck")
        os.system("mkdir -p " + self.pwd)

    def storeModify(state, store):
        for i in range(0,int(state[ec['realpn']-1])):
            # store[partPos('ax',i)] = 0.
            store[partPos('vy',i)] = 1e-4 * \
                math.sqrt(state[ec['gamma']-1] * store[partPos('p',i)] /
                store[partPos('den',i)]) * math.sin(4. * math.pi * store[partPos('rx',i)]/2.)
            # store[partPos('ay',i)] = 0.
            # store[partPos('az',i)] = 0.

            store[partPos('bx',i)] = 1e-11
            # store[partPos('bx',i)] = 0.
            store[partPos('by',i)] = 0.
            store[partPos('bz',i)] = 0.

        # print(store[partPos('type',1)])
        # print(store[partPos('rx',1)])
        # print(store[partPos('ry',1)])
        # print(store[partPos('rz',1)])

        state[ec['muzero']-1] = 1.
        state[ec['eqs']-1] = 303        #mhdd
        state[ec['ddw']-1] = 403        #2nw
        state[ec['disotropic']-1] = 1   #yes

        state[ec['tfinish']-1] = 1001.
        state[ec['dtprint']-1] = .0001
        return state, store

    def saveNewDump(t, state, store, resultfile, influencefile):
        with open(rargs.dumpmap[0], 'r') as fmap:
            fdumpnew = FortranFile(rargs.modifieddump[0], 'w')
            curmapline = fmap.readline()
            while (curmapline):
                splitted = curmapline.split()
                if len(splitted) != 4:
                    return False
                if splitted[0] == 't':
                    fdumpnew.write_record(numpy.array([t],dtype=numpy.float64))
                if splitted[0] == 'resultfile':
                    fdumpnew.write_record(
                        numpy.array(
                            (resultfile + influencefile).encode(encoding='utf_8')
                            )
                        )
                if splitted[0] == 'state':
                    fdumpnew.write_record(numpy.array(state,dtype=numpy.float64))
                if splitted[0] == 'store':
                    fdumpnew.write_record(numpy.array(store,dtype=numpy.float64))
                curmapline = fmap.readline()
            fdumpnew.close()
            return True
