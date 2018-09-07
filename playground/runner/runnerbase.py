#!/usr/bin/env python

import argparse
from scipy.io import FortranFile
import numpy
import time
import math

class Context:

    parser = argparse.ArgumentParser(description='Reformator for eva.')

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

    cState={}
    revcState={}
    cStore={}
    revcStore={}

    def ReadConsts(self):
        with open(self.args.constf90, 'r') as cf90:
            curline = cf90.readline()
            while(curline):
                splitted = curline.split()
                if (len(splitted)==3):
                    if((splitted[0][0:3] == "ec_") and (splitted[2][-2:] == ",&")):
                        self.cState[splitted[0][3:]] = int(splitted[2][:-2])
                        self.revcState[int(splitted[2][:-2])] = splitted[0][3:]
                    if((splitted[0][0:3] == "es_") and (splitted[2][-2:] == ",&")):
                        if (splitted[0][3:] != "total"):
                            self.cStore[splitted[0][3:]] = int(splitted[2][:-2])
                            self.revcStore[int(splitted[2][:-2])] = splitted[0][3:]
                curline = cf90.readline()

    def Read(self):
        self.ReadConsts()
        with open(self.args.dumpmap, 'r') as fmap:
            fdump = FortranFile(self.args.dumpfile, 'r')
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
            return True, t, state, store, resultfile, influencefile
        return False

    def saveNewDump(t, state, store, resultfile, influencefile):
        with open(args.dumpmap[0], 'r') as fmap:
            fdumpnew = FortranFile(args.modifieddump[0], 'w')
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

    def printState(self, state):
        print('# # State info')
        for i in range(1,len(state)+1):
            print('\t + ' + self.revcState[i] + ": " + str(state[i-1]))
        print('# # Store ids')
        for i in range(1,len(self.revcStore)+1):
            print('\t + ' + self.revcStore[i] + ": " + str(i))

    def storeModify(state, store):
        for i in range(0,int(state[cState['realpn']-1])):
            # store[partPos('ax',i)] = 0.
            store[partPos('vy',i)] = 1e-4 * \
                math.sqrt(state[cState['gamma']-1] * store[partPos('p',i)] /
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

        state[cState['muzero']-1] = 1.
        state[cState['eqs']-1] = 303        #mhdd
        state[cState['ddw']-1] = 403        #2nw
        state[cState['disotropic']-1] = 1   #yes

        state[cState['tfinish']-1] = 1001.
        state[cState['dtprint']-1] = .0001
        return state, store

    def partPos(propName, partIndx):
        return (partIndx*len(cStore) + cStore[propName]-1)
