import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from tap import Tap
from tqdm import tqdm
import logging


class SimpleArgParse(Tap):
    xyzdir: str
    """Directory with xyz files to compute"""
    basedir: str = os.path.dirname(os.path.realpath(__file__))
    """Defaults to file location directory"""
    xtb: bool = True
    stda: bool = True
    log: str = 'warning'


args = SimpleArgParse().parse_args()
levels = {
    'critical': logging.CRITICAL,
    'error': logging.ERROR,
    'warn': logging.WARNING,
    'warning': logging.WARNING,
    'info': logging.INFO,
    'debug': logging.DEBUG
}
logLevel = levels.get(args.log.lower())
basedir = args.basedir
xyzdir = os.path.join(basedir, args.xyzdir)
xyzdirname = xyzdir.split('/')[-1]
if len(xyzdirname) == 0:
    xyzdirname = xyzdir.split('/')[-2]
if ('_xyzfiles' in xyzdirname):
    xyzdirname = xyzdirname.replace('_xyzfiles', '')
elif ('xyzfiles_' in xyzdirname):
    xyzdirname = xyzdirname.replace('xyzfiles_', '')
resultsdir = os.path.join(basedir, 'xtb_results_' + xyzdirname)
try:
    os.mkdir(resultsdir)
except:
    pass
errFiles = []
logging.root.setLevel(logLevel)
logging.basicConfig(level=logLevel)


def moveFiles(currdir, basedir):
    logging.info('moving files')
    os.chdir(currdir)
    count = 0
    for dir in tqdm(sorted(os.listdir())):
        if (os.path.isfile(dir) and dir[0] != '.'):
            dirsplit = os.path.splitext(dir)
            dirname = dirsplit[0]
            dirext = dirsplit[1]
            if (dirext == '.xyz'):
                try:
                    os.mkdir(dirname)
                except:
                    pass
                os.system('mv -- \"' + dirname + '.xyz\" \"' +
                          dirname + '/' + dirname + '.xyz\"')
                count += 1
    logging.info('moved %s files', count)
    os.chdir(basedir)


def do_xtb(dir, currdir):
    os.chdir(currdir)
    recompute = True
    logging.info('\n')
    logging.info('starting xtb')
    if os.path.isdir(dir) and dir[0] != '.':
        os.chdir(dir)
        if os.path.isfile('output_xtb.out'):
            with open('output_xtb.out', 'r') as file:
                data = file.readlines()
                for line in data:
                    if 'finished run' in line:
                        logging.info("xtb already completed")
                        recompute = False
                        break

        chrg = '0'
        if recompute:
            logging.info('checking charge')
            os.system('xtb ./\"' + dir + ".xyz\" > output_xtb_setup.out 2>&1")
            with open('output_xtb_setup.out', 'r') as file:
                data = file.readlines()
            for line in data:
                if 'spin                       :                   0.5' in line:
                    chrg = '1'

        if recompute:
            logging.info('calculating triplet')
            os.system("xtb ./\"" + dir +
                      ".xyz\" --opt tight --gfn 2 --chrg " + chrg +
                      " --uhf 2 --gbsa acetonitrile > output_xtb_trip.out 2>&1")
            T1done = False
            with open('output_xtb_trip.out', 'r') as file:
                data = file.readlines()
                for line in data:
                    if 'finished run' in line:
                        T1done = True
            if not T1done:
                os.system("xtb ./\"" + dir + ".xyz\" --opt tight --gfn 1 --chrg "
                          + chrg + " --uhf 2 --gbsa acetonitrile > output_xtb_trip.out")
            os.system("grep 'TOTAL ENERGY' output_xtb_trip.out > T1ExXTB.dat")

        with open(os.path.join(currdir, dir, 'T1ExXTB.dat')) as data:
            datastr = data.read()
            temp = datastr.split()
            try:
                T1En_xtb = float(temp[3])
            except IndexError:
                if dir not in errFiles:
                    errFiles.append(dir)
                logging.info(dir)
                T1En_xtb = None
            except ValueError:
                if dir not in errFiles:
                    errFiles.append(dir)
                logging.info(dir)
                T1En_xtb = None

        if (recompute):
            logging.info('calculating ground')
            os.system("xtb ./\"" + dir +
                      ".xyz\" --opt tight --gfn 2 --chrg " + chrg + " --gbsa acetonitrile > output_xtb.out 2>&1")
            S0done = False
            with open('output_xtb.out', 'r') as file:
                data = file.readlines()
                for line in data:
                    if ('finished run' in line):
                        S0done = True
            if (not S0done):
                os.system("xtb ./\"" + dir +
                          ".xyz\" --opt tight --gfn 1 --chrg " + chrg + " --gbsa acetonitrile > output_xtb.out")
            os.system("grep 'TOTAL ENERGY' output_xtb.out > S0XTB.dat")

        with open(os.path.join(currdir, dir, 'S0XTB.dat')) as data:
            datastr = data.read()
            temp = datastr.split()
            try:
                S0En_xtb = float(temp[3])
            except IndexError:
                if dir not in errFiles:
                    errFiles.append(dir)
                S0En_xtb = None
            except ValueError:
                if dir not in errFiles:
                    errFiles.append(dir)
                S0En_xtb = None
        try:
            T1Ex_xtb = (T1En_xtb - S0En_xtb) * 27.2114
        except TypeError:
            if dir not in errFiles:
                errFiles.append(dir)
            T1Ex_xtb = None
    else:
        T1En_xtb = None
        S0En_xtb = None
        T1Ex_xtb = None
    return T1En_xtb, S0En_xtb, T1Ex_xtb


def do_stda(dir, currdir):
    os.chdir(currdir)
    os.chdir(dir)
    recompute = True
    logging.info('\n')
    logging.info('starting stda')
    if (os.path.isfile('wfn.xtb') and os.path.isfile('T1Ex_tddft.dat')):
        with open('T1Ex_tddft.dat', 'r') as file:
            data = file.readlines()
        if (len(data) >= 1):
            recompute = False
    if (not recompute):
        with open('S1Ex.dat') as data:
            datastr = data.read()
            temp = datastr.split()
            try:
                S1En_stda = float(temp[1])
            except:
                S1En_stda = None
        with open('S1Ex_tddft.dat', 'r') as data:
            datastr = data.read()
            temp = datastr.split()
            try:
                S1En_stddft = float(temp[1])
            except:
                S1En_stddft = None
        with open('T1Ex.dat') as data:
            datastr = data.read()
            temp = datastr.split()
            try:
                T1En_stda = float(temp[1])
            except:
                T1En_stda = None
        with open('T1Ex_tddft.dat') as data:
            datastr = data.read()
            temp = datastr.split()
            try:
                T1En_stddft = float(temp[1])
            except:
                T1En_stddft = None

    if recompute is False:
        logging.info('stda already completed')
        return S1En_stda, S1En_stddft, T1En_stda, T1En_stddft

    logging.info('setting up stda')
    os.system("xtb4stda xtbopt.xyz -gbsa acetonitrile > output_xtb4stda.out")
    logging.info('calculating s1 with stda')
    os.system("stda -xtb -e 10 > output_stda.out")  # stda
    os.system("grep -A1 state output_stda.out | tail -1 > S1Ex.dat")
    with open('S1Ex.dat') as data:
        datastr = data.read()
        temp = datastr.split()
        try:
            S1En_stda = float(temp[1])
        except:
            S1En_stda = None
    logging.info('calculating s1 with stddft')
    os.system("stda -xtb -e 10 -rpa > output_stddft.out")  # stddft
    os.system("grep -A1 state output_stddft.out | tail -1 > S1Ex_tddft.dat")
    with open('S1Ex_tddft.dat') as data:
        datastr = data.read()
        temp = datastr.split()
        try:
            S1En_stddft = float(temp[1])
        except:
            S1En_stddft = None
    logging.info('calculating t1 with stda')
    os.system("stda -xtb -t -e 10 > output_stdaxtb_trip.out")
    os.system("grep -A1 state output_stdaxtb_trip.out | tail -1 > T1Ex.dat")
    with open('T1Ex.dat') as data:
        datastr = data.read()
        temp = datastr.split()
        try:
            T1En_stda = float(temp[1])
        except:
            T1En_stda = None
    logging.info('calculating t1 with stddft')
    os.system("stda -xtb -t -e 10 -rpa > output_stddft_trip.out")
    os.system("grep -A1 state output_stddft_trip.out | tail -1 > T1Ex_tddft.dat")
    with open('T1Ex_tddft.dat') as data:
        datastr = data.read()
        temp = datastr.split()
        try:
            T1En_stddft = float(temp[1])
        except:
            T1En_stddft = None
    return S1En_stda, S1En_stddft, T1En_stda, T1En_stddft


def run_all(xyzdir, basedir):
    count = 0
    os.chdir(resultsdir)
    # try:
    #     with open('exData_' + xyzdirname + '.txt', 'r') as file:
    #         exData = json.load(file)
    # except:
    #     exData = {}
    if not os.path.isfile('exData_' + xyzdirname + '.csv'):
        with open('exData_' + xyzdirname + '.csv', 'w') as file:
            file.write('CID,T1xtb,T1stda,S1stda\n')
    try:
        with open('startInd_xtb_' + xyzdirname + '.txt', 'r') as file:
            data = file.readlines()
        startInd = int(data[0].replace('\n', ''))
    except:
        startInd = 0
    with open('PCQC_175k_dirs_sorted.txt', 'r') as file:
        listdirs = json.load(file)
    os.chdir(xyzdir)
    # listdirs = []
    # for dir in sorted(os.listdir()):
    #     if (os.path.isdir(dir) and dir[0] != '.'):
    #         listdirs.append(dir)
    for index, dir in enumerate(tqdm(listdirs)):
        if (index < startInd):
            continue
        isError = False

        # xtb
        # if (args.xtb):
        T1En_xtb, S0En_xtb, T1Ex_xtb = do_xtb(dir, xyzdir)
        if (T1Ex_xtb is not None):
            if (T1Ex_xtb < 1e-3 or T1Ex_xtb > 10):
                isError = True
            # exData[dir] = {}
            # exData[dir]["xtb"] = {}
            # exData[dir]["xtb"]["S0_abs"] = S0En_xtb
            # exData[dir]["xtb"]["T1_abs"] = T1En_xtb
            # exData[dir]["xtb"]["T1"] = T1Ex_xtb
        else:
            isError = True

        # stda
        # if (args.stda):
        S1En_stda, S1En_stddft, T1En_stda, T1En_stddft = do_stda(
            dir, xyzdir)

        # try:
        #     T1Ex_xtb = exData[dir]["xtb"]["T1"]
        # except KeyError:
        #     exData[dir] = {}
        # exData[dir]["std"] = {}
        # exData[dir]["std"]["T1_stda"] = T1En_stda
        # exData[dir]["std"]["T1_stddft"] = T1En_stddft
        # exData[dir]["std"]["S1_stda"] = S1En_stda
        # exData[dir]["std"]["S1_stddft"] = S1En_stddft
        if (S1En_stda is None or S1En_stddft is None or T1En_stda is None or T1En_stddft is None):
            isError = True
        elif (S1En_stda < 1e-3 or S1En_stda > 10):
            isError = True
        if (isError):
            errFiles.append(dir)
        logging.info('\n')
        logging.info('calculations complete for molecule %s exporting data',
                     index + 1)
        os.chdir(resultsdir)
        # if index % 10 == 0:
        # with open('exData_' + xyzdirname + '.txt', 'w') as file:
        #     json.dump(exData, file)
        if not isError:
            with open('exData_' + xyzdirname + '.csv', 'a') as file:
                file.write(dir + ',' + str(T1Ex_xtb) + ',' + str(T1En_stda) + ',' + str(S1En_stda) + '\n')
        with open('exData_errs_' + xyzdirname + '.txt', 'w') as file:
            json.dump(errFiles, file)
        with open('startInd_xtb_' + xyzdirname + '.txt', 'w') as file:
            file.write(str(index))
        os.chdir(xyzdir)
    os.chdir(basedir)
    # return exData


def main():
    # moveFiles(xyzdir, basedir)
    run_all(xyzdir, basedir)
    logging.info('succesfully completed')
    if (logging.root.level > 20):
        print('succesfully completed')


if __name__ == '__main__':
    main()
