with open('progress.out', 'w') as file:
    file.write('importing\n')
print('importing')
import ray
from ray.exceptions import GetTimeoutError
import os
import sys
import logging
import time
import json
# from filelock import FileLock
# from p_tqdm import p_map
from tap import Tap
from openbabel import pybel
from openbabel import openbabel as ob


class SimpleArgParse(Tap):
    # xyzdir: str
    # """Directory with xyz files to compute"""
    smiles_path: str
    # basedir: str = os.path.dirname(os.path.realpath(__file__))
    basedir: str = os.path.dirname('/rds/general/user/sv920/home/Databases/SCOP_DB/')
    ephemdir: str = os.path.dirname('/rds/general/user/sv920/ephemeral/')
    """Defaults to SCOP_DB directory"""
    log: str = 'warning'


def shell(cmd, shell=False):
    if shell:
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = cmd.split()
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate()
    return output


@ray.remote(num_cpus=32)
def gen_sdf(smiles, dir):
    if os.path.isfile(os.path.join(xyzdir, dir, dir + '.xyz')):
        return
    if os.path.isfile(os.path.join(xyzdir, dir + '.sdf')):
        pass
        # return
    mol2 = pybel.readstring('smi', smiles)
    gen3d = ob.OBOp.FindType("gen3D")
    gen3d.Do(mol2.OBMol, "--medium")
    mol2.write('sdf', filename=os.path.join(xyzdir, dir + '.sdf'), overwrite=True)
    if not os.path.isfile(os.path.join(xyzdir, dir + '.sdf')):
        mol2.write('sdf', filename=os.path.join(xyzdir, dir + '.sdf'), overwrite=True)
    with open(os.path.join(xyzdir, dir + '.sdf'), 'r') as file:
        data = file.readlines()
    doGen2D = False
    for line in data:
        if '0.0000    0.0000    0.0000' in line:
            doGen2D = True
            break
    if (doGen2D):
        mol2 = pybel.readstring('smi', smiles)
        mol2.addh()
        gen2d = ob.OBOp.FindType("gen2D")
        gen2d.Do(mol2.OBMol, "--medium")
        mol2.write('sdf', filename=os.path.join(xyzdir, dir + '.sdf'), overwrite=True)
    with open(os.path.join(xyzdir, dir + '.sdf'), 'r') as file:
        data = file.readlines()
    doLocalOpt = False
    for line in data:
        if '0.0000    0.0000    0.0000' in line:
            doLocalOpt = True
            break
    if (doLocalOpt):
        mol2 = pybel.readstring('smi', smiles)
        mol2.addh()
        mol2.make3D()
        mol2.write('sdf', filename=os.path.join(xyzdir, dir + '.sdf'), overwrite=True)


@ray.remote(num_cpus=32)
def gen_xyz(filename, isFile=True):
    if os.path.isfile(os.path.join(xyzdir, filename, filename + '.xyz')):
        return
    file = filename + ".sdf"
    recompute = True
    try:
        os.mkdir(os.path.join(xyzdir, filename))
    except:
        if os.path.isfile(os.path.join(xyzdir, filename, filename + '.xyz')):
            # recompute = True
            recompute = False
    if recompute:
        sdffile = os.path.join(xyzdir, filename, file)
        if isFile:
            try:
                os.mkdir(os.path.join(xyzdir, filename))
            except:
                pass
            os.system('mv ' + os.path.join(xyzdir, file) + ' ' + sdffile)
        mol = next(pybel.readfile('sdf', sdffile))
        chrg = str(mol.charge)
        os.system('(cd ' + os.path.join(xyzdir, filename) +
                  ' && xtb ' + file + ' --opt tight --gfn 2 --charge ' + chrg + ' --gbsa benzene --norestart > '
                  + 'output_xtb.out)')
        if not os.path.isfile(os.path.join(xyzdir, filename, 'xtbopt.sdf')):
            os.system('(cd ' + os.path.join(xyzdir, filename) +
                      ' && xtb ' + file + ' --opt tight --gfn 1 --charge ' + chrg + ' --gbsa benzene --norestart > '
                      + 'output_xtb.out)')
        os.system('(cd ' + os.path.join(xyzdir, filename) +
                  ' && obabel -isdf ' + 'xtbopt.sdf' + ' -O xtbopt.xyz)')

        os.system('(cd ' + os.path.join(xyzdir, filename) +
                  ' && cp xtbopt.xyz ' + filename + '.xyz)')
        os.system('(cd ' + os.path.join(xyzdir, filename) +
                  " && find . -type f -not -name 'ID*.xyz' -exec rm -rf {} \;)")
    else:
        if isFile:
            os.system('rm ' + os.path.join(xyzdir, file))
        os.system('(cd ' + os.path.join(xyzdir, filename) +
                  " && find . -type f -not -name 'ID*.xyz' -exec rm -rf {} \;)")


@ray.remote(num_cpus=32)
def run_xtb_stda(filename):
    xyzfile = os.path.join(xyzdir, filename, filename + '.xyz')
    if not os.path.isfile(xyzfile):
        return
    if os.path.isfile(os.path.join(xyzdir, filename, 'S1Ex.dat')) and \
            os.path.isfile(os.path.join(xyzdir, filename, 'T1Ex.dat')):
        return
    try:
        mol = next(pybel.readfile('xyz', xyzfile))
    except:
        return
    chrg = str(mol.charge)
    os.system('(cd ' + os.path.join(xyzdir, filename) + ' && xtb ' +
              filename + '.xyz' + ' --opt tight --gfn 2 --charge ' + chrg + ' --gbsa benzene --norestart > '
              + 'output_xtb.out)')
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              ' && xtb4stda ' + 'xtbopt.xyz -gbsa benzene > output_xtb4stda.out)')
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              ' && stda -xtb -e 10 > output_stda.out)')
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              ' && grep -A1 state output_stda.out | tail -1 > S1Ex.dat)')
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              ' && stda -xtb -t -e 10 > output_stda_trip.out)')
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              ' && grep -A1 state output_stda_trip.out | tail -1 > T1Ex.dat)')
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              ' && rm wfn.xtb && rm tda.dat)')
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              " && find . -type f -not -name '*.dat' -a -not -name 'ID*.xyz' -exec rm -rf {} \;)")


def compile_xtb_stda(i):
    exDataFile = 'xTB_' + smiles_filename
    if i == 0:
        with open(exDataFile, 'w') as file:
            file.write('SMILES,xtbS1,xtbT1\n')
    else:
        pass
    for val in SMI2ID:
        dir = SMI2ID[val]
        S1file = os.path.join(xyzdir, dir, 'S1Ex.dat')
        T1file = os.path.join(xyzdir, dir, 'T1Ex.dat')
        if not os.path.isfile(S1file) or not os.path.isfile(T1file):
            continue
        with open(S1file, 'r') as file:
            S1data = file.readline()
        with open(T1file, 'r') as file:
            T1data = file.readline()
        try:
            S1Ex = float(S1data.split()[1])
        except IndexError:
            S1Ex = None
        except ValueError:
            S1Ex = None
        try:
            T1Ex = float(T1data.split()[1])
        except IndexError:
            T1Ex = None
        except ValueError:
            T1Ex = None
        if S1Ex is not None and T1Ex is not None:
            with open(exDataFile, 'a') as file:
                file.write(val + ',' + str(S1Ex) + ',' + str(T1Ex) + '\n')


startTime = time.time()
with open('progress.out', 'a') as file:
    file.write('setting up\n')
print('setting up')
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
ephemdir = args.ephemdir
smiles_path = args.smiles_path
smiles_filename = os.path.basename(smiles_path)
smiles_file = os.path.splitext(smiles_filename)[0]
xyzdirname = 'xyzfiles_' + smiles_file
xtb_multi_dirname = os.path.join(ephemdir, 'xtb_multi')
xyzdir = os.path.join(ephemdir, xyzdirname)
if not os.path.isdir(xyzdir):
    os.mkdir(xyzdir)
    with open('current_numrun_' + smiles_file + '.txt', 'w') as file:
        file.write('0')
resultsdirname = xyzdirname.replace('xyzfiles_', 'xtbresults_')
resultsdir = os.path.join(ephemdir, resultsdirname)
if not os.path.isdir(resultsdir):
    os.mkdir(resultsdir)
errFiles = []
logging.root.setLevel(logLevel)
logging.basicConfig(level=logLevel)

with open(smiles_path, 'r') as file:
    data = file.readlines()
smiInd = -1
for index, val in enumerate(data[0].replace('\n', '').split(',')):
    if val.lower() == 'smiles':
        smiInd = index
num_runs = int((len(data) - 1) / 100) + 1
if os.path.isfile('current_numrun_' + smiles_file + '.txt'):
    with open('current_numrun_' + smiles_file + '.txt', 'r') as file:
        num_run = int(file.readline())
else:
    with open('current_numrun_' + smiles_file + '.txt', 'w') as file:
        num_run = 0
        file.write(str(num_run))

print('starting compiling previous xtb-stda')
with open('progress.out', 'a') as file:
    file.write('starting compiling previous xtb-stda\n')
for i in range(num_run):
    SMI2ID = {}
    ID2SMI = {}
    ID = i * 100

    num_added = 0
    for line in data[(ID + 1):]:
        line = line.replace('\n', '')
        line = line.split(',')
        smiles = line[smiInd]
        IDstr = 'ID' + '%09d' % ID
        SMI2ID[smiles] = IDstr
        ID2SMI[IDstr] = smiles
        ID += 1
        num_added += 1
        if num_added >= 100:
            break
    with open('progress.out', 'a') as file:
        file.write('run ' + str(i) + ' : ' + 'compiling data\n')
        file.write('elapsed time: ' + str(time.time() - startTime) + ' seconds\n')
    compile_xtb_stda(i)
print('done compiling previous xtb-stda')
with open('progress.out', 'a') as file:
    file.write('done compiling previous xtb-stda\n')

for i in range(num_run, num_runs):
    with open('current_numrun_' + smiles_file + '.txt', 'w') as file:
        file.write(str(i))
    # ray.shutdown()
    os.system("hostname -I | awk '{print $1}' > ip.txt")
    with open('ip.txt', 'r') as file:
        addressIP = file.readline().replace('\n', '')
    addressIP = addressIP + ":6379"
    with open('progress.out', 'a') as file:
        file.write('initializing ray\n')
    print('initializing ray')
    ray.init(address=addressIP, _redis_password='5241590000000000', _temp_dir='/rds/general/user/sv920/home/tmp/ray')
    # ray.init(address="auto", _redis_password='5241590000000000', _temp_dir='/rds/general/user/sv920/home/tmp/ray')
    with open('progress.out', 'a') as file:
        file.write('still initializing\n')
    availAttempt = 0
    while True:
        time.sleep(5)
        try:
            avail_res = ray.available_resources()
        except:
            availAttempt += 1
            with open('progress.out', 'a') as file:
                file.write('attempt ' + str(availAttempt) + ': re-attempting connection to ray\n')
            if availAttempt >= 10:
                print('ray failed, exiting')
                with open('progress.out', 'a') as file:
                    file.write('ray failed, exiting\n')
                os.chdir('/rds/general/user/sv920/home')
                os.system("find . -mindepth 1 -maxdepth 1 -name '*core*' -exec rm -rf {} \;")
                os.chdir('/rds/general/user/sv920/home/Databases/SCOP_DB/xtb_multi')
                # os.system('qsub job_xtb_run_multi')
                time.sleep(60)
                sys.exit()
            continue
        if avail_res is not None:
            if 'CPU' in avail_res.keys():
                if avail_res['CPU'] >= 32.0:
                    with open('avail_res.json', 'w') as file:
                        json.dump(avail_res, file)
                    break
                else:
                    availAttempt += 1
                    with open('progress.out', 'a') as file:
                        file.write('not enough CPU\n')
            else:
                availAttempt += 1
                with open('progress.out', 'a') as file:
                    file.write('no CPU detected\n')
        else:
            availAttempt += 1
            with open('progress.out', 'a') as file:
                file.write('no resources detected\n')
        if availAttempt >= 10:
            print('ray failed, exiting')
            with open('progress.out', 'a') as file:
                file.write('ray failed, exiting\n')
            time.sleep(60)
            sys.exit()
    with open('ray_info.txt', 'w') as file:
        file.write(str(len(ray.nodes())) + '\n')
        file.write(str(ray.nodes()) + '\n')
    with open('progress.out', 'a') as file:
        file.write('done initializing\n')
    print('done initializing')
    SMI2ID = {}
    ID2SMI = {}
    ID = i * 100
    num_added = 0
    for line in data[(ID + 1):]:
        line = line.replace('\n', '')
        line = line.split(',')
        smiles = line[smiInd]
        IDstr = 'ID' + '%09d' % ID
        SMI2ID[smiles] = IDstr
        ID2SMI[IDstr] = smiles
        ID += 1
        num_added += 1
        if num_added >= 100:
            break

    print('starting generating sdf')
    with open('progress.out', 'a') as file:
        file.write('run ' + str(i) + ' : ' + 'generating sdf files\n')
        file.write('elapsed time: ' + str(time.time() - startTime) + ' seconds\n')
    countAttempt = 0
    while True:
        try:
            ray.get([gen_sdf.remote(x, SMI2ID[x]) for x in SMI2ID], timeout=100)
            break
        except GetTimeoutError:
            if countAttempt >= 10:
                print('ray failed, exiting')
                with open('progress.out', 'a') as file:
                    file.write('ray failed, exiting\n')
                # os.chdir('/rds/general/user/sv920/home')
                # os.system("find . -mindepth 1 -maxdepth 1 -name '*core*' -exec rm -rf {} \;")
                # os.chdir('/rds/general/user/sv920/home/Databases/SCOP_DB/xtb_multi')
                # os.system('qsub job_xtb_run_multi')
                time.sleep(10)
                sys.exit()
            countAttempt += 1
            print('timed out, restarting')
            with open('progress.out', 'a') as file:
                file.write('run ' + str(i) + ', attempt ' + str(countAttempt) + ' : ' + 'timed out, restarting\n')
                file.write('elapsed time: ' + str(time.time() - startTime) + ' seconds\n')
    print('done generating sdf')

    print('starting generating xyz')
    with open('progress.out', 'a') as file:
        file.write('run ' + str(i) + ' : ' + 'generating xyz files\n')
        file.write('elapsed time: ' + str(time.time() - startTime) + ' seconds\n')
    countAttempt = 0
    while True:
        try:
            ray.get([gen_xyz.remote(SMI2ID[x]) for x in SMI2ID], timeout=200)
            break
        except GetTimeoutError:
            if countAttempt >= 10:
                print('ray failed, exiting')
                with open('progress.out', 'a') as file:
                    file.write('ray failed, exiting\n')
                # os.chdir('/rds/general/user/sv920/home')
                # os.system("find . -mindepth 1 -maxdepth 1 -name '*core*' -exec rm -rf {} \;")
                # os.chdir('/rds/general/user/sv920/home/Databases/SCOP_DB/xtb_multi')
                # os.system('qsub job_xtb_run_multi')
                time.sleep(10)
                sys.exit()
            countAttempt += 1
            print('timed out, restarting')
            with open('progress.out', 'a') as file:
                file.write('run ' + str(i) + ', attempt ' + str(countAttempt) + ' : ' + 'timed out, restarting\n')
                file.write('elapsed time: ' + str(time.time() - startTime) + ' seconds\n')
    print('done generating xyz')

    print('starting running xtb-stda')
    with open('progress.out', 'a') as file:
        file.write('run ' + str(i) + ' : ' + 'running xtb\n')
        file.write('elapsed time: ' + str(time.time() - startTime) + ' seconds\n')
    countAttempt = 0
    while True:
        try:
            ray.get([run_xtb_stda.remote(SMI2ID[x]) for x in SMI2ID], timeout=300)
            break
        except GetTimeoutError:
            if countAttempt >= 10:
                print('ray failed, exiting')
                with open('progress.out', 'a') as file:
                    file.write('ray failed, exiting\n')
                os.chdir('/rds/general/user/sv920/home')
                os.system("find . -mindepth 1 -maxdepth 1 -name '*core*' -exec rm -rf {} \;")
                os.chdir('/rds/general/user/sv920/home/Databases/SCOP_DB/xtb_multi')
                # os.system('qsub job_xtb_run_multi')
                time.sleep(10)
                sys.exit()
            countAttempt += 1
            print('timed out, restarting')
            with open('progress.out', 'a') as file:
                file.write('run ' + str(i) + ', attempt ' + str(countAttempt) + ' : ' + 'timed out, restarting\n')
                file.write('elapsed time: ' + str(time.time() - startTime) + ' seconds\n')
    print('done running xtb-stda')

    print('starting compiling xtb-stda')
    with open('progress.out', 'a') as file:
        file.write('run ' + str(i) + ' : ' + 'compiling data\n')
        file.write('elapsed time: ' + str(time.time() - startTime) + ' seconds\n')
    compile_xtb_stda(i)
    print('done compiling xtb-stda')

    with open('progress.out', 'a') as file:
        file.write('run ' + str(i) + ' : ' + 'successfully completed\n')
        file.write('in ' + str(time.time() - startTime) + ' seconds\n')

    with open('smiles_filename.txt', 'w') as file:
        file.write(smiles_filename)

    print('successfully completed!')
    ray.shutdown()

with open('current_numrun_' + smiles_file + '.txt', 'w') as file:
    file.write(str(i))

with open('progress.out', 'a') as file:
    file.write('successfully completed all ' + str(len(data) - 1) + ' molecules\n')
    file.write('in ' + str(time.time() - startTime) + ' seconds\n')
