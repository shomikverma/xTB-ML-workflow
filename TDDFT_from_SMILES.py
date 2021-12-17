print('importing')
from openbabel import pybel
from openbabel import openbabel as ob
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
import ray
import time
import sys
import json
import re

f = open('progress.out', 'a')
sys.stdout = f


# utils
def shell(cmd, shell=False):
    if shell:
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = cmd.split()
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate()
    return output


# generate 3D coordinates with openbabel
def gen_sdf(smiles, filename):
    if os.path.isfile(os.path.join(xyzdir, filename, filename + '.xyz')):
        return
    if os.path.isfile(os.path.join(xyzdir, filename + '.sdf')):
        return
    mol2 = pybel.readstring('smi', smiles)
    gen3d = ob.OBOp.FindType("gen3D")
    gen3d.Do(mol2.OBMol, "--medium")
    mol2.write('sdf', filename=os.path.join(xyzdir, filename + '.sdf'), overwrite=True)
    if not os.path.isfile(os.path.join(xyzdir, filename + '.sdf')):
        mol2.write('sdf', filename=os.path.join(xyzdir, filename + '.sdf'), overwrite=True)
    with open(os.path.join(xyzdir, filename + '.sdf'), 'r') as file:
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
        mol2.write('sdf', filename=os.path.join(xyzdir, filename + '.sdf'), overwrite=True)
    with open(os.path.join(xyzdir, filename + '.sdf'), 'r') as file:
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
        mol2.write('sdf', filename=os.path.join(xyzdir, filename + '.sdf'), overwrite=True)


# optimize initial 3D coordinates with xTB
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


# generate S0 comfiles for DFT
def gen_S0_comfiles(smiles, filename):
    # if os.path.isfile(os.path.join(xyzdir, filename, filename + '_S0.com')):
    #     return
    with open(os.path.join(xyzdir, filename, 'smiles'), 'w') as file:
        file.write(smiles + '\n')
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              " && obabel -ixyz " + filename + ".xyz -O " + filename + "_prelim.com)")
    with open(os.path.join(xyzdir, filename, filename + '_prelim.com'), 'r') as file:
        comdata = file.readlines()
    with open(os.path.join(xyzdir, filename, filename + '_S0.com'), 'w') as file:
        file.write("%chk=gauss_S0.chk\n")
        file.write("%mem=100GB\n")
        file.write("%nproc=48\n")
        file.write("#b3lyp 6-31G* opt \n\n")
        file.write(' ' + filename + "_S0\n\n")
        mol = Chem.MolFromSmiles(smiles)
        chrg = Chem.GetFormalCharge(mol)
        rads = Descriptors.NumRadicalElectrons(mol)
        if rads % 2 == 0:
            mult = 1
        else:
            mult = 2
        file.write(str(chrg) + ' ' + str(mult) + '\n')
        # f.write("0 1\n")
        writeLine = False
        for index, line in enumerate(comdata):
            if index <= 5:
                continue
            file.write(line)
            # if writeLine:
            #     file.write(line)
            # if '0  1' in line:
            #     writeLine = True
        file.write('\n')
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              " && rm " + filename + "_prelim.com)")


# generate S0 jobfile
# autosubmit jobfile

# run S0 DFT
def run_S0_DFT(filename):
    if os.path.isfile(os.path.join(xyzdir, filename, filename + '_S0.log')):
        with open(os.path.join(xyzdir, filename, filename + '_S0.log'), 'r') as file:
            data = file.readlines()
        if 'Normal termination' in data[-1]:
            pass
            # return
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              " && g16 " + filename + "_S0.com)")


# get 3D coordinates from DFT output
def gen_coord_for_S1T1(filename):
    if not os.path.isfile(os.path.join(xyzdir, filename, filename + '_S0.log')):
        print('no log file for ' + filename)
        return
    # if os.path.isfile(os.path.join(xyzdir, filename, filename + '_S0.xyz')):
    #     return
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              " && g092xyz.pl " + filename + "_S0.log)")
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              " && mv g09-result.xyz " + filename + "_S0.xyz)")


# generate S1,T1 comfiles for TD-DFT
def gen_S1T1_comfiles(smiles, filename):
    # if os.path.isfile(os.path.join(xyzdir, filename, filename + '_S1.com')) and \
    #         os.path.isfile(os.path.join(xyzdir, filename, filename + '_T1.com')):
    #     return
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              " && obabel -ixyz " + filename + "_S0.xyz -O " + filename + "_prelim.com)")
    with open(os.path.join(xyzdir, filename, filename + '_prelim.com'), 'r') as file:
        comdata = file.readlines()
    with open(os.path.join(xyzdir, filename, filename + '_S1.com'), 'w') as file:
        file.write("%chk=gauss_S1.chk\n")
        file.write("%mem=100GB\n")
        file.write("%nproc=48\n")
        file.write("#b3lyp/6-31+g*  units=ang scf=(tight) nosym "
                   "td(nstates=10) int=(grid=ultrafine) \n\n")
        file.write(' ' + filename + "_S1\n\n")
        mol = Chem.MolFromSmiles(smiles)
        chrg = Chem.GetFormalCharge(mol)
        rads = Descriptors.NumRadicalElectrons(mol)
        if rads % 2 == 0:
            mult = 1
        else:
            mult = 2
        file.write(str(chrg) + ' ' + str(mult) + '\n')
        for index, line in enumerate(comdata):
            if index <= 5:
                continue
            file.write(line)
        # writeLine = False
        # for line in comdata:
        #     if writeLine:
        #         file.write(line)
        #     if '0  1' in line:
        #         writeLine = True
        file.write('\n')
    with open(os.path.join(xyzdir, filename, filename + '_T1.com'), 'w') as file:
        file.write("%chk=gauss_T1.chk\n")
        file.write("%mem=100GB\n")
        file.write("%nproc=48\n")
        file.write("#b3lyp/6-31+g*  units=ang scf=(tight) nosym "
                   "td(triplets,nstates=10) int=(grid=ultrafine) \n\n")
        file.write(' ' + filename + "_T1\n\n")
        mol = Chem.MolFromSmiles(smiles)
        chrg = Chem.GetFormalCharge(mol)
        rads = Descriptors.NumRadicalElectrons(mol)
        if rads % 2 == 0:
            mult = 1
        else:
            mult = 2
        file.write(str(chrg) + ' ' + str(mult) + '\n')
        for index, line in enumerate(comdata):
            if index <= 5:
                continue
            file.write(line)
        # writeLine = False
        # for line in comdata:
        #     if writeLine:
        #         file.write(line)
        #     if '0  1' in line:
        #         writeLine = True
        file.write('\n')
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              " && rm " + filename + "_prelim.com)")


# generate S1,T1 comfiles
# autosubmit comfiles

# run S1 TDDFT
def run_S1_TDDFT(filename):
    if os.path.isfile(os.path.join(xyzdir, filename, filename + '_S1.log')):
        with open(os.path.join(xyzdir, filename, filename + '_S1.log'), 'r') as file:
            data = file.readlines()
        if 'Normal termination' in data[-1]:
            # pass
            return
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              " && g16 " + filename + "_S1.com)")


# run T1 TDDFT
def run_T1_TDDFT(filename):
    if os.path.isfile(os.path.join(xyzdir, filename, filename + '_T1.log')):
        with open(os.path.join(xyzdir, filename, filename + '_T1.log'), 'r') as file:
            data = file.readlines()
        if 'Normal termination' in data[-1]:
            # pass
            return
    os.system('(cd ' + os.path.join(xyzdir, filename) +
              " && g16 " + filename + "_T1.com)")


# compile results
def extract_data(filename):
    if os.path.isfile(os.path.join(xyzdir, filename, filename + '_S1.log')) and \
            os.path.isfile(os.path.join(xyzdir, filename, filename + '_T1.log')):
        os.system('(cd ' + os.path.join(xyzdir, filename) +
                  " && grep 'Excited State   1' " + filename + "_S1.log | tail -1 > excitedstate_S1.out)")
        os.system('(cd ' + os.path.join(xyzdir, filename) +
                  " && grep 'Excited State   1' " + filename + "_T1.log | tail -1 > excitedstate_T1.out)")
    else:
        with open('excitedstate_S1.out', 'w') as file:
            file.write('')
        with open('excitedstate_T1.out', 'w') as file:
            file.write('')


def compile_data():
    allData = {}
    for smiles in SMI2ID:
        filename = SMI2ID[smiles]
        if os.path.isfile(os.path.join(xyzdir, filename, 'excitedstate_S1.out')) and \
                os.path.isfile(os.path.join(xyzdir, filename, 'excitedstate_T1.out')):
            with open(os.path.join(xyzdir, filename, 'excitedstate_S1.out'), 'r') as file:
                S1data = file.readline()
            with open(os.path.join(xyzdir, filename, 'excitedstate_T1.out'), 'r') as file:
                T1data = file.readline()
            try:
                S1Ex = float(re.search('(.+?)eV', S1data).group(1).split()[-1])
            except:
                continue
            try:
                T1Ex = float(re.search('(.+?)eV', T1data).group(1).split()[-1])
            except:
                continue
            if S1Ex is not None and T1Ex is not None:
                if 0 < S1Ex < 10 and 0 < T1Ex < 10:
                    allData[smiles] = {}
                    allData[smiles]['S1'] = S1Ex
                    allData[smiles]['T1'] = T1Ex
    with open('../TDDFT_' + smiles_file + '.csv', 'w') as file:
        file.write('SMILES,S1,T1\n')
        for smiles in allData:
            file.write(smiles + ',' + str(allData[smiles]['S1']) + ',' + str(allData[smiles]['T1']) + '\n')


print('setting up')
# smiles_path = '/rds/general/user/sv920/home/Databases/SCOP_DB/TDDFT_auto/pot_sens_mopssam.csv'
smiles_path = sys.argv[2]
basedir = '/rds/general/user/sv920/ephemeral'
smiles_filename = os.path.basename(smiles_path)
smiles_file = os.path.splitext(smiles_filename)[0]
xyzdir = os.path.join(basedir, 'xyzfiles_TDDFT_' + smiles_file)
comdir = os.path.join(basedir, 'comfiles_TDDFT_' + smiles_file)
if not os.path.isdir(xyzdir):
    os.mkdir(xyzdir)
if not os.path.isdir(comdir):
    os.mkdir(comdir)

with open(smiles_path, 'r') as file:
    data = file.readlines()
smiInd = -1
for index, val in enumerate(data[0].replace('\n', '').split(',')):
    if val.lower() == 'smiles':
        smiInd = index

print('loading data')
SMI2ID = {}
ID2SMI = {}
ID = 0
for line in data:
    if 'smiles' in line:
        continue
    line = line.replace('\n', '')
    line = line.split(',')
    smiles = line[smiInd]
    IDstr = 'ID' + '%09d' % ID
    SMI2ID[smiles] = IDstr
    ID2SMI[IDstr] = smiles
    ID += 1

doSerial = False
doRay = False
if doSerial:
    for smiles in SMI2ID:
        filename = SMI2ID[smiles]
        print(filename)
        print('generating sdf')
        gen_sdf(smiles, filename)
        print('generating xyz')
        gen_xyz(filename)
        print('generating S0 comfile')
        gen_S0_comfiles(smiles, filename)
        print('running ground state DFT')
        # run_S0_DFT(filename)
        print('getting ground state coords')
        gen_coord_for_S1T1(filename)
        print('generating S1/T1 comfiles')
        gen_S1T1_comfiles(smiles, filename)
        print('running S1 TDDFT')
        # run_S1_TDDFT(filename)
        print('running T1 TDDFT')
        # run_T1_TDDFT(filename)
        print('extracting data')
        extract_data(filename)


def auto_TDDFT(smiles, filename):
    print('starting calcs for ' + filename)
    with open('progress.out', 'a') as file:
        file.write('starting calcs for ' + filename + '\n')
    print('generating sdf')
    gen_sdf(smiles, filename)
    print('generating xyz')
    gen_xyz(filename)
    print('generating S0 comfile')
    gen_S0_comfiles(smiles, filename)
    print('running ground state DFT')
    run_S0_DFT(filename)
    print('getting ground state coords')
    gen_coord_for_S1T1(filename)
    print('generating S1/T1 comfiles')
    gen_S1T1_comfiles(smiles, filename)
    print('running S1 TDDFT')
    run_S1_TDDFT(filename)
    print('running T1 TDDFT')
    run_T1_TDDFT(filename)
    print('extracting data')
    extract_data(filename)
    with open('progress.out', 'a') as file:
        file.write('finished calcs for ' + filename + '\n')


@ray.remote(num_cpus=40)
def parallel_TDDFT(smiles, filename):
    print('starting calcs for ' + filename)
    with open('progress.out', 'a') as file:
        file.write('starting calcs for ' + filename + '\n')
    print('generating sdf')
    gen_sdf(smiles, filename)
    print('generating xyz')
    gen_xyz(filename)
    print('generating S0 comfile')
    gen_S0_comfiles(smiles, filename)
    print('running ground state DFT')
    run_S0_DFT(filename)
    print('getting ground state coords')
    gen_coord_for_S1T1(filename)
    print('generating S1/T1 comfiles')
    gen_S1T1_comfiles(smiles, filename)
    print('running S1 TDDFT')
    run_S1_TDDFT(filename)
    print('running T1 TDDFT')
    run_T1_TDDFT(filename)
    print('extracting data')
    extract_data(filename)
    with open('progress.out', 'a') as file:
        file.write('finished calcs for ' + filename + '\n')


if not doSerial:
    index = int(sys.argv[1])
    smiles = list(SMI2ID.keys())[index]
    filename = SMI2ID[smiles]
    auto_TDDFT(smiles, filename)

if doRay:
    print('setting up ray')
    os.system("hostname -I | awk '{print $1}' > ip.txt")
    with open('ip.txt', 'r') as file:
        addressIP = file.readline().replace('\n', '')
    addressIP = addressIP + ":6379"
    print('initializing ray')
    with open('progress.out', 'a') as file:
        file.write('initializing ray\n')
    time.sleep(60)
    # ray.init(address=addressIP, _redis_password='5241590000000000',
    #          _temp_dir='/rds/general/user/sv920/home/tmp/ray')
    ray.init(address="auto", _redis_password='5241590000000000', _temp_dir='/rds/general/user/sv920/home/tmp/ray')
    print('checking initialization')
    with open('progress.out', 'a') as file:
        file.write('checking initialization\n')
    availAttempt = 0
    while True:
        time.sleep(5)
        try:
            avail_res = ray.available_resources()
        except:
            availAttempt += 1
            print('attempt ' + str(availAttempt) + ': re-attempting connection to ray')
            with open('progress.out', 'a') as file:
                file.write('attempt ' + str(availAttempt) + ': re-attempting connection to ray\n')
            continue
        if avail_res is not None:
            if 'CPU' in avail_res.keys():
                # if avail_res['CPU'] >= 36.0:
                with open('avail_res.json', 'w') as file:
                    json.dump(avail_res, file)
                break
    with open('ray_info.txt', 'w') as file:
        file.write(str(len(ray.nodes())) + '\n')
        file.write(str(ray.nodes()) + '\n')
    print('done initializing')
    time.sleep(60)
    print('starting TDDFT calculations')
    with open('progress.out', 'a') as file:
        file.write('starting TDDFT calculations\n')
    rs = []
    for x in SMI2ID:
        rs.append(parallel_TDDFT.remote(x, SMI2ID[x]))
        time.sleep(5)
    # ray.get([parallel_TDDFT.remote(x, SMI2ID[x]) for x in SMI2ID])
    ray.get(rs)

compile_data()

print('successfully completed!')
with open('progress.out', 'a') as file:
    file.write('successfully completed!\n')
f.close()
