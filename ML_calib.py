import os
import time
import numpy as np

with open('smiles_filename.txt', 'r') as file:
    smiles_filename = file.readline().replace('\n', '')


def ML_calib():
    os.system('chemprop_predict --test_path xTB_' + smiles_filename +
              ' --checkpoint_dir /rds/general/user/sv920/home/Databases/SCOP_DB/QMsym_SCOP_ML '
              '--preds_path xTB_ML_errs_' + smiles_filename + ' > chemprop_progress.out 2>&1')
    with open('xTB_ML_errs_' + smiles_filename, 'r') as file:
        errData = file.readlines()
    xtbS1Ind = -1
    xtbT1Ind = -1
    S1errInd = -1
    T1errInd = -1
    line = errData[0]
    line = line.replace('\n', '')
    line = line.split(',')
    for index, val in enumerate(line):
        if val == 'xtbS1':
            xtbS1Ind = index
        if val == 'xtbT1':
            xtbT1Ind = index
        if val == 'S1err':
            S1errInd = index
        if val == 'T1err':
            T1errInd = index
    with open('xTB_ML_calib_' + smiles_filename, 'w') as file:
        file.write('SMILES,S1_xTB_ML,T1_xTB_ML,Sens_FOM,Emit_FOM,Sens_good,Emit_good\n')
        for line in errData[1:]:
            line = line.replace('\n', '')
            line = line.split(',')
            smiles = line[0]
            S1calib = float(line[xtbS1Ind]) + float(line[S1errInd])
            T1calib = float(line[xtbT1Ind]) + float(line[T1errInd])
            if T1calib == 0 or S1calib == 0:
                continue
            Sens_FOM = 0 if T1calib > S1calib else np.exp(-(1 - T1calib / S1calib))
            Emit_FOM = 0 if S1calib > 2 * T1calib else np.exp(-(2 - S1calib / T1calib))
            Sens_Good = 1 if 0.95 < T1calib / S1calib < 1 else 0
            Emit_Good = 1 if 1.9 < S1calib / T1calib < 2 else 0
            file.write(smiles + ',' +
                       str(S1calib) + ',' + str(T1calib) + ',' +
                       str(Sens_FOM) + ',' + str(Emit_FOM) + ',' +
                       str(Sens_Good) + ',' + str(Emit_Good) + '\n')


with open('progress_ML.out', 'w') as file:
    file.write('calibrating with ML\n')

startTime = time.time()
ML_calib()

with open('progress_ML.out', 'a') as file:
    file.write('ML calib completed\n')
    file.write('in ' + str(time.time() - startTime) + ' seconds\n')

with open('done_xTB_ML_calib_' + smiles_filename, 'w') as file:
    file.write('done')
