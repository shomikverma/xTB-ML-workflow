import matplotlib.pyplot as plt
import numpy as np
import os
import json
import sys
from tqdm import tqdm
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn import linear_model
from cycler import cycler
import string
from itertools import cycle
import pandas as pd


def label_axes(fig, labels=None, loc=None, **kwargs):
    if labels is None:
        labels = string.ascii_lowercase
    labels = cycle(labels)
    if loc is None:
        loc = (-0.1, 1.1)
    axes = [ax for ax in fig.axes if ax.get_label() != '<colorbar>']
    for ax, lab in zip(axes, labels):
        ax.annotate('(' + lab + ')', size=14, xy=loc,
                    xycoords='axes fraction',
                    **kwargs)


plt.style.use(['science', 'grid'])
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colors = [colors[0], colors[2], colors[1]] + colors[3:]
plt.rcParams['axes.prop_cycle'] = cycler(color=colors)


def gen_train_data(runType):
    with open('../xtb_results_qmsym_10k/QMsym_10k_xtb.csv', 'r') as file:
        dataQMsym = pd.read_csv(file)

    with open('../xtb_results_TTA_SF/TTA_SF_xtb.csv', 'r') as file:
        dataSCOP = pd.read_csv(file)

    with open('../xtb_results_AL_PCQC_T1/AL_PCQC_T1_0_xtb.csv', 'r') as file:
        dataALT1 = pd.read_csv(file)

    with open('../xtb_results_AL_PCQC_S1/AL_PCQC_S1_xtb.csv', 'r') as file:
        dataALS1 = pd.read_csv(file)

    with open('../xtb_results_qmsym/QMsymex_xtb.csv', 'r') as file:
        dataQMsymex = pd.read_csv(file)

    trainData = pd.DataFrame()
    if runType == '':
        trainData = dataQMsym.merge(dataSCOP, how='outer')
    elif runType == 'scop':
        trainData = dataSCOP
    elif runType == 'scop_qm':
        trainData = dataSCOP.merge(dataQMsym, how='outer')
    elif runType == 'scop_AL':
        trainData = dataSCOP.merge(dataALT1, how='outer').merge(dataALS1, how='outer')
    elif runType == 'scop_qm_AL':
        trainData = dataQMsym.merge(dataSCOP, how='outer').merge(dataALT1, how='outer').merge(dataALS1, how='outer')
    elif runType == 'scop_qm_AL_ex':
        trainData = dataQMsym.merge(dataSCOP, how='outer').merge(dataALT1, how='outer').merge(dataALS1, how='outer'). \
            merge(dataQMsymex, how='outer')
    elif runType == '300k':
        trainData = dataQMsym.merge(dataSCOP, how='outer').merge(dataALT1, how='outer').merge(dataALS1, how='outer'). \
            merge(dataQMsymex, how='outer')

    trainData.to_csv('all_train_data_' + runType + '.csv', columns=['SMILES', 'T1err', 'S1err'], index=False)
    trainData.to_csv('all_featu_data_' + runType + '.csv', columns=['xtb_T1', 'xtb_S1'], index=False)


def train_ML(runType):
    os.system('chemprop_train --data_path all_train_data_' + runType +
              '.csv --features_path all_featu_data_' + runType +
              '.csv --dataset_type regression --save_dir xTB_ML_model_nrg_' + runType +
              ' --split_type cv --save_smiles_splits --num_folds 10 --target_columns S1err T1err')


def gen_test_data(testType):
    with open('xTB_' + testType + '.csv', 'r') as file:
        dataXTB = pd.read_csv(file)
    dataXTB.rename({'xtbT1': 'xtb_T1', 'xtbS1': 'xtb_S1'}, axis='columns', inplace=True)
    dataXTB.to_csv('all_preds_smiles_' + testType + '.csv', columns=['SMILES'], index=False)
    dataXTB.to_csv('all_preds_featur_' + testType + '.csv', columns=['xtb_T1', 'xtb_S1'], index=False)


def pred_ML(runType, testType):
    if os.path.isfile(testType + '_preds_nrg_' + runType + '.csv'):
        return
    os.system('chemprop_predict --test_path all_preds_smiles_' + testType + '.csv --features_path all_preds_featur_'
              + testType + '.csv --checkpoint_dir xTB_ML_model_nrg_' + runType
              + ' --preds_path ' + testType + '_preds_nrg_' + runType + '.csv --drop_extra_columns')


def analyze_ML(runType, testType):
    with open('xTB_' + testType + '.csv', 'r') as file:
        dataXTB = pd.read_csv(file)
    with open(testType + '_preds_nrg_' + runType + '.csv', 'r') as file:
        dataPreds = pd.read_csv(file)
    with open('TDDFT_' + testType + '.csv', 'r') as file:
        dataTDDFT = pd.read_csv(file)
        dataTDDFT.rename({'smiles': 'SMILES'}, axis='columns', inplace=True)

    dataAll = dataXTB.merge(dataPreds).merge(dataTDDFT)
    dataAll.rename({'S1': 'TDDFT_S1', 'T1': 'TDDFT_T1'}, axis='columns', inplace=True)

    # ML calib
    dataAll['xTB_ML_S1'] = dataAll['xtbS1'] + dataAll['S1err']
    dataAll['xTB_ML_T1'] = dataAll['xtbT1'] + dataAll['T1err']

    # lin calib
    regr = linear_model.LinearRegression()
    regr.fit(np.array(dataAll['xtbS1']).reshape(-1, 1), np.array(dataAll['TDDFT_S1']).reshape(-1, 1))
    dataAll['xTB_Lin_S1'] = [x[0] for x in regr.predict(np.array(dataAll['xtbS1']).reshape(-1, 1))]
    regr = linear_model.LinearRegression()
    regr.fit(np.array(dataAll['xtbT1']).reshape(-1, 1), np.array(dataAll['TDDFT_T1']).reshape(-1, 1))
    dataAll['xTB_Lin_T1'] = [x[0] for x in regr.predict(np.array(dataAll['xtbT1']).reshape(-1, 1))]
    # dataAll.info()
    dataAll.to_csv('all_data_' + runType + '_' + testType + '.csv', index=False)

    # get metrics
    r2_S1 = r2_score(dataAll['TDDFT_S1'], dataAll['xtbS1'])
    r2_T1 = r2_score(dataAll['TDDFT_S1'], dataAll['xtbT1'])
    MAE_S1 = mean_absolute_error(dataAll['TDDFT_S1'], dataAll['xtbS1'])
    MAE_T1 = mean_absolute_error(dataAll['TDDFT_T1'], dataAll['xtbT1'])
    r2_lin_S1 = r2_score(dataAll['TDDFT_S1'], dataAll['xTB_Lin_S1'])
    r2_lin_T1 = r2_score(dataAll['TDDFT_T1'], dataAll['xTB_Lin_T1'])
    MAE_lin_S1 = mean_absolute_error(dataAll['TDDFT_S1'], dataAll['xTB_Lin_S1'])
    MAE_lin_T1 = mean_absolute_error(dataAll['TDDFT_T1'], dataAll['xTB_Lin_T1'])
    r2_ML_S1 = r2_score(dataAll['TDDFT_S1'], dataAll['xTB_ML_S1'])
    r2_ML_T1 = r2_score(dataAll['TDDFT_T1'], dataAll['xTB_ML_T1'])
    MAE_ML_S1 = mean_absolute_error(dataAll['TDDFT_S1'], dataAll['xTB_ML_S1'])
    MAE_ML_T1 = mean_absolute_error(dataAll['TDDFT_T1'], dataAll['xTB_ML_T1'])
    print(testType, runType)
    print('stda')
    print(r2_S1, r2_T1)
    print(MAE_S1, MAE_T1)
    print('lin')
    print(r2_lin_S1, r2_lin_T1)
    print(MAE_lin_S1, MAE_lin_T1)
    print('ML')
    print(r2_ML_S1, r2_ML_T1)
    print(MAE_ML_S1, MAE_ML_T1)

    fig = plt.figure(num=1, figsize=[7, 4], dpi=300, clear=True)
    ax = fig.add_subplot(1, 2, 1)
    ax.set_axisbelow(True)
    if testType == 'AL_PCQC_T1_xTB_test':
        plt.plot(dataAll['xtbS1'], dataAll['TDDFT_S1'], '.', color=colors[3], markersize=0.1, label='orig')
        plt.plot(dataAll['xTB_Lin_S1'], dataAll['TDDFT_S1'], '.', color=colors[2], markersize=0.1, label='lin calib')
        plt.plot(dataAll['xTB_ML_S1'], dataAll['TDDFT_S1'], '.', color=colors[0], markersize=0.1, label='ML calib')
    else:
        plt.plot(dataAll['xtbS1'], dataAll['TDDFT_S1'], '.', color=colors[3], label='orig')
        plt.plot(dataAll['xTB_Lin_S1'], dataAll['TDDFT_S1'], '.', color=colors[2], label='lin calib')
        plt.plot(dataAll['xTB_ML_S1'], dataAll['TDDFT_S1'], '.', color=colors[0], label='ML calib')
    x = np.linspace(0, 9, 100)
    plt.plot(x, x, 'k--')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.xlabel('stda S1 (eV)')
    plt.ylabel('TDDFT S1 (eV)')
    plt.title('TD-DFT vs. stda S1 comparison')
    plt.xlim(0, 9)
    plt.ylim(0, 9)
    plt.annotate('R2 orig: %0.2f\n' % r2_S1 +
                 'R2 lin: %0.2f\n' % r2_lin_S1 +
                 'R2 ML: %0.2f\n' % r2_ML_S1 +
                 'MAE orig: %0.2f\n' % MAE_S1 +
                 'MAE lin: %0.2f\n' % MAE_lin_S1 +
                 'MAE ML: %0.2f' % MAE_ML_S1,
                 (8.5, 0.5),
                 bbox=dict(facecolor='white', alpha=0.5),
                 ha='right')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()

    ax = fig.add_subplot(1, 2, 2)
    ax.set_axisbelow(True)
    if testType == 'AL_PCQC_T1_xTB_test':
        plt.plot(dataAll['xtbT1'], dataAll['TDDFT_T1'], '.', color=colors[3], markersize=0.1, label='orig')
        plt.plot(dataAll['xTB_Lin_T1'], dataAll['TDDFT_T1'], '.', color=colors[2], markersize=0.1, label='lin calib')
        plt.plot(dataAll['xTB_ML_T1'], dataAll['TDDFT_T1'], '.', color=colors[0], markersize=0.1, label='ML calib')
    else:
        plt.plot(dataAll['xtbT1'], dataAll['TDDFT_T1'], '.', color=colors[3], label='orig')
        plt.plot(dataAll['xTB_Lin_T1'], dataAll['TDDFT_T1'], '.', color=colors[2], label='lin calib')
        plt.plot(dataAll['xTB_ML_T1'], dataAll['TDDFT_T1'], '.', color=colors[0], label='ML calib')
    x = np.linspace(0, 9, 100)
    plt.plot(x, x, 'k--')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.xlabel('stda T1 (eV)')
    plt.ylabel('TDDFT T1 (eV)')
    plt.title('TD-DFT vs. stda T1 comparison')
    plt.xlim(0, 9)
    plt.ylim(0, 9)
    plt.annotate('R2 orig: %0.2f\n' % r2_T1 +
                 'R2 lin: %0.2f\n' % r2_lin_T1 +
                 'R2 ML: %0.2f\n' % r2_ML_T1 +
                 'MAE orig: %0.2f\n' % MAE_T1 +
                 'MAE lin: %0.2f\n' % MAE_lin_T1 +
                 'MAE ML: %0.2f' % MAE_ML_T1,
                 (8.5, 0.5),
                 bbox=dict(facecolor='white', alpha=0.5),
                 ha='right')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    label_axes(fig, ha='left')
    plt.savefig('mopssam_calib_nrg_' + runType + '_' + testType + '.png')


# runType = ''
# testType = 'MOPSSAM_143'
# for runType in ['scop', 'scop_qm', 'scop_AL', 'scop_qm_AL', 'scop_qm_AL_ex']:
for runType in ['scop_qm', 'scop_qm_AL_ex']:
    # gen_train_data(runType)
    # train_ML(runType)
    for testType in ['MOPSSAM_143', 'MOPSSAM_1k', 'verde_smiles', 'INDT_smiles', 'AL_PCQC_T1_xTB_test']:
        # for testType in ['AL_PCQC_T1_xTB_test']:
        gen_test_data(testType)
        pred_ML(runType, testType)
        analyze_ML(runType, testType)
