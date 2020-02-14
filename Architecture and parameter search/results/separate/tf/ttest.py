import numpy as np
import pandas as pd
import sys
from scipy import stats
import re


def ttest(tMetric, pMetric, metric):
    t, p = stats.ttest_ind(tMetric, pMetric)
    print("\n\n\#\#\#\#\#\#")
    print('t-value for ' + metric + ' = ' + str(t))
    print('p-value for ' + metric + ' = ' + str(p))
    if t > p:
        print('Thus, '+ metric +' is statistically Significant')
    else:
        print('Thus, '+ metric +' is statistically INsignificant')

def fixDfToPd(df, col):
    tempList = df.values.tolist()
    newList = []
    for j in tempList:
        newList.append(float(re.sub("[^0-9\.]", "", j)))
    return pd.DataFrame(newList).to_numpy()

def main():

    startEpoch = 1

    deepPsychName = sys.argv[1]
    deepSeaName = sys.argv[2]
    epochs = int(sys.argv[3])

    trainLoss = []
    trainTpr= []
    trainTnr = []
    trainUncer = []
    trainRoc = []
    trainAcc = []

    testLoss = []
    testTpr= []
    testTnr = []
    testUncer = []
    testRoc = []
    testAcc = []

    print("running T-Test between " + deepPsychName.upper() + ' and ' + deepSeaName.upper() + '\n\n')
    
    dataDeepPysch = pd.read_csv(deepPsychName, names=[ 'Epoch','Batch','Loss','Uncertainty',
                                            'nOnesTarget','nOnesPred',
                                            'tpr', 'tnr','roc','acc'])

    dataDeepSea = pd.read_csv(deepSeaName, names=[ 'Epoch','Batch','Loss','Uncertainty',
                                            'nOnesTarget','nOnesPred',
                                            'tpr', 'tnr','roc','acc'])
                                            
    dataDeepPysch['Epoch'] = dataDeepPysch['Epoch'].str.split('|').str[1]
    dataDeepSea['Epoch'] = dataDeepSea['Epoch'].str.split('|').str[1]


    testDeepPsychDf = dataDeepPysch[dataDeepPysch['Epoch'].str.contains('VALIDATION')]
    testDeepSeaDf = dataDeepSea[dataDeepSea['Epoch'].str.contains('VALIDATION')]

    for i in range(startEpoch, epochs+1):
        if i == epochs:
            epochDeepPsychTestDf = testDeepPsychDf[testDeepPsychDf['Epoch'].str.endswith('Epoch: '+str(i))]
            epochDeepSeaTestDf = testDeepSeaDf[testDeepSeaDf['Epoch'].str.endswith('Epoch: '+str(i))]
            
            # Testing Metrics
            dpLoss = fixDfToPd(epochDeepPsychTestDf['Loss'], 'Loss')
            dcLoss = fixDfToPd(epochDeepSeaTestDf['Loss'], 'Loss') 

            dpTpr = fixDfToPd(epochDeepPsychTestDf['tpr'], 'tpr') 
            dcTpr = fixDfToPd(epochDeepSeaTestDf['tpr'], 'tpr')

            dpAcc = fixDfToPd(epochDeepPsychTestDf['acc'], 'acc') 
            dcAcc = fixDfToPd(epochDeepSeaTestDf['acc'], 'acc')

            dpCert = fixDfToPd(epochDeepPsychTestDf['Uncertainty'], 'Uncertainty')
            dcCert = fixDfToPd(epochDeepSeaTestDf['Uncertainty'], 'Uncertainty')

            dpRoc = fixDfToPd(epochDeepPsychTestDf['roc'], 'roc')
            dcRoc = fixDfToPd(epochDeepSeaTestDf['roc'], 'roc')
            
        i += 1

    ttest(dpLoss, dcLoss, 'LOSS')
    ttest(dpTpr, dcTpr, 'TPR')
    ttest(dpTpr, dcTpr, 'TPR')
    ttest(dpCert, dcCert, 'CERTAINTY')
    ttest(dpRoc, dcRoc, 'ROC')

main()