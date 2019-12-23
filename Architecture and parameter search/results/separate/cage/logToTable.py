import pandas as pd 
import sys
import statistics

'''
	python logToTable.py FILE.log numOfTotalEpochsInLogFile
'''

def calcMedian(df, col):

	list = df[col].values.tolist()
	if col in ['Loss','Uncertainty']:
		list = [float(i.split(' ', 3)[2]) for i in list]
	else:
		list = [float(i.split(' ', 4)[3]) for i in list]
	return statistics.median(list)

startEpoch = 1

filename = sys.argv[1]
epochs = int(sys.argv[2])

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

data = pd.read_csv(sys.argv[1], names=[ 'Epoch','Batch','Loss','Uncertainty',
										'nOnesTarget','nOnesPred',
										'tpr', 'tnr','roc','acc'])
data['Epoch'] = data['Epoch'].str.split('|').str[1]

trainDf = data[data['Epoch'].str.contains('Epoch: ')]
testDf = data[data['Epoch'].str.contains('VALIDATION')]

for i in range(startEpoch, epochs+1):

	epochTrainDf = trainDf[trainDf['Epoch'].str.endswith('Epoch: '+str(i))]
	epochTestDf = testDf[testDf['Epoch'].str.endswith('Epoch: '+str(i))]
	
	# Training metrics
	trainLoss.append(calcMedian(epochTrainDf, 'Loss'))
	trainUncer.append(calcMedian(epochTrainDf, 'Uncertainty'))
	trainTpr.append(calcMedian(epochTrainDf, 'tpr'))
	trainTnr.append(calcMedian(epochTrainDf, 'tnr'))
	trainRoc.append(calcMedian(epochTrainDf, 'roc'))
	trainAcc.append(calcMedian(epochTrainDf, 'acc'))
	
	# Testing Metrics
	testLoss.append(calcMedian(epochTestDf, 'Loss'))
	testUncer.append(calcMedian(epochTestDf, 'Uncertainty'))
	testTpr.append(calcMedian(epochTestDf, 'tpr'))
	testTnr.append(calcMedian(epochTestDf, 'tnr'))
	testRoc.append(calcMedian(epochTestDf, 'roc'))
	testAcc.append(calcMedian(epochTestDf, 'acc'))
	
	i += 1

dfTrain = pd.DataFrame(list(zip(trainLoss, trainAcc, trainUncer, trainTpr, trainTnr, trainRoc,
								testLoss,  testAcc,  testUncer,  testTpr,  testTnr,  testRoc)), 
						columns =['trainLoss', 'trainAcc', 'trainUncer', 'trainTpr', 'trainTnr', 'trainRoc',
								  'testLoss',  'testAcc',  'testUncer',  'testTpr',  'testTnr',  'testRoc'])

dfTrain.to_csv(filename[:-4]+'_Processed.csv')
print('DONE')