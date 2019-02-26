import pandas as pd
import numpy as np
import sys
import os


def makeMemmap(setType):
    df = pd.read_csv(setType + ".dat")#,
                     # usecols = ["dna.seq",
                     #            "normalized_log2_tags_BA9_81_April2015_LR",
                     #            "normalized_log2_tags_Vermis_62_Mar2015_LR"],
                     # dtype = {"normalized_log2_tags_BA9_81_April2015_LR": np.uint8,
                     #          "normalized_log2_tags_Vermis_62_Mar2015_LR": np.uint8})#,
                     # nrows = 1000000)
    print("1")

    rows = df.shape[0]
    bpLength = len(df["dna.seq"][0])
    numClasses = 131
    channels = 4

    inputFileName = setType + "Input.memmap"
    labelFileName = setType + "Label.memmap"

    inputMemmap = np.memmap(inputFileName, mode = "w+", shape = (1, channels, bpLength))
    labelMemmap = np.memmap(labelFileName, mode = "w+", shape = (1, numClasses))

    batchSize = 100000

    accumSamples = 0
    for i in range(0, rows, batchSize):
        endIndex = i + batchSize if i + batchSize < rows else rows
        inputs = df["dna.seq"][i:endIndex].apply(lambda x: pd.Series(list(x)))
        print("2")
        classes = df.iloc[:, 6:][i:endIndex]
        print("3")

        goodRows = ~inputs.isin(["N"]).any(1)
        print("4")

        inputs = inputs[goodRows].astype("category", categories = ["A", "C", "G", "T"], ordered = True)
        print("5")
        classes = classes[goodRows]
        print("6")

        inputNumpy = pd.get_dummies(inputs).to_numpy().reshape(-1, 1000, 4).swapaxes(1, 2)
        print("7")
        classesNumpy = classes.to_numpy(np.uint8)
        print("8")

        samples = classesNumpy.shape[0]

        # inputMemmap.resize((accumSamples + samples, channels, bpLength))
        # labelMemmap.resize((accumSamples + samples, classes))
        inputMemmap.flush()
        labelMemmap.flush()
        inputMemmap = np.memmap(inputFileName, mode = "r+", shape = (accumSamples + samples, channels, bpLength))
        labelMemmap = np.memmap(labelFileName, mode = "r+", shape = (accumSamples + samples, numClasses))

        inputMemmap[accumSamples:accumSamples + samples, :, :] = inputNumpy[:]
        labelMemmap[accumSamples:accumSamples + samples, :] = classesNumpy[:]


        accumSamples += samples



    # inputFileName = setType + "Input.memmap"
    # labelFileName = setType + "Label.memmap"

    # inputMemmap = np.memmap(inputFileName, mode = "w+", shape = inputNumpy.shape)
    # inputMemmap[:] = inputNumpy[:]
    # del inputMemmap

    # labelMemmap = np.memmap(labelFileName, mode = "w+", shape = classesNumpy.shape)
    # labelMemmap[:] = classesNumpy[:]
    # del labelMemmap

def main():
    os.chdir(sys.argv[1])	# first argument is the data directory
    rawInputFileName=sys.argv[2]	# second argument is the raw input file name
    # makeMemmap("H3K27ac_rnaSeq.Pos.Neg.tfSpecific")
    makeMemmap(rawInputFileName)

if __name__ == '__main__':
    main()
