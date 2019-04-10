"""
Author: Dr. A. K. M. Azad
@vafaeelabe, UNSW

Partially reproduced from Cameron Stewart's code (https://github.com/castewart/deepBrain/blob/master/scripts/csv2memmap.py)
"""

import argparse
import pandas as pd
import numpy as np
import sys
import os

parser = argparse.ArgumentParser(description='Splitting DeepBrain data')
parser.add_argument('--datadir', required=True)
parser.add_argument('--datafilename', required=True)
parser.add_argument('--valid_chr_id', required=True)


def split(rawInputFileName, valid_chr, args):
    df = pd.read_csv(filename + ".dat")
    chrName = 'chr' + valid_chr
    train_df = df[df['chr'] != chrName]
    makeMemmap(train_df, "DeepBrain_Training_wo_Chr_"+valid_chr)
    val_df = df[df['chr'] == chrName]
    makeMemmap(train_df, "DeepBrain_Validation_Chr_"+valid_chr)


def makeMemmap(df, filename, args):
    rows = df.shape[0]  # read the number of records
    bpLength = len(df["dna.seq"][0])    # read the DNA sequence length
    nAcc_classes = 2
    nRNAseq_classes = 1
    nTF_classes = 128
#     numClasses = df.shape[1]-6    # first 6 columns are ['chr', 'start', 'end', 'dna.seq', 'id', 'strand'] and remainings are classes
    channels = 4    # for A,T,G,C

    # all single-type features
    valueFileName = filename + "_Data.memmap"
    valueMemmap = np.memmap(valueFileName, mode = "w+", shape = (1, channels, bpLength))

    labelFileName_acc = filename + "_ACConly_Label.memmap"
    labelMemmap_acc = np.memmap(labelFileName_acc, mode = "w+", shape = (1, nAcc_classes))
    labelFileName_rnaSeq = filename + "_RNAseqonly_Label.memmap"
    labelMemmap_rnaSeq = np.memmap(labelFileName_rnaSeq, mode = "w+", shape = (1, nRNAseq_classes))
    labelFileName_tf = filename + "_TFonly_Label.memmap"
    labelMemmap_tf = np.memmap(labelFileName_tf, mode = "w+", shape = (1, nTF_classes))
    
    # two-types of features
    labelFileName_acc_rnaSeq = filename + "_ACC_RNAseqonly_Label.memmap"
    labelMemmap_acc_rnaSeq = np.memmap(labelFileName_acc_rnaSeq, mode = "w+", shape = (1, nAcc_classes+nRNAseq_classes))    
    labelFileName_rnaSeq_tf = filename + "_RNAseq_TFonly_Label.memmap"
    labelMemmap_rnaSeq_tf = np.memmap(labelFileName_rnaSeq_tf, mode = "w+", shape = (1, nRNAseq_classes+nTF_classes))
    labelFileName_acc_tf = filename + "_ACC_TFonly_Label.memmap"
    labelMemmap_acc_tf = np.memmap(labelFileName_acc_tf, mode = "w+", shape = (1, nAcc_classes+nTF_classes))
    
    # all features combined
    labelFileName_acc_rnaSeq_tf = filename + "_ACC_RNAseq_TF_Label.memmap"      
    labelMemmap_acc_rnaSeq_tf = np.memmap(labelFileName_acc_rnaSeq_tf, mode = "w+", shape = (1, nAcc_classes+nRNAseq_classes+nTF_classes))

    batchSize = 100000

    accumSamples = 0
    for i in range(0, rows, batchSize):
        endIndex = i + batchSize if i + batchSize < rows else rows  # until which row this iteration would work on
        values = df["dna.seq"][i:endIndex].apply(lambda x: pd.Series(list(x))) # make a matrix of neucleotides

        classes_acc = df.iloc[:, 6:8][i:endIndex]
        classes_rnaSeq = df.iloc[:, 8:9][i:endIndex]
        classes_tf = df.iloc[:, 9:][i:endIndex]
        classes_acc_rnaSeq = df.iloc[:, 6:9][i:endIndex]
        classes_rnaSeq_tf = df.iloc[:, 8:][i:endIndex]
        classes_all = df.iloc[:, 6:][i:endIndex]
        classes_acc_tf = classes_all.loc[:, classes_all.columns != "stringTie.Transcript.SpikeIns_filtered"]    # skip the RNA_seq label
        
        goodRows = ~values.isin(["N"]).any(1)

        values = values[goodRows].astype("category", categories = ["A", "C", "G", "T"], ordered = True)
        classes_acc = classes_acc[goodRows]
        classes_rnaSeq = classes_rnaSeq[goodRows]
        classes_tf = classes_tf[goodRows]
        classes_acc_rnaSeq = classes_acc_rnaSeq[goodRows]
        classes_rnaSeq_tf = classes_rnaSeq_tf[goodRows]
        classes_acc_tf = classes_acc_tf[goodRows]
        classes_all = classes_all[goodRows]

        # getting outputed in the memmap
        valueNumpy = pd.get_dummies(values).to_numpy().reshape(-1, 1000, 4).swapaxes(1, 2)
        classes_accNumpy = classes_acc.to_numpy(np.uint8)
        classes_rnaSeqNumpy = classes_rnaSeq.to_numpy(np.uint8)
        classes_tfNumpy = classes_tf.to_numpy(np.uint8)
        classes_acc_rnaSeqNumpy = classes_acc_rnaSeq.to_numpy(np.uint8)
        classes_rnaSeq_tfNumpy = classes_rnaSeq_tf.to_numpy(np.uint8)
        classes_acc_tfNumpy = classes_acc_tf.to_numpy(np.uint8)
        classes_allNumpy = classes_all.to_numpy(np.uint8)

        samples = classes_allNumpy.shape[0]     # number of rows, effectively the batch size

        valueMemmap.flush()
        labelMemmap_acc.flush()
        labelMemmap_rnaSeq.flush()
        labelMemmap_tf.flush()
        labelMemmap_acc_rnaSeq.flush()
        labelMemmap_rnaSeq_tf.flush()
        labelMemmap_acc_tf.flush()
        labelMemmap_acc_rnaSeq_tf.flush()
        
        valueMemmap = np.memmap(valueFileName, mode = "r+", shape = (accumSamples + samples, channels, bpLength))
        labelMemmap_acc = np.memmap(labelFileName_acc, mode = "r+", shape = (accumSamples + samples, nAcc_classes))
        labelMemmap_rnaSeq = np.memmap(labelFileName_rnaSeq, mode = "r+", shape = (accumSamples + samples, nRNAseq_classes))
        labelMemmap_tf = np.memmap(labelFileName_tf, mode = "r+", shape = (accumSamples + samples, nTF_classes))
        labelMemmap_acc_rnaSeq = np.memmap(labelFileName_acc_rnaSeq, mode = "r+", shape = (accumSamples + samples, nAcc_classes+nRNAseq_classes))
        labelMemmap_rnaSeq_tf = np.memmap(labelFileName_rnaSeq_tf, mode = "r+", shape = (accumSamples + samples, nRNAseq_classes+nTF_classes))
        labelMemmap_acc_tf = np.memmap(labelFileName_acc_tf, mode = "r+", shape = (accumSamples + samples, nAcc_classes+nTF_classes))
        labelMemmap_acc_rnaSeq_tf = np.memmap(labelFileName_acc_rnaSeq_tf, mode = "r+", shape = (accumSamples + samples, nAcc_classes+nRNAseq_classes+nTF_classes))

        valueMemmap[accumSamples:accumSamples + samples, :, :] = valueNumpy[:]
        labelMemmap_acc[accumSamples:accumSamples + samples, :] = classes_accNumpy[:]
        labelMemmap_rnaSeq[accumSamples:accumSamples + samples, :] = classes_rnaSeqNumpy[:]
        labelMemmap_tf[accumSamples:accumSamples + samples, :] = classes_tfNumpy[:]
        labelMemmap_acc_rnaSeq[accumSamples:accumSamples + samples, :] = classes_acc_rnaSeqNumpy[:]
        labelMemmap_rnaSeq_tf[accumSamples:accumSamples + samples, :] = classes_rnaSeq_tfNumpy[:]
        labelMemmap_acc_tf[accumSamples:accumSamples + samples, :] = classes_acc_tfNumpy[:]
        labelMemmap_acc_rnaSeq_tf[accumSamples:accumSamples + samples, :] = classes_acc_rnaSeq_tfNumpy[:]

        accumSamples += samples     # accumulative rows

        
def main():
    args = parser.parse_args()
    os.chdir(args.datadir)	# first argument is the data directory
    rawInputFileName=args.datafilename	# second argument is the raw input file name
    valid_chr = args.valid_chr_id
    split(rawInputFileName, valid_chr)

if __name__ == '__main__':
    main()
