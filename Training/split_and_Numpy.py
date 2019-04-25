import pandas as pd
import numpy as np
import sys
import os
import argparse
import logging

parser = argparse.ArgumentParser(description='Splitting DeepBrain data')
parser.add_argument('--datadir', required=True)
parser.add_argument('--datafilename', required=True)
parser.add_argument('--valid_chr_id', required=True)


def split(filename, args, logger):
    df = pd.read_csv(filename + ".dat", low_memory=True)
    chrName = 'chr' + args.valid_chr_id

    train_df = df[df['chr'] != chrName]
    # train_df.to_csv(filename+"_trainingData_" + chrName + ".dat", sep=",", index=False)
    logger.info("training split and numpy conversion has been started")
    makeNumpyArrAlCombination(filename+"_trainingData_" + chrName, train_df, logger)
    logger.info("Ended")
    del train_df

    logger.info("validation split and numpy conversion has been started")
    val_df = df[df['chr'] == chrName]
    # val_df.to_csv(filename+"_validationData_" + chrName + ".dat", sep=",", index=False)
    makeNumpyArrAlCombination(filename + "_validationData_" + chrName, val_df, logger)
    logger.info("Ended")
    del val_df


def makeNumpyArrAlCombination(fileName, df, logger):
    # df = pd.read_csv(fileName + ".dat")  # read the file as a dataframe
    logger.info("nRow df: {}".format(df.shape[0]))
    df = df[~df["dna.seq"].str.contains('N')]   # filter out the noisy records
    logger.info("nRow filtered df: {}".format(df.shape[0]))


    # Save the 'value' numpy
    inputs = df["dna.seq"].apply(lambda x: pd.Series(list(x)))  # creating the neucleotide array from the string
    inputs = inputs.astype("category", categories=["A", "C", "G", "T"], ordered=True)  # categorize them
    logger.info("nRow values: {}".format(inputs.shape[0]))
    inputNumpy = pd.get_dummies(inputs).to_numpy().reshape(-1, 1000, 4).swapaxes(1, 2)  # one-hot encoding of inputs
    np.save(fileName + "_value", inputNumpy[:])  # save the numpy nDarray into the file
    logger.info("value saved")
    del inputNumpy
    del inputs


    # Save all combinations of labels (i.e. ACC only, RNA_seq only, TFs only, or their combinations)
    classes_acc = df.iloc[:, 6:8]
    classes_rnaSeq = df.iloc[:, 8:9]
    classes_tf = df.iloc[:, 9:]
    classes_acc_rnaSeq = df.iloc[:, 6:9]
    classes_rnaSeq_tf = df.iloc[:, 8:]
    classes_all = df.iloc[:, 6:]
    classes_acc_tf = classes_all.loc[:, classes_all.columns != "stringTie.Transcript.SpikeIns_filtered"]  # skip the RNA_seq label

    del df  # free objects

    classes_accNumpy = classes_acc.to_numpy(np.uint8); del classes_acc   # labels are already binary,
    classes_rnaSeqNumpy = classes_rnaSeq.to_numpy(np.uint8); del classes_rnaSeq  # labels are already binary,
    classes_tfNumpy = classes_tf.to_numpy(np.uint8); del classes_tf  # labels are already binary,
    classes_acc_rnaSeqNumpy = classes_acc_rnaSeq.to_numpy(np.uint8); del classes_acc_rnaSeq  # labels are already binary,
    classes_rnaSeq_tfNumpy = classes_rnaSeq_tf.to_numpy(np.uint8); del classes_rnaSeq_tf  # labels are already binary,
    classes_acc_tfNumpy = classes_acc_tf.to_numpy(np.uint8); del classes_acc_tf  # labels are already binary,
    classes_allNumpy = classes_all.to_numpy(np.uint8); del classes_all  # labels are already binary,

    np.save(fileName + "_label_ACC", classes_accNumpy[:]); del classes_accNumpy # save the numpy nDarray into the file
    np.save(fileName + "_label_rnaSeq", classes_rnaSeqNumpy[:]); del classes_rnaSeqNumpy    # save the numpy nDarray into the file
    np.save(fileName + "_label_tf", classes_tfNumpy[:]); del classes_tfNumpy  # save the numpy nDarray into the file
    np.save(fileName + "_label_acc_rnaSeq", classes_acc_rnaSeqNumpy[:]); del classes_acc_rnaSeqNumpy  # save the numpy nDarray into the file
    np.save(fileName + "_label_rnaSeq_tf", classes_rnaSeq_tfNumpy[:]); del classes_rnaSeq_tfNumpy  # save the numpy nDarray into the file
    np.save(fileName + "_label_acc_tf", classes_acc_tfNumpy[:]); del classes_acc_tfNumpy  # save the numpy nDarray into the file
    np.save(fileName + "_label_all", classes_allNumpy[:]); del classes_allNumpy  # save the numpy nDarray into the file

    logger.info("labels are all saved")


def get_logger(file_path):
    """ Make python logger """
    # [!] Since tensorboardX use default logger (e.g. logging.info()), we should use custom logger
    logger = logging.getLogger('db2')
    log_format = '%(asctime)s | %(message)s'
    formatter = logging.Formatter(log_format, datefmt='%m/%d %I:%M:%S %p')
    file_handler = logging.FileHandler(file_path)
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(logging.INFO)

    return logger

def main():
    args = parser.parse_args()
    logger = get_logger(os.path.join(args.datadir, "split_and_numpy.log"))
    os.chdir(args.datadir)  # first argument is the data directory
    rawInputFileName = args.datafilename  # second argument is the raw input file name

    split(rawInputFileName, args, logger)


if __name__ == '__main__':
    main()
