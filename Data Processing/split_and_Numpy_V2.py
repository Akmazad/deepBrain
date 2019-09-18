import pandas as pd
import numpy as np
import sys
import os
import argparse
import logging

parser = argparse.ArgumentParser(description='Splitting DeepBrain data V2')
parser.add_argument('--datadir', required=True)
parser.add_argument('--datafilename', required=True)
parser.add_argument('--valid_chr_id', required=True)



def split_util(fileName, df, logger):
    # process data
    logger.info("nRow df: {}".format(df.shape[0]))
    df = df[~df["dna.seq"].str.contains('N')]  # filter out the noisy records
    logger.info("nRow filtered df: {}".format(df.shape[0]))
    # Save the 'value' numpy
    inputs = df["dna.seq"].apply(lambda x: pd.Series(list(x)))  # creating the neucleotide array from the string
    inputs = inputs.astype("category", categories=["A", "C", "G", "T"], ordered=True)  # categorize them
    logger.info("nRow values: {}".format(inputs.shape[0]))
    inputNumpy = pd.get_dummies(inputs).to_numpy().reshape(-1, inputs.shape[0], 4).swapaxes(1, 2)  # one-hot encoding of inputs
    np.save(fileName + "_value", inputNumpy[:])  # save the numpy nDarray into the file
    logger.info("value saved")
    del inputNumpy
    del inputs

    # process Label
    classes_all = df.iloc[:, 5:]
    del df  # free objects
    classes_allNumpy = classes_all.to_numpy(np.uint8); del classes_all  # labels are already binary,
    np.save(fileName + "_label_all", classes_allNumpy[:]); del classes_allNumpy  # save the numpy nDarray into the file
    logger.info("labels are all saved")


def split(filename, args, logger):
    df = pd.read_csv(filename + ".bed", low_memory=True, sep='\t')
    chrName = 'chr' + args.valid_chr_id

    train_df = df[df['chr'] != chrName]
    logger.info("training split and numpy conversion has been started")
    try:
        split_util(filename+"_trainingData_" + chrName, train_df, logger)
        logger.info("Ended")
    except Exception as e:
        print(e)
    del train_df

    logger.info("validation split and numpy conversion has been started")
    val_df = df[df['chr'] == chrName]
    split_util(filename + "_validationData_" + chrName, val_df, logger)
    logger.info("Ended")
    del val_df


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
    logger = get_logger(os.path.join(args.datadir, "split_and_numpy_V2.log"))
    os.chdir(args.datadir)  # first argument is the data directory
    rawInputFileName = args.datafilename  # second argument is the raw input file name

    split(rawInputFileName, args, logger)


if __name__ == '__main__':
    main()
