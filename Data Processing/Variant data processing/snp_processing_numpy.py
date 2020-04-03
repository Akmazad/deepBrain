import pandas as pd
import numpy as np
import sys
import os
import argparse
import logging

parser = argparse.ArgumentParser(description='Splitting DeepBrain data V2')
parser.add_argument('--datadir', required=True)
parser.add_argument('--datafilename', required=True)

def util(fileName, df, logger):
    # process data
    logger.info("nRow inside df: {}".format(df.shape[0]))
    # df = df[~df.iloc[:,0].str.contains('N')]  # filter out the noisy records
    # logger.info("nRow filtered df: {}".format(df.shape[0]))
    # Save the 'value' numpy
    try:
        inputs = df.apply(lambda x: pd.Series(list(x)))  # creating the neucleotide array from the string

        logger.info(inputs.shape)
        inputs = inputs.astype("category", categories=["A", "C", "G", "T"], ordered=True)  # categorize them
        logger.info("nRow values: {}".format(inputs.shape[0]))
        inputNumpy = pd.get_dummies(inputs).to_numpy().reshape(-1, inputs.shape[1], 4).swapaxes(1,2)  # one-hot encoding of inputs
        np.save(fileName + "_value", inputNumpy[:])  # save the numpy nDarray into the file
        logger.info("value saved")
        del inputNumpy
        del inputs
    except Exception as e:
        print(e)

def process(filename, args, logger):
    df = pd.read_csv(filename + ".csv", low_memory=True, sep=',')

    # ref seq processing
    logger.info("nRow df: {}".format(df.shape[0]))
    df = df[~df["refDNAseq"].str.contains('N')]  # filter out the noisy records
    df = df[~df["varDNAseq"].str.contains('N')]  # filter out the noisy records

    logger.info("nRow filtered df: {}".format(df.shape[0]))
    # Save the 'refSeq' numpy
    inputs_ref = df["refDNAseq"].apply(lambda x: pd.Series(list(x)))  # creating the neucleotide array from the string
    inputs_ref = inputs_ref.astype("category", categories=["A", "C", "G", "T"], ordered=True)  # categorize them
    logger.info("nRow values: {}".format(inputs_ref.shape[0]))
    inputNumpy = pd.get_dummies(inputs_ref).to_numpy().reshape(-1, inputs_ref.shape[1], 4).swapaxes(1, 2)  # one-hot encoding of inputs
    np.save(filename + "_refSeq_data", inputNumpy[:])  # save the numpy nDarray into the file
    logger.info("value saved")
    del inputNumpy
    del inputs_ref

    # Save the 'varSeq' numpy
    inputs_var = df["varDNAseq"].apply(lambda x: pd.Series(list(x)))  # creating the neucleotide array from the string
    inputs_var = inputs_var.astype("category", categories=["A", "C", "G", "T"], ordered=True)  # categorize them
    logger.info("nRow values: {}".format(inputs_var.shape[0]))
    inputNumpy = pd.get_dummies(inputs_var).to_numpy().reshape(-1, inputs_var.shape[1], 4).swapaxes(1,2)  # one-hot encoding of inputs
    np.save(filename + "_varSeq_data", inputNumpy[:])  # save the numpy nDarray into the file
    logger.info("value saved")
    del inputNumpy
    del inputs_var

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
    logger = get_logger(os.path.join(args.datadir, "Seq_Numpy.log"))
    os.chdir(args.datadir)  # first argument is the data directory
    rawInputFileName = args.datafilename  # second argument is the raw input file name

    process(rawInputFileName, args, logger)


if __name__ == '__main__':
    main()
