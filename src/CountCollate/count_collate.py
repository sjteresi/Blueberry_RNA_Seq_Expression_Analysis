#!/usr/bin/env python3

"""
Collate the expression counts together
"""

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np
import os
import logging
import coloredlogs
import argparse
from functools import reduce


def load_and_collate_data(input_dir):
    """Import the count data"""
    only_files = []
    for my_file in os.listdir(input_dir):
        if os.path.isfile(os.path.join(input_dir, my_file)):
            only_files.append((os.path.join(input_dir, my_file), my_file))

    all_df = []
    for my_file_pair in only_files:
        df = pd.read_csv(my_file_pair[0], sep="\t", names=["Gene", my_file_pair[1]])
        df.Gene = df["Gene"].str.split("-mRNA-1", n=1, expand=True)[0]
        all_df.append(df)

    df_merged = reduce(
        lambda left, right: pd.merge(left, right, on=["Gene"], how="outer"), all_df
    )
    return df_merged


def save_data(output_dir, merged_data):
    merged_data.to_csv(
        os.path.join(output_dir, "AllCounts_Blueberry.tsv"), sep="\t", index=False
    )


def validate_args(args, logger):
    """Raise error if input argument is invalid"""
    if not os.path.isdir(args.input_dir):
        logger.critical("argument 'input_dir' is not a directory")
        raise ValueError("%s is not a directory" % (abs_path))
    if not os.path.isdir(args.output_dir):
        logger.critical("argument 'output_dir' is not a directory")
        raise ValueError("%s is not a directory" % (abs_path))


if __name__ == "__main__":
    """Command line interface to collate the counts"""
    parser = argparse.ArgumentParser(description="Collate Counts")
    path_main = os.path.abspath(__file__)
    parser.add_argument("input_dir", type=str, help="Parent directory of input counts")
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=os.path.join(
            path_main, "../../../../Projects/Blueberry_Data/Counts/Collate/"
        ),
        help="Parent directory of output data",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debgging level to DEBUG"
    )

    args = parser.parse_args()
    args.input_dir = os.path.abspath(args.input_dir)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    logger.info("Start processing directory '%s'" % (args.input_dir))
    logger.info("Importing input counts")
    all_merged = load_and_collate_data(args.input_dir)
    logger.info("Saving data")
    save_data(args.output_dir, all_merged)
    logger.info("Done!")
