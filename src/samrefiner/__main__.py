#!/usr/bin/env python3
from __future__ import annotations

import itertools
import os
from multiprocessing import Pool
from typing import TYPE_CHECKING, Union

from samrefiner import chimera_removal as chimera
from samrefiner import cli, utils
from samrefiner import sam_parsing as sam

if TYPE_CHECKING:
    import argparse


def run(settings: Union[argparse.Namespace, utils.RunSettings]) -> None:  # noqa: UP007
    """
    Handle the flow of data through samrefiner's functions, handling cases where a
    reference may not have been provided.

    To control the run and supply inputs, either use the samrefiner command line
    interface or instantiate a RunSettings first.
    """

    if settings.ref:
        ref = utils.get_ref(
            settings,
        )  # get the reference ID and sequence from the FASTA file
        if ref[1] == "":
            print(
                "Reference not recognized as a Fasta or Genbank format, skipping SAM parsing",
            )
        else:
            # collect SAM files to process, either from the arguments or the working directory
            SAMs = []
            try:
                settings.Sam_files[0]
            except:
                for file in os.listdir(os.getcwd()):
                    if (file.lower()).endswith(".sam"):
                        SAMs.append(file)
            else:
                for files in settings.Sam_files:
                    for file in files:
                        if os.path.isfile(file):
                            SAMs.append(file)
                        else:
                            print(f"Can't find {file}, skipping")

            settings.ref = ""
            if ref[2] == "fasta":
                with Pool(processes=settings.mp) as pool:
                    pool.starmap(
                        sam.fa_sam_parse,
                        zip(itertools.repeat(settings), itertools.repeat(ref), SAMs),
                    )
            elif ref[2] == "gb":
                with Pool(processes=settings.mp) as pool:
                    pool.starmap(
                        sam.gb_sam_parse,
                        zip(itertools.repeat(settings), itertools.repeat(ref), SAMs),
                    )
            print("End Sam Parsing Output")
    else:
        print("No reference provided, skipping SAM parsing")

    seq_files = []
    # Begin chimera removal if enabled
    if settings.chim_rm == 1 or settings.deconv == 1:
        for file in os.listdir(os.getcwd()):
            if file.endswith(
                "_seqs.tsv",
            ):  # get unique sequence files for chimera removal
                seq_files.append(file[0:-16])
        with Pool(processes=settings.mp) as pool:
            pool.starmap(
                chimera.chim_process,
                zip(itertools.repeat(settings), seq_files),
            )

    # begin collection of sample outputs
    if settings.collect == 1:
        utils.collection(settings)


def main() -> None:
    """
    Main functional logic
    Calls functions to get argument values, then as relevant, get reference, process sam
    files, perform chimera removal and collect sample's output files.
    """
    # get command line argements if running as an executable and parse them into a
    # settings instance
    args = cli.arg_parser()
    settings = cli.args_to_settings(args)

    # run the workflow
    run(settings)


if __name__ == "__main__":
    main()
