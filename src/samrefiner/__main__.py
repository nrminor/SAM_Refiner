#!/usr/bin/env python3

import itertools
import os
from multiprocessing import Pool

from . import chimera_removal as chim
from . import cli, utils
from . import sam_parsing as sam


def main():
    """
    Main functional logic
    Calls functions to get argument values, then as relevant, get reference, process sam files, perform chimera removal and collect sample's
    """
    args = cli.arg_parser()  # getting command line arguments

    if args.ref:
        ref = cli.get_ref(args)  # get the reference ID and sequence from the FASTA file
        if ref[1] == "":
            print(
                "Reference not recognized as a Fasta or Genebank format, skipping SAM parsing",
            )
        else:
            # collect SAM files to process, either from the arguments or the working directory
            SAMs = []
            try:
                args.Sam_files[0]
            except:
                for file in os.listdir(os.getcwd()):
                    if (file.lower()).endswith(".sam"):
                        SAMs.append(file)
            else:
                for files in args.Sam_files:
                    for file in files:
                        if os.path.isfile(file):
                            SAMs.append(file)
                        else:
                            print(f"Can't find {file}, skipping")

            args.ref = ""
            if ref[2] == "fasta":
                with Pool(processes=args.mp) as pool:
                    pool.starmap(
                        sam.fa_sam_parse,
                        zip(itertools.repeat(args), itertools.repeat(ref), SAMs),
                    )
            elif ref[2] == "gb":
                with Pool(processes=args.mp) as pool:
                    pool.starmap(
                        sam.gb_sam_parse,
                        zip(itertools.repeat(args), itertools.repeat(ref), SAMs),
                    )
            print("End Sam Parsing Output")
    else:
        print("No reference provided, skipping SAM parsing")

    seq_files = []
    # Begin chimera removal if enabled
    if args.chim_rm == 1 or args.deconv == 1:
        for file in os.listdir(os.getcwd()):
            if file.endswith(
                "_seqs.tsv",
            ):  # get unique sequence files for chimera removal
                seq_files.append(file[0:-16])
        with Pool(processes=args.mp) as pool:
            pool.starmap(chim.chim_process, zip(itertools.repeat(args), seq_files))

    # begin collection of sample outputs
    if args.collect == 1:
        utils.collection(args)


if __name__ == "__main__":
    main()
