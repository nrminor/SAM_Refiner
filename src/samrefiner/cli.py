from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from io import TextIOWrapper
    from pathlib import Path


def arg_parser() -> argparse.Namespace:
    """
    Called to get the arguments passed by the command line and process them
    Functionality:
    uses argparse module to set up argument values based on the command line then does some conflict checking
    to make sure incompatible values aren't assigned or values aren't out of usable bounds
    Returns stored argument values
    """

    parser = argparse.ArgumentParser(
        description="process Sam files for variant information",
    )

    parser.add_argument(
        "-r",
        "-reference",
        type=argparse.FileType("r"),
        dest="ref",
        help="reference fasta or genbank file.  Only chimera removal and collections will be performed if omitted.",
    )
    parser.add_argument(
        "-S",
        "--Sam_files",
        nargs="*",
        dest="Sam_files",
        action="append",
        help='optional .sam files, can use multiple files i.e. "-S Sample1.sam -S Sample2.sam" or "-S Sample1.sam Sample2.sam"',
    )
    parser.add_argument(
        "--use_count",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) use of counts in sequence IDs, default enabled (--use_count 1) (default: 1)",
    )
    parser.add_argument(
        "--min_count",
        type=int,
        default=10,
        help="Minimum observations required to be included in sample reports; >= 1 occurance count; < 1 %% observed (.1 = 10%%), (default: .001)",
    )
    parser.add_argument(
        "--min_samp_abund",
        type=float,
        default=0.001,
        help="Minimum abundance required for inclusion in sample reports; %% observed (.1 = 10%%), (default: .001)",
    )
    parser.add_argument(
        "--min_col_abund",
        type=float,
        default=0.01,
        help="Minimum abundance required for variants to be included in collection reports; must be non-negative and  < 1, %% observed (.1 = 10%%), (default: .01)",
    )
    parser.add_argument(
        "--ntabund",
        type=float,
        default=0.001,
        help="Minimum abundance relative to total reads required for a position to be reported in the nt call output; must be non-negative and < 1, %% observed (.1 = 10%%), (default: .001)",
    )
    parser.add_argument(
        "--ntcover",
        type=int,
        default=5,
        help="Minimum coverage at a position to be reported in the nt call output. (default: 5)",
    )
    parser.add_argument(
        "--max_dist",
        type=int,
        default=40,
        help="Maximum number of variances from the reference a sequence can have to be consider in covars processing (default: 40)",
    )
    parser.add_argument(
        "--max_covar",
        type=int,
        default=8,
        help="Maximum number of variances from the reference to be reported in covars (default: 8)",
    )
    parser.add_argument(
        "--covar_tile_coverage",
        type=int,
        default=0,
        choices=[0, 1],
        help="Enable/Disable (1/0) using tiles covering positions instead of minimum nt coverage to calculate abundance of covariants (--covar_tile_coverage), require --wgs 1 (default: 0)",
    )
    parser.add_argument(
        "--AAreport",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) amino acid reporting, default enabled (--AAreport 1) (default: 1)",
    )
    parser.add_argument(
        "--AAcodonasMNP",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) reporting multiple nt changes in a single codon as one polymorphism, default enabled (--AAcodonasMNP 1), requires AAreport enabled (default: 1)",
    )
    parser.add_argument(
        "--AAcentered",
        type=int,
        default=0,
        choices=[0, 1],
        help="Enable/Disable (1/0) amino acid centered seq and covar outputs for .gb processing (--AAcentered 0), requires AAreport enabled (default: 0)",
    )
    parser.add_argument(
        "--chim_in_abund",
        type=float,
        default=0.001,
        help="Minimum abundance a unique sequence must have to be considered in chimera removal / deconvolution (default: .001)",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=1.2,
        help="Modifier for chim_rm chimera checking, default 1.2.  Higher = more sensitive, more false chimeras removed; lower = less sensitive, fewer chimeras removed",
    )
    parser.add_argument(
        "--foldab",
        type=float,
        default=1.8,
        help="Threshold for potential parent / chimera abundance ratio for chim_rm; default is 1.8",
    )
    parser.add_argument(
        "--redist",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) redistribution of chimera counts for chim_rm, default enabled (--redist 1)",
    )
    parser.add_argument(
        "--max_cycles",
        type=int,
        default=100,
        help="Max number of times chimera removal will be performed for chim_rm; default is 100",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=1,
        help="Modifier for covar pass checking, default 1.  Higher = more sensitive, more failed checks; lower = less sensitive, fewer failed checks",
    )
    parser.add_argument(
        "--autopass",
        type=float,
        default=0.3,
        help="threshold for a sequence to automatically pass the covar pass checking (default: 0.3)",
    )
    parser.add_argument(
        "--colID",
        type=str,
        default="",
        help="ID to prepend collections",
    )
    parser.add_argument(
        "--collect",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) collection step, default enabled (--collect 1) (default: 1)",
    )
    parser.add_argument(
        "--read",
        type=int,
        default=0,
        choices=[0, 1],
        help="Enable/Disable (1/0) reads output, default disabled (--read 0) (default: 0)",
    )
    parser.add_argument(
        "--nt_call",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) nt_call output, default enabled (--nt_call 1) (default: 1)",
    )
    parser.add_argument(
        "--ntvar",
        type=int,
        default=0,
        choices=[0, 1],
        help="Enable/Disable (1/0) nt_call output, default enabled (--ntvar 1) (default: 0)",
    )
    parser.add_argument(
        "--indel",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) indel output, default enabled (--indel 1) (default: 1)",
    )
    parser.add_argument(
        "--seq",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) unique seq output, default enabled (--seq 1) (default: 1)",
    )
    parser.add_argument(
        "--covar",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) covar output, default enabled (--covar 1) (default: 1)",
    )
    parser.add_argument(
        "--pass_out",
        type=int,
        default=0,
        choices=[0, 1],
        help="Enable/Disable (1/0) covar_pass output, default disabled (--pass_out 0) (default: 0)",
    )
    parser.add_argument(
        "--chim_rm",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) chim_rm output, default enabled (--chim_rm 1) (default: 1)",
    )
    parser.add_argument(
        "--deconv",
        type=int,
        default=1,
        choices=[0, 1],
        help="Enable/Disable (1/0) covar deconv, default enabled (--deconv 1) (default: 1)",
    )
    parser.add_argument(
        "--wgs",
        type=int,
        default=0,
        choices=[0, 1],
        help="Enable/Disable (1/0) covar deconv, default enabled (--wgs 1)(default: 0)",
    )
    parser.add_argument(
        "--mp",
        type=int,
        default=3,
        choices=range(1, 21),
        help="set number of processes SAM Refiner will run in parallel, default = 4 (--mp 4)",
    )

    args = parser.parse_args()

    # checking for proper range of some parameters and consistency/compatibility
    if args.wgs == 1 and (args.deconv == 1 or args.chim_rm == 1):
        args.deconv = 0
        args.chim_rm = 0
        print("WGS mode enabled, disabling chimera removal methods")

    if args.min_count < 1:
        print("--min_count must be 1 or greter, setting to 1")
        args.min_count = 1

    if args.min_samp_abund < 0 or args.min_samp_abund >= 1:
        print("--min_samp_abund must be non-negative and < 1, defaulting to .001")
        args.min_samp_abund = 0.001

    if args.min_col_abund < 0 or args.min_col_abund >= 1:
        print("--min_col_abund must be non-negative and < 1, defaulting to .01")
        args.min_col_abund = 0.01

    if args.ntabund < 0 or args.ntabund >= 1:
        print("--ntabund must be non-negative and < 1, defaulting to .001")
        args.ntabund = 0.001

    if args.ntcover < 1:
        print("--ntcover must be positive, defaulting to 5")
        args.ntcover = 5

    if args.max_dist < 0:
        print("--max_dist must be non-negative, defaulting to 40")
        args.max_dist = 40

    if args.max_covar < 0:
        print("--max_covar must be non-negative, defaulting to 8")
        args.max_covar = 8

    if args.chim_in_abund < 0 or args.chim_in_abund >= 1:
        print("--chim_in_abund must be non-negative, defaulting to 0.001")
        args.chim_in_abund = 0.001

    if args.alpha <= 0:
        print("--alpha must be positive, defaulting to 1.2")
        args.alpha = 1.2

    if args.foldab <= 0:
        print("--foldab must be positive, defaulting to 1.8")
        args.foldab = 1.8

    if args.max_cycles <= 0:
        print("--max_cycles must be positive, defaulting to 100")
        args.max_cycles = 100

    if args.beta <= 0:
        print("--beta must be positive, defaulting to 1")
        args.beta = 1

    if args.autopass <= 0:
        print("--autopass must be positive, defaulting to .3")
        args.autopass = 0.3

    return args


@dataclass
class RunSettings:
    """
    A dataclass container for all command line arguments to be used in library mode.
    """

    ref: Union[str, Path, TextIOWrapper, None] = None  # noqa: UP007
    Sam_files: Union[str, Path, list[str], None] = None  # noqa: UP007
    use_count: int = 1
    min_count: int = 10
    min_samp_abund: float = 0.01
    min_col_abund: float = 0.01
    ntabund: float = 0.001
    ntcover: int = 5
    max_dist: int = 40
    max_covar: int = 8
    covar_tile_coverage: int = 1
    AAreport: int = 1
    AAcodonasMNP: int = 1
    AAcentered: int = 0
    chim_in_abund: float = 0.001
    alpha: float = 1.2
    foldab: float = 1.8
    redist: float = 1
    max_cycles: int = 100
    beta: float = 1.0
    autopass: float = 0.3
    colID: str = ""  # noqa: N815
    collect: int = 1
    read: int = 0
    nt_call: int = 1
    ntvar: int = 0
    indel: int = 1
    seq: int = 1
    covar: int = 1
    pass_out: int = 0
    chim_rm: int = 1
    deconv: int = 1
    wgs: int = 0
    mp: int = 3


def args_to_settings(args: argparse.Namespace) -> RunSettings:
    """
    Create an instance of RunSettings using values from an argparse.Namespace.
    """
    return RunSettings(**vars(args))
