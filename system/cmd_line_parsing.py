import argparse
import sys

from algorithms import available_algorithms
from algorithms.alignment.alignment_algorithm import AlignmentAlgorithm
from maths import matrix_types
from system import commands
from system import messages

SUBPARSERS_NAME = "algorithms"
STORE_TRUE = "store_true"

def start():
    """Starts parsing the command line parameters."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title=SUBPARSERS_NAME)
    __parse_gotoh(subparsers)
    __parse_needleman_wunsch(subparsers)
    __parse_needleman_wunsch_3d(subparsers)
    __parse_nussinov(subparsers)

    args = parser.parse_args()

    return __input(args)


def __parse_gotoh(subparsers):
    """"
    Parses Needleman-Wunsch parametes.

    Args:
        subparsers: set of parsers which stores the parser for Gotoh
    """
    gotoh = subparsers.add_parser(available_algorithms.GOTOH)

    matrix_type = gotoh.add_mutually_exclusive_group(required=True)
    matrix_type.add_argument(matrix_types.BLOSUM62_SHORT, matrix_types.BLOSUM62, action=STORE_TRUE)
    matrix_type.add_argument(matrix_types.PAM250_SHORT, matrix_types.PAM250, action=STORE_TRUE)

    gotoh.add_argument(commands.ALL_SHORT, commands.ALL, default=False, action=STORE_TRUE)

    gotoh.add_argument(commands.PATH_1, type=str)
    gotoh.add_argument(commands.PATH_2, type=str)
    gotoh.add_argument(commands.GAP_OPEN, type=int)
    gotoh.add_argument(commands.GAP_EXTENSION, type=int)


def __parse_needleman_wunsch(subparsers):
    """"
    Parses Needleman-Wunsch parametes.

    Args:
        subparsers: set of parsers which stores the parser for Needleman-Wunsch
    """
    needleman_wunsch = subparsers.add_parser(available_algorithms.NEEDLEMAN_WUNSCH)

    matrix_type = needleman_wunsch.add_mutually_exclusive_group(required=True)
    matrix_type.add_argument(matrix_types.BLOSUM62_SHORT, matrix_types.BLOSUM62, action=STORE_TRUE)
    matrix_type.add_argument(matrix_types.PAM250_SHORT, matrix_types.PAM250, action=STORE_TRUE)

    needleman_wunsch.add_argument(commands.ALL_SHORT, commands.ALL, default=False, action=STORE_TRUE)

    needleman_wunsch.add_argument(commands.PATH_1, type=str)
    needleman_wunsch.add_argument(commands.PATH_2, type=str)
    needleman_wunsch.add_argument(commands.GAP_OPEN, type=int)


def __parse_needleman_wunsch_3d(subparsers):
    """"
    Parses Needleman-Wunsch 3d parametes.

    Args:
        subparsers: set of parsers which stores the parser for Needleman-Wunsch
    """
    needleman_wunsch_3d = subparsers.add_parser(available_algorithms.NEEDLEMAN_WUNSCH_3D)

    matrix_type = needleman_wunsch_3d.add_mutually_exclusive_group(required=True)
    matrix_type.add_argument(matrix_types.BLOSUM62_SHORT, matrix_types.BLOSUM62, action=STORE_TRUE)
    matrix_type.add_argument(matrix_types.PAM250_SHORT, matrix_types.PAM250, action=STORE_TRUE)

    needleman_wunsch_3d.add_argument(commands.PATH, type=str)
    needleman_wunsch_3d.add_argument(commands.GAP_OPEN, type=int)


def __parse_nussinov(subparsers):
    nussinov = subparsers.add_parser(available_algorithms.NUSSINOV)
    nussinov.add_argument(commands.PATH, type=str)


def __input(args):
    """
    Creates the input which an algorithm needs.

    Args:
        args: structure containing parsed input parameters
    """
    alg_input = []

    if len(sys.argv) > 1:
        algorithm = sys.argv[1]

        if algorithm == available_algorithms.GOTOH:
            alg_input = __input_gotoh(args)
        elif algorithm == available_algorithms.NEEDLEMAN_WUNSCH:
            alg_input = __input_needleman_wunsch(args)
        elif algorithm == available_algorithms.NEEDLEMAN_WUNSCH_3D:
            alg_input = __input_needleman_wunsch_3d(args)
        elif algorithm == available_algorithms.NUSSINOV:
            alg_input = __input_nussinov(args)
        else:
            print(messages.DOES_NOT_EXIST)
    else:
        print(messages.HELP)

    return alg_input


def __input_gotoh(args):
    """
    Creates the input for Gotoh.

    Args:
        args: structure containing parsed input parameters
    """
    input = []
    input.append(available_algorithms.GOTOH)

    if args.blosum62 == True:
        input.append(matrix_types.BLOSUM62)
    else:
        input.append(matrix_types.PAM250)

    input.append(args.gap_open)
    input.append(args.gap_extension)
    input.append(args.seq_1_path)
    input.append(args.seq_2_path)

    if args.all == True:
        input.append(AlignmentAlgorithm.TRACEBACK_ALL)

    return input


def __input_needleman_wunsch(args):
    """
    Creates the input for Needleman-Wunsch.

    Args:
        args: structure containing parsed input parameters
    """
    input = []
    input.append(available_algorithms.NEEDLEMAN_WUNSCH)

    if args.blosum62 == True:
        input.append(matrix_types.BLOSUM62)
    else:
        input.append(matrix_types.PAM250)

    input.append(args.gap_open)
    input.append(args.seq_1_path)
    input.append(args.seq_2_path)

    if args.all == True:
        input.append(AlignmentAlgorithm.TRACEBACK_ALL)

    return input


def __input_needleman_wunsch_3d(args):
    """
    Creates the input for Needleman-Wunsch 3d.

    Args:
        args: structure containing parsed input parameters
    """
    input = []
    input.append(available_algorithms.NEEDLEMAN_WUNSCH_3D)

    if args.blosum62 == True:
        input.append(matrix_types.BLOSUM62)
    else:
        input.append(matrix_types.PAM250)

    input.append(args.gap_open)
    input.append(args.seq_path)

    return input

def __input_nussinov(args):
    """
    Creates the input for Needleman-Wunsch 3d.

    Args:
        args: structure containing parsed input parameters
    """
    input = []
    input.append(available_algorithms.NUSSINOV)
    input.append(args.seq_path)

    return input