from prakt.nw import NeedlemanWunschBase

from algorithms.alignment.alignment_algorithm import AlignmentAlgorithm
from algorithms.backtracking.backtracking import Backtracking
from data.alignment_input_data import AlignmentInputData
from data.alignment_output_data import AlignmentOutputData
from formats.fasta import Fasta
from maths.cost_function import CostFunction
from maths.vector import Vector
from system import messages

import copy
import os
import system.string_symbols as strings

@NeedlemanWunschBase.register
class NeedlemanWunsch(NeedlemanWunschBase, AlignmentAlgorithm):
    """
    Computes optimal, global alignment.
    """

    def compute(self, input):
        """
        Checks input and if it is correct,
        it runs the algorithm.

        Args:
            input: input from console
        """
        if self.__input_ok(input):
            gap_cost = input[2]
            sequence_a_path = input[3]
            sequence_b_path = input[4]

            if len(input) >= 6:
                traceback_mode = input[5]
            else:
                traceback_mode = strings.EMPTY

            self.run(sequence_a_path,
                     sequence_b_path,
                     CostFunction(input[1]).get_path(), gap_cost,
                     traceback_mode == self.TRACEBACK_ALL)
        else:
            print(messages.WRONG_PATHS)

    def __input_ok(self, input):
        """
        Checks whether the input is correct
        or not.

        Args:
            input: input from console
        """
        if len(input) >= 5:
            string_a_path = input[3]
            string_b_path = input[4]

            if not os.path.isfile(string_a_path):
                return False
            if not os.path.isfile(string_b_path):
                return False

            return True

        return False

    def run(self, seq1_fasta_fn, seq2_fasta_fn, subst_matrix_fn, cost_gap_open, complete_traceback):
        """
        Calculate optimal alignment(s) with Needleman-Wunsch algorithm and returns outputs of this algorithm
        for testing purposes.

        Args:
            seq1_fasta_fn: path to fasta file containing first sequence
            seq2_fasta_fn: path to fasta file containing second sequence
            subst_matrix_fn: path to substitution matrix
            cost_gap_open: cost to open a gap
            complete_traceback: If True, return all optimal alignments. Otherwise choose a random alignment.

        Returns:
            tuple of
            (id_seq1: fasta id of first sequence,
            seq1: first sequence,
            id_seq2: fasta id of second sequence,
            seq2: second sequence,
            score: score of optimal alignment,
            [(aln_string_seq1, aln_string_seq2), ...]: list of tuples containing optimal alignments)
        """
        self.__evaluate_parameters(seq1_fasta_fn, seq2_fasta_fn, subst_matrix_fn, cost_gap_open, complete_traceback)
        self.__initialize_global()
        self._compute_alignments()

        backtracking = Backtracking()
        if self._traceback_mode == self.TRACEBACK_ALL:
            AlignmentOutputData.paths = self._traceback(Vector(self._tbl_len_x - 1, self._tbl_len_y - 1), True,
                                                        backtracking)
            self._create_alignments()
        else:
            AlignmentOutputData.paths = self._traceback(Vector(self._tbl_len_x - 1, self._tbl_len_y - 1), False,
                                                        backtracking)
            self._create_alignments()

        self._output()

        return (self._data.ids[0],
                self._data.sequence_a,
                self._data.ids[1],
                self._data.sequence_b,
                AlignmentOutputData.table_values[self._tbl_len_y - 1][self._tbl_len_x - 1],
                AlignmentOutputData.alignments)

    def __evaluate_parameters(self, seq1_fasta_fn, seq2_fasta_fn, subst_matrix_fn, cost_gap_open, complete_traceback):
        """
        Stores parameters needed to run the algorithm.

        Args:
            seq1_fasta_fn: path to fasta file containing first sequence
            seq2_fasta_fn: path to fasta file containing second sequence
            subst_matrix_fn: path to substitution matrix
            cost_gap_open: cost to open a gap
            complete_traceback: If True, return all optimal alignments. Otherwise choose a random alignment.
        """
        cost_function = CostFunction()
        cost_function.set_matrix(subst_matrix_fn)

        gap_cost = cost_gap_open
        sequence_a = Fasta().get_sequence(seq1_fasta_fn)
        sequence_b = Fasta().get_sequence(seq2_fasta_fn)
        sequence_id_a = Fasta().get_sequence_id(seq1_fasta_fn)
        sequence_id_b = Fasta().get_sequence_id(seq2_fasta_fn)

        if complete_traceback:
            self._traceback_mode = self.TRACEBACK_ALL

        self.__set_parameters(cost_function, gap_cost, sequence_a, sequence_b, sequence_id_a, sequence_id_b)

    def __set_parameters(self, cost_function, gap_cost, sequence_a, sequence_b, sequence_id_a, sequence_id_b):
        """
        Sets the parameters needed for the calculation with the  Needleman-Wunsch algorithm.

        Args:
               cost_function: initialized evaluation function used to evaluate alignment
               gap_cost: costs for a gap
               sequence_a: first sequence
               sequence_b: second sequence
               sequence_id_a: first sequence id
               sequence_id_b: second sequence id
        """
        self._tbl_len_x = len(sequence_b) + 1
        self._tbl_len_y = len(sequence_a) + 1
        self._data = AlignmentInputData().init_needleman(cost_function,
                                                         gap_cost,
                                                         sequence_a,
                                                         sequence_b,
                                                         ids=[sequence_id_a, sequence_id_b])

    def __initialize_global(self):
        """
        Initializes the Needleman-Wunsch matrix.
        """
        AlignmentOutputData.table_values = [[0 for x in range(self._tbl_len_x)] for x in range(self._tbl_len_y)]

        for k in range(1, self._tbl_len_y):
            for y in range(0, k):
                AlignmentOutputData.table_values[k][0] = AlignmentOutputData.table_values[y][0] + self._data.gap_cost

        for k in range(1, self._tbl_len_x):
            for x in range(0, k):
                AlignmentOutputData.table_values[0][k] = AlignmentOutputData.table_values[0][x] + self._data.gap_cost

    def _value(self, char_a, char_b, x, y):
        """
        Maximum scoring function which returns the score of a specific position.

        Args:
            char_a: char from string a at position x, y
            char_b: char from string b at position x, y
            x: current position x in the Needleman-Wunsch matrix
            y: current position y in the Needleman-Wunsch matrix
        """

        value = max(AlignmentOutputData.table_values[y][x - 1] + self._data.gap_cost,
                    AlignmentOutputData.table_values[y - 1][x - 1] + self._data.cost_function.get_value(char_a, char_b),
                    AlignmentOutputData.table_values[y - 1][x] + self._data.gap_cost)

        return value

    def get_new_table(self, cost_function, gap_cost, sequence_a, sequence_b):
        """
        Returns a Needleman-Wunsch 2D table needed for the 3D version of Needleman-Wunsch.

        Args:
            cost_function: initialized evaluation function used to evaluate alignment
            gap_cost: costs for a gap
            sequence_a: first sequence
            sequence_b: second sequence

        Returns:
            a deepcopy of a a 2D Needleman-Wunsch matrix
        """

        self.__set_parameters(cost_function, gap_cost, sequence_a, sequence_b, strings.EMPTY, strings.EMPTY)
        self.__initialize_global()
        self._compute_alignments()
        return copy.deepcopy(AlignmentOutputData.table_values)

    def get_data(self, cost_function, gap_cost, sequence_a, sequence_b):
        """
        Returns relevant data needed for example for the Feng-Doolittle algorithm.

            Args:
                cost_function: initialized evaluation function used to evaluate alignment
                gap_cost: costs for a gap
                sequence_a: first sequence
                sequence_b: second sequence

        Returns:
            a tuple with the alignment score, the number of gaps and the length of the alignment
        """

        self.__set_parameters(cost_function, gap_cost, sequence_a, sequence_b)
        self.__initialize_global()
        self._compute_alignments()
        backtracking = Backtracking()
        paths = self._traceback(Vector(self._tbl_len_x - 1, self._tbl_len_y - 1), False, backtracking)
        gaps = self._count_gaps(paths[0])
        return (AlignmentOutputData.table_values[self._tbl_len_y - 1][self._tbl_len_x - 1], gaps, len(paths[0]) - 1)

if __name__ == '__main__':
    # run Needleman-Wunsch with some parameters
    nw = NeedlemanWunsch()
    nw.run("INPUT/sequence1.fasta", "INPUT/sequence2.fasta", "INPUT/pam250.txt", -1, True)
    '''
    Optimal Alignment Score: 31
    Optimal Alignments:
    Seq1   TACGCAGA
           * * * **
    Seq2   T-C-C-GA
    '''