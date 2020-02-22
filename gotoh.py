from prakt.gt import GotohBase

import system.string_symbols as strings
from algorithms.alignment.alignment_algorithm import AlignmentAlgorithm
from algorithms.backtracking.multi_table_backtracking import MultiTableBacktracking
from data.alignment_input_data import AlignmentInputData
from data.alignment_output_data import AlignmentOutputData
from formats.fasta import Fasta
from maths.cost_function import CostFunction
from maths.matrix_types import MATRICES
from maths.vector import Vector
from system import messages

import os

@GotohBase.register
class Gotoh(GotohBase, AlignmentAlgorithm):
    """
    Computes optimal, affine, global alignment.
    """

    def compute(self, input):
        """
        Checks input from console and if it is correct,
        it runs the algorithm.

        Args:
            input: input from console
        """
        if self.__input_ok(input):
            gap_alpha = input[2]
            gap_beta = input[3]
            sequence_a_path = input[4]
            sequence_b_path = input[5]

            if len(input) >= 7:
                traceback_mode = input[6]
            else:
                traceback_mode = strings.EMPTY

            self.run(sequence_a_path,
                     sequence_b_path,
                     CostFunction(input[1]).get_path(),
                     gap_alpha,
                     gap_beta,
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
        if len(input) >= 6:
            string_a_path = input[4]
            string_b_path = input[5]

            if not os.path.isfile(string_a_path):
                return False
            if not os.path.isfile(string_b_path):
                return False

            return True

        else:
            return False

    def run(self, seq1_fasta_fn, seq2_fasta_fn, subst_matrix_fn, affine_cost_gap_open, affine_cost_gap_extend,
            complete_traceback):
        """
        Calculate optimal alignment(s) with Needleman-Wunsch algorithm.

        Args:
            seq1_fasta_fn: path to fasta file containing first sequence
            seq2_fasta_fn: path to fasta file containing second sequence
            subst_matrix_fn: path to substitution matrix
            cost_gap_open: cost to open a gap
            cost_gap_extend: cost to extend a gap
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

        self.__evaluate_parameters(seq1_fasta_fn,
                                   seq2_fasta_fn,
                                   subst_matrix_fn,
                                   affine_cost_gap_open,
                                   affine_cost_gap_extend,
                                   complete_traceback)
        self.__initialize_tables()
        self._compute_alignments()
        # print(self._data.cost_function.get_value("O", "Z"))

        backtracking = MultiTableBacktracking()
        main_lbl = backtracking.MATRIX_LBL_MAIN

        if self._traceback_mode == self.TRACEBACK_ALL:
            AlignmentOutputData.paths = \
                self._traceback(Vector(self._tbl_len_x - 1, self._tbl_len_y - 1).create(main_lbl), True, backtracking)
            self._create_alignments()
        else:
            AlignmentOutputData.paths = \
                self._traceback(Vector(self._tbl_len_x - 1, self._tbl_len_y - 1).create(main_lbl), False, backtracking)
            self._create_alignments()

        self._output()

        return (self._data.ids[0],
                self._data.sequence_a,
                self._data.ids[1],
                self._data.sequence_b,
                AlignmentOutputData.table_values[self._tbl_len_y - 1][self._tbl_len_x - 1],
                AlignmentOutputData.alignments)

    def __evaluate_parameters(self,
                              seq1_fasta_fn,
                              seq2_fasta_fn,
                              subst_matrix_fn,
                              affine_cost_gap_open,
                              affine_cost_gap_extend,
                              complete_traceback):
        """
        Calculate optimal alignment(s) with Needleman-Wunsch algorithm.

        Args:
            seq1_fasta_fn: path to fasta file containing first sequence
            seq2_fasta_fn: path to fasta file containing second sequence
            subst_matrix_fn: path to substitution matrix
            affine_cost_gap_open: cost to open a gap
            affine_cost_gap_extend: cost to extend a gap
            complete_traceback: If True, return all optimal alignments. Otherwise choose a random alignment.
        """
        cost_function = CostFunction()
        cost_function.set_matrix(subst_matrix_fn)

        cost_gap_open = affine_cost_gap_open
        cost_gap_extend = affine_cost_gap_extend
        sequence_a = Fasta().get_sequence(seq1_fasta_fn)
        sequence_b = Fasta().get_sequence(seq2_fasta_fn)
        sequence_id_a = Fasta().get_sequence_id(seq1_fasta_fn)
        sequence_id_b = Fasta().get_sequence_id(seq2_fasta_fn)

        if complete_traceback:
            self._traceback_mode = self.TRACEBACK_ALL

        self.__set_parameters(cost_function, cost_gap_open, cost_gap_extend, sequence_a, sequence_b, sequence_id_a,
                              sequence_id_b)

    def __set_parameters(self, cost_function, affine_cost_gap_open, affine_cost_gap_extend, sequence_a, sequence_b,
                         sequence_id_a,
                         sequence_id_b):
        """
        Sets the parameters needed for the calculation with the  Needleman-Wunsch algorithm.

        Args:
               cost_function: initialized evaluation function used to evaluate alignment
               affine_cost_gap_open: cost to open a gap
               affine_cost_gap_extend: cost to extend a gap
               sequence_a: first sequence
               sequence_b: second sequence
               sequence_id_a: first sequence id
               sequence_id_b: second sequence id
        """
        self._tbl_len_x = len(sequence_b) + 1
        self._tbl_len_y = len(sequence_a) + 1
        self._data = AlignmentInputData().init_gotoh(cost_function, affine_cost_gap_open, affine_cost_gap_extend,
                                                     sequence_a,
                                                     sequence_b,
                                                     ids=[sequence_id_a, sequence_id_b])

    def __initialize_tables(self):
        """
        Initializes the Gotoh matrices.
        """
        self.__init_similarity_table()
        self.__init_horizontal_gap_cost_table()
        self.__init_vertical_gap_cost_table()

    def __init_similarity_table(self):
        """
        Initializes the Gotoh main matrix.
        """
        AlignmentOutputData.table_values = [[0 for x in range(self._tbl_len_x)] for x in range(self._tbl_len_y)]

        for k in range(1, self._tbl_len_y):
            for y in range(0, k):
                AlignmentOutputData.table_values[k][0] = self._data.gap_alpha + self._data.gap_beta * k

        for k in range(1, self._tbl_len_x):
            for x in range(0, k):
                AlignmentOutputData.table_values[0][k] = self._data.gap_alpha + self._data.gap_beta * k

    def __init_horizontal_gap_cost_table(self):
        """
        Initializes the Gotoh matrix for horizontal costs.
        """
        AlignmentOutputData.table_horizontal_gaps = [[0 for x in range(self._tbl_len_x)] for x in
                                                     range(self._tbl_len_y)]

        for x in range(1, self._tbl_len_x):
            AlignmentOutputData.table_horizontal_gaps[0][x] = strings.GAP

        for y in range(1, self._tbl_len_y):
            AlignmentOutputData.table_horizontal_gaps[y][0] = strings.NEGATIVE_INFINITY

    def __init_vertical_gap_cost_table(self):
        """
        Initializes the Gotoh matrix for vertical costs.
        """
        AlignmentOutputData.table_vertical_gaps = [[0 for x in range(self._tbl_len_x)] for x in range(self._tbl_len_y)]

        for x in range(1, self._tbl_len_x):
            AlignmentOutputData.table_vertical_gaps[0][x] = strings.NEGATIVE_INFINITY

        for y in range(1, self._tbl_len_y):
            AlignmentOutputData.table_vertical_gaps[y][0] = strings.GAP

    def _value(self, char_a, char_b, x, y):
        """
        Maximum scoring function which returns the score of a specific position.

        Args:
            char_a: char from string a at position x, y
            char_b: char from string b at position x, y
            x: current position x in the Gotoh matrix
            y: current position y in the Gotoh matrix
        """
        AlignmentOutputData.table_horizontal_gaps[y][x] = max(
            AlignmentOutputData.table_horizontal_gaps[y][x - 1] + self._data.gap_beta,
            AlignmentOutputData.table_values[y][x - 1] + self._data.gap_opening)

        AlignmentOutputData.table_vertical_gaps[y][x] = max(
            AlignmentOutputData.table_vertical_gaps[y - 1][x] + self._data.gap_beta,
            AlignmentOutputData.table_values[y - 1][x] + self._data.gap_opening)

        value = max(AlignmentOutputData.table_horizontal_gaps[y][x],
                    AlignmentOutputData.table_values[y - 1][x - 1] + self._data.cost_function.get_value(char_a, char_b),
                    AlignmentOutputData.table_vertical_gaps[y][x])

        return value

if __name__ == '__main__':
    # run Needleman-Wunsch with some parameters
    gt = Gotoh()
    gt.run("INPUT/sequence1.fasta", "INPUT/sequence2.fasta", "INPUT/pam250.txt", -8, -1, True)
    '''
    Optimal Alignment Score: 9
    Optimal Alignments:
    Seq1   TACGCAGA
           *   *:**
    Seq2   T---CCGA

    Seq1   TACGCAGA
           *:*   **
    Seq2   TCC---GA

    Seq1   TACGCAGA
           *:**   *
    Seq2   TCCG---A
    '''
