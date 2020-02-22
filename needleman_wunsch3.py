from prakt.nw3 import NeedlemanWunsch3Base

from algorithms.alignment.alignment_algorithm_3d import AlignmentAlgorithm3D
from algorithms.backtracking.backtracking_3d import Backtracking3D
from data.alignment_input_data import AlignmentInputData
from data.alignment_output_data import AlignmentOutputData
from formats.fasta import Fasta
from formats.multi_fasta_format import MultiFasta
from maths.cost_function import CostFunction
from maths.vector import Vector
from needleman_wunsch import NeedlemanWunsch

import os
import system.messages as messages
import system.string_symbols as strings

@NeedlemanWunsch3Base.register
class NeedlemanWunsch3(NeedlemanWunsch3Base, AlignmentAlgorithm3D):
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
            sequences_path = input[3]

            self.run(sequences_path,
                     CostFunction(input[1]).get_path(),
                     gap_cost,
                     False)
        else:
            print(messages.WRONG_PATHS)

    def __input_ok(self, input):
        """
        Checks whether the input is correct
        or not.

        Args:
            input: input from console
        """
        if len(input) >= 4:
            string_sequences_path = input[3]

            if not os.path.isfile(string_sequences_path):
                return False

            return True
        else:
            return False

    def run(self, seq_fasta_fn, subst_matrix_fn, cost_gap_open, complete_traceback):
        """
        Calculate optimal alignment with Needleman-Wunsch 3D algorithm and returns outputs of this algorithm
        for testing purposes.

        Args:
            seq_fasta_fn: path to fasta file containing three sequences
            subst_matrix_fn: path to substitution matrix
            cost_gap_open: cost to open a gap
            complete_traceback: If True, return all optimal alignments. Otherwise choose a random alignment.

        Returns:
            tuple of
            (id_seq1: fasta id of first sequence,
            seq1: first sequence,
            id_seq2: fasta id of second sequence,
            seq2: second sequence,
            id_seq3: fasta id of third sequence,
            seq3: third sequence,
            score: score of optimal alignment,
            [(aln_string_seq1, aln_string_seq2), ...]: list of tuples containing optimal alignments)
        """

        self.__evaluate_parameters(seq_fasta_fn, subst_matrix_fn, cost_gap_open)
        self._input()
        self.__initialize_global()
        self._compute_alignments()

        backtracking = Backtracking3D()
        AlignmentOutputData.paths = \
            self._traceback(Vector(self._tbl_len_x - 1, self._tbl_len_y - 1, self._tbl_len_z - 1), False, backtracking)
        self._create_alignments()
        self._output()

        return (self._data.ids[0],
                self._data.sequence_a,
                self._data.ids[1],
                self._data.sequence_b,
                self._data.ids[2],
                self._data.sequence_c,
                AlignmentOutputData.table_values_3d[self._tbl_len_z - 1][self._tbl_len_y - 1][self._tbl_len_x - 1],
                AlignmentOutputData.alignments)

    def __evaluate_parameters(self, seq_fasta_fn, subst_matrix_fn, cost_gap_open):
        """
        Stores parameters needed to run the algorithm.

        Args:
            seq_fasta_fn: path to fasta file containing the sequences
            subst_matrix_fn: path to substitution matrix
            cost_gap_open: cost to open a gap
            complete_traceback: If True, return all optimal alignments. Otherwise choose a random alignment.
        """
        cost_function = CostFunction()
        cost_function.set_matrix(subst_matrix_fn)
        gap_cost = cost_gap_open

        sequences = MultiFasta().get_3_sequences(seq_fasta_fn)
        if len(sequences) == 3:
            sequence_a = Fasta().get_sequence_from_string(sequences[0])
            sequence_b = Fasta().get_sequence_from_string(sequences[1])
            sequence_c = Fasta().get_sequence_from_string(sequences[2])
            sequence_id_a = Fasta().get_id(sequences[0])
            sequence_id_b = Fasta().get_id(sequences[1])
            sequence_id_c = Fasta().get_id(sequences[2])
        else:
            sequence_a = strings.EMPTY
            sequence_b = strings.EMPTY
            sequence_c = strings.EMPTY
            sequence_id_a = strings.EMPTY
            sequence_id_b = strings.EMPTY
            sequence_id_c = strings.EMPTY

        self._tbl_len_x = len(sequence_c) + 1
        self._tbl_len_y = len(sequence_b) + 1
        self._tbl_len_z = len(sequence_a) + 1
        self._data = AlignmentInputData().\
            init_needleman(cost_function, gap_cost, sequence_a, sequence_b, sequence_c,
                           ids=[sequence_id_a, sequence_id_b, sequence_id_c])

    def __initialize_global(self):
        """
        Initializes the three-dimensional Needleman-Wunsch matrix.
        """
        AlignmentOutputData.table_values_3d = \
            [[[0 for x in range(self._tbl_len_x)] for x in range(self._tbl_len_y)] for x in range(self._tbl_len_z)]

        self.__init_2d_needleman_tables()

        for x in range(0, self._tbl_len_x):
            for y in range(0, self._tbl_len_y):
                AlignmentOutputData.table_values_3d[0][y][x] \
                    = AlignmentOutputData.table_values_xy[x][y] + (x + y) * self._data.gap_cost

        for x in range(0, self._tbl_len_x):
            for z in range(0, self._tbl_len_z):
                AlignmentOutputData.table_values_3d[z][0][x] \
                    = AlignmentOutputData.table_values_xz[x][z] + (x + z) * self._data.gap_cost

        for y in range(0, self._tbl_len_y):
            for z in range(0, self._tbl_len_z):
                AlignmentOutputData.table_values_3d[z][y][0] \
                    = AlignmentOutputData.table_values_yz[y][z] + (y + z) * self._data.gap_cost

    def __init_2d_needleman_tables(self):
        """
        Computes two-dimensional Needleman-Wunsch to create the faces of the three-dimensional Needleman-Wunsch matrix.
        """
        AlignmentOutputData.table_values_xy = NeedlemanWunsch(). \
            get_new_table(self._data.cost_function, self._data.gap_cost, self._data.sequence_c, self._data.sequence_b)

        AlignmentOutputData.table_values_xz = NeedlemanWunsch(). \
            get_new_table(self._data.cost_function, self._data.gap_cost, self._data.sequence_c, self._data.sequence_a)

        AlignmentOutputData.table_values_yz = NeedlemanWunsch(). \
            get_new_table(self._data.cost_function, self._data.gap_cost, self._data.sequence_b, self._data.sequence_a)

        # print(OutputData.table_values_xy)
        # print(OutputData.table_values_xz)
        # print(OutputData.table_values_yz)

if __name__ == '__main__':
    nw3 = NeedlemanWunsch3()
    nw3.run("INPUT/sequences.fasta", "INPUT/pam250.txt", -8, False)