from prakt.nussinov import NussinovBase

from algorithms.backtracking.nussinov_backtracking import NussinovBacktracking
from algorithms.prediction.structure_prediction_algorithm import StructurePredictionAlgorithm
from data.prediction_input_data import PredictionInputData
from data.prediction_output_data import PredictionOutputData
from formats.fasta import Fasta
from maths.vector import Vector

import algorithms.prediction.base_pairs as base_pairs
import os
import system.messages as messages

@NussinovBase.register
class Nussinov(NussinovBase, StructurePredictionAlgorithm):
    """
    Computes optimal, global alignment.
    """

    MINIMUM_LOOP_LENGTH = 0 # minimum distance between two sequences

    def compute(self, input):
        """
        Checks input and if it is correct,
        it runs the algorithm.

        Args:
            input: input from console
        """
        if self.__input_ok(input):
            sequence_path = input[1]
            self.run(sequence_path)
        else:
            print(messages.WRONG_PATHS)

    def __input_ok(self, input):
        """
        Checks whether the input is correct
        or not.

        Args:
            input: input from console
        """
        if len(input) >= 1:
            sequence_path = input[1]

            if not os.path.isfile(sequence_path):
                return False

            return True
        return False

    def run(self, seq_fasta_fn):
        """
        Fold RNA with Nussinov algorithm.
        Args:
            seq_fasta_fn: path to fasta file containing sequence
        Returns:
            tuple of
            (id_seq: fasta id of sequence,
             seq: sequence,
             structure: dot-bracket string of optimal folding)
        """
        self.__evaluate_parameters(seq_fasta_fn)
        self.__initialize()
        self._compute_structures()

        backtracking = NussinovBacktracking()
        PredictionOutputData.structures = self._traceback(Vector(self._tbl_len_x - 1, 0), backtracking)

        self._create_structures()
        self._output()

        return (self._data.id, self._data.sequence, PredictionOutputData.dot_bracket_structures[0])

    def __evaluate_parameters(self, seq_fasta_fn):
        """
        Stores parameters needed to run the algorithm.

        Args:
            seq_fasta_fn: path to fasta file containing the sequence
        """
        sequence = Fasta().get_sequence(seq_fasta_fn)
        sequence_id = Fasta().get_sequence_id(seq_fasta_fn)

        self._tbl_len_x = len(sequence) + 1
        self._tbl_len_y = len(sequence)
        self._data = PredictionInputData().init_nussinov(self.MINIMUM_LOOP_LENGTH, sequence, sequence_id)

    def __initialize(self):
        """
        Initializes the Nussinov matrix.
        """
        PredictionOutputData.table_values = [[0 for x in range(self._tbl_len_x)] for x in range(self._tbl_len_y)]

    def _value(self, x, y):
        """
        Recursion function of Nussinov.

        Args:
            x: current position x in the Nussinov matrix
            y: current position y in the Nussinov matrix
        """
        value_a = PredictionOutputData.table_values[y][x - 1]
        value_b = 0

        current_base = self._data.sequence[x - 1]  # -1 because the x-value of the table has size: len(sequence) + 1

        for k in range(y, x - 1 - self._data.loop_length):
            base = self._data.sequence[k]

            if base_pairs.complementary(base, current_base):
                value = PredictionOutputData.table_values[y][k] + PredictionOutputData.table_values[k + 1][x - 1] + 1
                if value > value_b:
                    value_b = value

        return max(value_a, value_b)

if __name__ == '__main__':
    # run Nussinov with some parameters
    nussinov = Nussinov()
    nussinov.run("INPUT/rna_sequence.fasta")