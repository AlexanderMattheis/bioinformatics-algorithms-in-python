from data.prediction_output_data import PredictionOutputData
from formats.character_format import CharacterFormat

import system.messages as messages
import system.string_symbols as symbols
import system.words as words

class StructurePredictionAlgorithm:
    """
    Declares all functions which are shared by alignment algorithms.
    """

    global _data
    global _tbl_len_x
    global _tbl_len_y

    def compute(self, input_data):
        """
        Abstract function which has to be overriden to avoid an error message.
        It checks input from console and if it is correct,
        it should run the algorithm.

        Args:
            input: input from console
        """
        raise NotImplementedError(messages.OVERWRITE_METHOD)

    def _compute_structures(self):
        """
        Computes the prediction matrix by using the recursion function.
        """
        for i in range(0, self._tbl_len_x - 1):
            #print("---")
            for x, y in zip(range(2, self._tbl_len_x), range(0, self._tbl_len_y-1 - i)):
                PredictionOutputData.table_values[y][x + i] = self._value(x + i, y)
                #print((y, x + i))

    def _create_structures(self):
        """Creates the structures by going through the traceback paths."""
        structures_data = PredictionOutputData.structures[0]
        PredictionOutputData.dot_bracket_structures = []

        self._create_structure(structures_data)

    def _create_structure(self, structure):
        """
        Creates one structure by going through a traceback path.

        Args:
            structure: from which you want get the traceback
        """
        dot_structure = []

        for i in range(0, self._tbl_len_y):
            dot_structure.append(symbols.DOT)
        for points in structure:
            dot_structure[points[0]-1] = symbols.BRACKET_LEFT
            dot_structure[points[1]-1] = symbols.BRACKET_RIGHT

        string_dot_structure = symbols.EMPTY.join(dot_structure)
        PredictionOutputData.dot_bracket_structures.append(string_dot_structure)

    def _traceback(self, vec, backtracking):
        """
        Starts the traceback.

        Args:
            vec: from which you begin with the traceback
            backtracking: a Backtracking instance
        """
        path = []
        path.append(vec)

        backtracking.traceback_one(path, self._data)
        if len(PredictionOutputData.table_values) > 0:
            return [backtracking.base_pairs, PredictionOutputData.table_values[0][self._tbl_len_x - 1]]
        else:
            return [backtracking.base_pairs, 0]

    def _output(self):
        """
        Outputs a nice string on the console.
        """
        if len(PredictionOutputData.table_values) > 0:
            print(words.OPTIMAL_BASEPAIRS_NUMBER + str(PredictionOutputData.table_values[0][self._tbl_len_x - 1]))
        else:
            print(words.OPTIMAL_BASEPAIRS_NUMBER + str(0))

        print(words.OPTIMAL_STRUCTURE)
        formatted_lines = CharacterFormat().lines_output(PredictionOutputData.dot_bracket_structures,
                                                         self._data.sequence)

        for structure in formatted_lines:
            print(structure)



