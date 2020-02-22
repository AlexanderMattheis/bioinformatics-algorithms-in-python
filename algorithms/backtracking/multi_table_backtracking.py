from algorithms.backtracking.backtracking import Backtracking
from data.alignment_output_data import AlignmentOutputData
from maths.vector import Vector

import system.string_symbols as strings

class MultiTableBacktracking(Backtracking):

    MATRIX_LBL_MAIN = "S"
    MATRIX_LBL_HORIZONTAL_GAPS = "-Q"
    MATRIX_LBL_VERTICAL_GAPS = "-P"

    def __init__(self):
        """Initializes multi table backtracking variables."""
        self.paths = []

    def _get_neighbours(self, pos, input_data):
        """
        Returns all neighbours
        to which you can go from a matrix
        cell.

        Args:
            pos: the position from which you want get all neighbours
            input_data: the input for which you create the traceback

        Returns:
            neighbours of the cell at position pos
        """
        neighbours = []

        if pos.matrix_label == self.MATRIX_LBL_HORIZONTAL_GAPS:
            return self.__get_neighbours_in_Q(pos, input_data)
        elif pos.matrix_label == self.MATRIX_LBL_VERTICAL_GAPS:
            return self.__get_neighbours_in_P(pos, input_data)

        start = AlignmentOutputData.table_values[pos.y][pos.x]
        diagonal = float(strings.NAN)
        matrix_p = float(strings.NAN)
        matrix_q = float(strings.NAN)
        up = float(strings.NAN)
        left = float(strings.NAN)

        cur_char_seq_1 = strings.EMPTY
        cur_char_seq_2 = strings.EMPTY

        if pos.y - 1 >= 0 and pos.x - 1 >= 0:
            diagonal = AlignmentOutputData.table_values[pos.y - 1][pos.x - 1]

        if pos.y > 0:
            matrix_p = AlignmentOutputData.table_vertical_gaps[pos.y][pos.x]

        if pos.x > 0:
            matrix_q = AlignmentOutputData.table_horizontal_gaps[pos.y][pos.x]

        # marginal case
        if pos.y - 1 >= 0 and pos.x == 0:
            up = AlignmentOutputData.table_values[pos.y - 1][pos.x]

        # marginal case
        if pos.x - 1 >= 0 and pos.y == 0:
            left = AlignmentOutputData.table_values[pos.y][pos.x - 1]

        if pos.y - 1 >= 0:
            cur_char_seq_1 = input_data.sequence_a[pos.y - 1]
        if pos.x - 1 >= 0:
            cur_char_seq_2 = input_data.sequence_b[pos.x - 1]

        if cur_char_seq_1 != strings.EMPTY and cur_char_seq_2 != strings.EMPTY:
            matching = start == diagonal + input_data.cost_function.get_value(cur_char_seq_1, cur_char_seq_2)
        else:
            matching = False

        #matching = start == diagonal + (0 if cur_char_seq_1 == cur_char_seq_2 else -1)
        change_to_p = start == matrix_p
        change_to_q = start == matrix_q

        # marginal cases
        deletion = start == up + input_data.gap_beta
        insertion = start == left + input_data.gap_beta

        if matching:
            neighbours.append(Vector(pos.x - 1, pos.y - 1).create(self.MATRIX_LBL_MAIN))

        if change_to_p:
            neighbours.append(Vector(pos.x, pos.y).create(self.MATRIX_LBL_VERTICAL_GAPS))

        if change_to_q:
            neighbours.append(Vector(pos.x, pos.y).create(self.MATRIX_LBL_HORIZONTAL_GAPS))

        # marginal case
        if insertion:
            neighbours.append(Vector(pos.x - 1, pos.y).create(self.MATRIX_LBL_MAIN))

        # marginal case
        if deletion:
            neighbours.append(Vector(pos.x, pos.y - 1).create(self.MATRIX_LBL_MAIN))

        # marginal case
        if not (matching or change_to_p or change_to_q or insertion or deletion) and (pos.x != 0 or pos.y != 0):
            neighbours.append(Vector(0, 0).create(self.MATRIX_LBL_MAIN))

        return neighbours

    def __get_neighbours_in_P(self, pos, input_data):
        """
        Returns all neighbours
        to which you can go from a matrix
        cell.

        Args:
            pos: the position from which you want get all neighbours
            input_data: the input for which you create the traceback

        Returns:
            neighbours of the cell at position pos
        """
        neighbours = []

        start = AlignmentOutputData.table_vertical_gaps[pos.y][pos.x]

        if pos.y >= 0:
            matrix_p = AlignmentOutputData.table_vertical_gaps[pos.y - 1][pos.x]
            matrix_d = AlignmentOutputData.table_values[pos.y - 1][pos.x]

        up_in_p = start == matrix_p + input_data.gap_beta
        up_in_d = start == matrix_d + input_data.gap_opening

        if up_in_p:
            neighbours.append(Vector(pos.x, pos.y - 1).create(self.MATRIX_LBL_VERTICAL_GAPS))

        if up_in_d:
            neighbours.append(Vector(pos.x, pos.y - 1).create(self.MATRIX_LBL_MAIN))

        return neighbours

    def __get_neighbours_in_Q(self, pos, input_data):
        """
        Returns all neighbours
        to which you can go from a matrix
        cell.

        Args:
            pos: the position from which you want get all neighbours
            input_data: the input for which you create the traceback

        Returns:
            neighbours of the cell at position pos
        """
        neighbours = []

        start = AlignmentOutputData.table_horizontal_gaps[pos.y][pos.x]

        if pos.y >= 0:
            matrix_q = AlignmentOutputData.table_horizontal_gaps[pos.y][pos.x - 1]
            matrix_d = AlignmentOutputData.table_values[pos.y][pos.x - 1]

        left_in_q = start == matrix_q + input_data.gap_beta
        left_in_d = start == matrix_d + input_data.gap_opening

        if left_in_q:
            neighbours.append(Vector(pos.x - 1, pos.y).create(self.MATRIX_LBL_HORIZONTAL_GAPS))

        if left_in_d:
            neighbours.append(Vector(pos.x - 1, pos.y).create(self.MATRIX_LBL_MAIN))

        return neighbours
