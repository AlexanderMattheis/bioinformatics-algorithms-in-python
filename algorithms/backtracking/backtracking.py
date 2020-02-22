from data.alignment_output_data import AlignmentOutputData
from maths.vector import Vector

import copy
import random
import system.string_symbols as strings

class Backtracking:
    """
    Allows you to do a traceback.
    """
    global paths

    def __init__(self):
        """Initializes backtracking variables."""
        self.paths = []

    def traceback_all(self, path, input_data):
        """
        Returns all tracebacks for given input
        by returning all allowed paths from the bottom right corner
        to the left upper corner of the matrix.

        Args:
            path: list in which you store the traceback
            input_data: the input for which you create the traceback
        """
        current = path[-1]
        neighbours = self._get_neighbours(current, input_data)

        for vec in neighbours:
            if (vec.x == 0 and vec.y == 0 and vec.z == 0):
                path.append(vec)
                self.paths.append(copy.deepcopy(path))
                path.pop()
            else:
                path.append(vec)
                self.traceback_all(path, input_data)
                path.pop()

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

        start = AlignmentOutputData.table_values[pos.y][pos.x]
        diagonal = float(strings.NAN)
        up = float(strings.NAN)
        left = float(strings.NAN)

        cur_char_seq_1 = strings.EMPTY
        cur_char_seq_2 = strings.EMPTY

        if pos.y - 1 >= 0 and pos.x - 1 >= 0:
            diagonal = AlignmentOutputData.table_values[pos.y - 1][pos.x - 1]

        if pos.y - 1 >= 0:
            up = AlignmentOutputData.table_values[pos.y - 1][pos.x]

        if pos.x - 1 >= 0:
            left = AlignmentOutputData.table_values[pos.y][pos.x - 1]

        if pos.y - 1 >= 0:
            cur_char_seq_1 = input_data.sequence_a[pos.y - 1]
        if pos.x - 1 >= 0:
            cur_char_seq_2 = input_data.sequence_b[pos.x - 1]

        matching = start == diagonal + input_data.cost_function.get_value(cur_char_seq_1, cur_char_seq_2)
        deletion = start == up + input_data.gap_cost
        insertion = start == left + input_data.gap_cost

        if matching:
            neighbours.append(Vector(pos.x - 1, pos.y - 1))

        if insertion:
            neighbours.append(Vector(pos.x - 1, pos.y))

        if deletion:
            neighbours.append(Vector(pos.x, pos.y - 1))

        return neighbours

    def traceback_one(self, path, input_data):
        """
        Returns randomly one traceback for given input
        by randomly choosing a successor for given node.

        Args:
            path: list in which you store the traceback
            input_data: the input for which you create the traceback
        """
        current = path[-1]
        neighbours = self._get_neighbours(current, input_data)

        if len(neighbours) != 0:
            random_number = int(random.random() * len(neighbours))
            successor = neighbours[random_number]

            path.append(successor)
            self.traceback_one(path, input_data)
        elif len(neighbours) == 0:
            self.paths.append(path)