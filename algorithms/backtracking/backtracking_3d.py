from algorithms.backtracking.backtracking import Backtracking
from data.alignment_output_data import AlignmentOutputData
from maths.vector import Vector

import system.string_symbols as strings

class Backtracking3D(Backtracking):
    '''
    Allows you to do a traceback in three dimensions.
    '''

    def _get_neighbours(self, pos, input_data):
        '''
        Returns all neighbours
        to which you can go from a three-dimensional matrix
        cell.

        Args:
            pos: the position from which you want get all neighbours
            input_data: the input for which you create the traceback

        Returns:
            neighbours of the cell at position pos
        '''
        neighbours = []

        start = AlignmentOutputData.table_values_3d[pos.z][pos.y][pos.x]
        case_1 = float(strings.NAN)
        case_2 = float(strings.NAN)
        case_3 = float(strings.NAN)
        case_4 = float(strings.NAN)
        case_5 = float(strings.NAN)
        case_6 = float(strings.NAN)
        case_7 = float(strings.NAN)

        cur_char_seq_1 = strings.EMPTY
        cur_char_seq_2 = strings.EMPTY
        cur_char_seq_3 = strings.EMPTY

        # two matches
        if pos.z - 1 >= 0 and pos.y - 1 >= 0 and pos.x - 1 >= 0:
            case_1 = AlignmentOutputData.table_values_3d[pos.z - 1][pos.y - 1][pos.x - 1]

        # one match
        if pos.z - 1 >= 0 and pos.y - 1 >= 0:
            case_2 = AlignmentOutputData.table_values_3d[pos.z - 1][pos.y - 1][pos.x]

        if pos.z - 1 >= 0 and pos.x - 1 >= 0:
            case_3 = AlignmentOutputData.table_values_3d[pos.z - 1][pos.y][pos.x - 1]

        if pos.y - 1 >= 0 and pos.x - 1 >= 0:
            case_4 = AlignmentOutputData.table_values_3d[pos.z][pos.y - 1][pos.x - 1]

        # no matches
        if pos.z - 1 >= 0:
            case_5 = AlignmentOutputData.table_values_3d[pos.z - 1][pos.y][pos.x]

        if pos.y - 1 >= 0:
            case_6 = AlignmentOutputData.table_values_3d[pos.z][pos.y - 1][pos.x]

        if pos.x - 1 >= 0:
            case_7 = AlignmentOutputData.table_values_3d[pos.z][pos.y][pos.x - 1]

        # get strings
        if pos.z - 1 >= 0:
            cur_char_seq_1 = input_data.sequence_a[pos.z - 1]
        if pos.y - 1 >= 0:
            cur_char_seq_2 = input_data.sequence_b[pos.y - 1]
        if pos.x - 1 >= 0:
            cur_char_seq_3 = input_data.sequence_c[pos.x - 1]

        value_yz = input_data.cost_function.get_value(cur_char_seq_1, cur_char_seq_2)
        value_xz = input_data.cost_function.get_value(cur_char_seq_1, cur_char_seq_3)
        value_xy = input_data.cost_function.get_value(cur_char_seq_2, cur_char_seq_3)

        triple_matching = start == case_1 + value_yz + value_xz + value_xy
        match_ab = start == case_2 + value_yz + 2*input_data.gap_cost
        match_ac = start == case_3 + value_xz + 2*input_data.gap_cost
        match_bc = start == case_4 + value_xy + 2*input_data.gap_cost
        mismatch_1 = start == case_5 + 2*input_data.gap_cost
        mismatch_2 = start == case_6 + 2*input_data.gap_cost
        mismatch_3 = start == case_7 + 2*input_data.gap_cost

        # two matches
        if triple_matching:
            neighbours.append(Vector(pos.x - 1, pos.y - 1, pos.z - 1))

        # one match
        if match_ab:
            neighbours.append(Vector(pos.x, pos.y - 1, pos.z - 1))

        if match_ac:
            neighbours.append(Vector(pos.x - 1, pos.y, pos.z - 1))

        if match_bc:
            neighbours.append(Vector(pos.x - 1, pos.y - 1, pos.z))

        # no matches
        if mismatch_1:
            neighbours.append(Vector(pos.x, pos.y, pos.z - 1))

        if mismatch_2:
            neighbours.append(Vector(pos.x, pos.y - 1, pos.z))

        if mismatch_3:
            neighbours.append(Vector(pos.x - 1, pos.y, pos.z))

        return neighbours