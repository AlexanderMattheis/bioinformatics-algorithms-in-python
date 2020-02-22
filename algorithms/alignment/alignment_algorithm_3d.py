from algorithms.alignment.alignment_algorithm import AlignmentAlgorithm
from data.alignment_output_data import AlignmentOutputData
from formats.three_blocks_format import ThreeBlocksFormat

import os
import system.messages as messages
import system.string_symbols as strings
import system.words as words

class AlignmentAlgorithm3D(AlignmentAlgorithm):
    """
    Declares all functions which are shared by three-dimensional alignment algorithms.
    """

    __COMMA = ", "
    _OUTPUT_MAX_LINE_LENGTH = 77

    global _cost_function_name
    global _data
    global _tbl_len_z

    def _compute_alignments(self):
        """
        Computes the three-dimensional alignment matrix by using the scoring-function.

        Args:
            input: input from console
        """

        for z in range(1, self._tbl_len_z):
            cur_char_seq1 = self._data.sequence_a[z - 1]

            for y in range(1, self._tbl_len_y):
                cur_char_seq2 = self._data.sequence_b[y - 1]

                for x in range(1, self._tbl_len_x):
                    cur_char_seq3 = self._data.sequence_c[x - 1]

                    value_yz = self._data.cost_function.get_value(cur_char_seq1, cur_char_seq2)
                    value_xz = self._data.cost_function.get_value(cur_char_seq1, cur_char_seq3)
                    value_xy = self._data.cost_function.get_value(cur_char_seq2, cur_char_seq3)

                    case_1 = AlignmentOutputData.table_values_3d[z - 1][y - 1][x - 1] + value_yz + value_xz + value_xy

                    case_2 = AlignmentOutputData.table_values_3d[z - 1][y - 1][x] + value_yz + 2 * self._data.gap_cost
                    case_3 = AlignmentOutputData.table_values_3d[z - 1][y][x - 1] + value_xz + 2 * self._data.gap_cost
                    case_4 = AlignmentOutputData.table_values_3d[z][y - 1][x - 1] + value_xy + 2 * self._data.gap_cost

                    case_5 = AlignmentOutputData.table_values_3d[z - 1][y][x] + 2 * self._data.gap_cost
                    case_6 = AlignmentOutputData.table_values_3d[z][y - 1][x] + 2 * self._data.gap_cost
                    case_7 = AlignmentOutputData.table_values_3d[z][y][x - 1] + 2 * self._data.gap_cost

                    AlignmentOutputData.table_values_3d[z][y][x] \
                        = max(case_1, case_2, case_3, case_4, case_5, case_6, case_7)

    def _create_alignment(self, path):
        """Creates the alignments by going through the traceback paths."""
        path.reverse()
        start = path[0]

        alignment_a = strings.EMPTY
        alignment_b = strings.EMPTY
        alignment_c = strings.EMPTY

        cur_char_align_1 = 0
        cur_char_align_2 = 0
        cur_char_align_3 = 0

        for i in range(0, len(path) - 1):
            if path[i + 1].x - path[i].x > 0 and path[i + 1].y - path[i].y > 0 and path[i + 1].z - path[i].z > 0:
                alignment_a += self._data.sequence_a[start.z + cur_char_align_1]
                alignment_b += self._data.sequence_b[start.y + cur_char_align_2]
                alignment_c += self._data.sequence_c[start.x + cur_char_align_3]
                cur_char_align_1 += 1
                cur_char_align_2 += 1
                cur_char_align_3 += 1

            elif path[i + 1].z - path[i].z > 0 and path[i + 1].y - path[i].y > 0:
                alignment_a += self._data.sequence_a[start.z + cur_char_align_1]
                alignment_b += self._data.sequence_b[start.y + cur_char_align_2]
                alignment_c += strings.GAP
                cur_char_align_1 += 1
                cur_char_align_2 += 1
            elif path[i + 1].z - path[i].z > 0 and path[i + 1].x - path[i].x > 0:
                alignment_a += self._data.sequence_a[start.z + cur_char_align_1]
                alignment_b += strings.GAP
                alignment_c += self._data.sequence_c[start.x + cur_char_align_3]
                cur_char_align_1 += 1
                cur_char_align_3 += 1
            elif path[i + 1].y - path[i].y > 0 and path[i + 1].x - path[i].x > 0:
                alignment_a += strings.GAP
                alignment_b += self._data.sequence_b[start.y + cur_char_align_2]
                alignment_c += self._data.sequence_c[start.x + cur_char_align_3]
                cur_char_align_2 += 1
                cur_char_align_3 += 1

            elif path[i + 1].z - path[i].z > 0:
                alignment_a += self._data.sequence_a[start.z + cur_char_align_1]
                alignment_b += strings.GAP
                alignment_c += strings.GAP
                cur_char_align_1 += 1
            elif path[i + 1].y - path[i].y > 0:
                alignment_a += strings.GAP
                alignment_b += self._data.sequence_b[start.y + cur_char_align_2]
                alignment_c += strings.GAP
                cur_char_align_2 += 1
            elif path[i + 1].x - path[i].x > 0:
                alignment_a += strings.GAP
                alignment_b += strings.GAP
                alignment_c += self._data.sequence_c[start.x + cur_char_align_3]
                cur_char_align_3 += 1
        return [alignment_a, alignment_b, alignment_c]

    def _input(self):
        """Test function - only for testing purposes during the development."""
        sequence_a = self._data.sequence_a
        sequence_b = self._data.sequence_b
        sequence_c = self._data.sequence_c
        max_len = self._OUTPUT_MAX_LINE_LENGTH

        print(words.INPUT)
        print(words.GAP_COST + str(self._data.gap_cost))
        print(words.SEQUENCE_1 + sequence_a[:max_len] + (strings.ABBREVATION if len(sequence_a) > max_len else strings.EMPTY))
        print(words.SEQUENCE_2 + sequence_b[:max_len] + (strings.ABBREVATION if len(sequence_b) > max_len else strings.EMPTY))
        print(words.SEQUENCE_3 + sequence_c[:max_len] + (strings.ABBREVATION if len(sequence_c) > max_len else strings.EMPTY))

    def _output(self):
        """
        Outputs a nice string on the console.
        """
        print(os.linesep + words.OUTPUT)
        print(words.OPTIMAL_SCORE + str(AlignmentOutputData.table_values_3d[self._tbl_len_z - 1][self._tbl_len_y - 1][self._tbl_len_x - 1]))

        print(words.OPTIMAL_ALIGNMENTS)
        formatted_lines = ThreeBlocksFormat().lines_output(AlignmentOutputData.alignments, self._data.ids)

        for line in formatted_lines:
            print(line)

        print(words.OTHER_PARAMETERS, end=strings.SPACE)
        print(str(AlignmentOutputData.table_values_xy[self._tbl_len_x - 1][self._tbl_len_y - 1]) + self.__COMMA +
              str(AlignmentOutputData.table_values_xz[self._tbl_len_x - 1][self._tbl_len_z - 1]) + self.__COMMA +
              str(AlignmentOutputData.table_values_yz[self._tbl_len_y - 1][self._tbl_len_z - 1]))

        print(os.linesep + messages.HINT)