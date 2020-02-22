from data.alignment_output_data import AlignmentOutputData
from formats.fasta import Fasta
from system import messages

import system.string_symbols as strings
import system.words as words

class AlignmentAlgorithm:
    """
    Declares all functions which are shared by alignment algorithms.
    """

    TRACEBACK_ALL = "all"

    global _data
    global _tbl_len_x
    global _tbl_len_y
    global _traceback_mode

    def __init__(self):
        """
        Initializes parameters used by alignment algorithms.
        """
        self._traceback_mode = strings.EMPTY

    def compute(self, input_data):
        """
        Abstract function which has to be overriden to avoid an error message.
        It checks input from console and if it is correct,
        it should run the algorithm.

        Args:
            input: input from console
        """
        raise NotImplementedError(messages.OVERWRITE_METHOD)

    def _compute_alignments(self):
        """
        Computes the alignment matrix by using the scoring-function.
        """
        for y in range(1, self._tbl_len_y):
            cur_char_seq1 = self._data.sequence_a[y - 1]

            for x in range(1, self._tbl_len_x):
                cur_char_seq2 = self._data.sequence_b[x - 1]
                AlignmentOutputData.table_values[y][x] = self._value(cur_char_seq1, cur_char_seq2, x, y)

    def _create_alignments(self):
        """Creates the alignments by going through the traceback paths."""
        alignments = []

        for i in range(0, len(AlignmentOutputData.paths)):
            alignments.append(self._create_alignment(AlignmentOutputData.paths[i]))

        AlignmentOutputData.alignments = alignments

    def _create_alignment(self, path):
        """
        Creates one alignment by going through a traceback path.

        Args:
            path: from which you want get the traceback
        """
        path.reverse()
        start = path[0]

        alignment_a = strings.EMPTY
        alignment_b = strings.EMPTY

        cur_char_align_1 = 0
        cur_char_align_2 = 0

        for i in range(0, len(path) - 1):
            if path[i + 1].x - path[i].x > 0 and path[i + 1].y - path[i].y > 0:
                alignment_a += self._data.sequence_a[start.y + cur_char_align_1]
                alignment_b += self._data.sequence_b[start.x + cur_char_align_2]
                cur_char_align_1 += 1
                cur_char_align_2 += 1
            elif path[i + 1].x - path[i].x > 0:
                alignment_a += strings.GAP
                alignment_b += self._data.sequence_b[start.x + cur_char_align_2]
                cur_char_align_2 += 1
            elif path[i + 1].y - path[i].y > 0:
                alignment_a += self._data.sequence_a[start.y + cur_char_align_1]
                alignment_b += strings.GAP
                cur_char_align_1 += 1

        return [alignment_a, alignment_b]

    def _traceback(self, vec, all, backtracking):
        """
        Starts the traceback.

        Args:
            vec: from which you begin with the alignment
            all: if this is True then all tracebacks are returned
            backtracking: a Backtracking instance
        """
        path = []
        path.append(vec)

        if all:
            backtracking.traceback_all(path, self._data)
        else:
            backtracking.traceback_one(path, self._data)

        return backtracking.paths

    def _output(self):
        """
        Outputs a nice string on the console.
        """
        print(words.OPTIMAL_SCORE + str(AlignmentOutputData.table_values[self._tbl_len_y - 1][self._tbl_len_x - 1]))

        print(words.OPTIMAL_ALIGNMENTS)
        formatted_lines = Fasta().lines_output(AlignmentOutputData.alignments, self._data.ids)

        for alignment in formatted_lines:
            print(alignment)

    def _count_gaps(self, path):
        """
        Counts the number of gaps given a traceback path.

        Args:
            path: traceback path from which you calculate the number of gaps in the final alignment
        """
        path.reverse()
        gaps = 0

        for i in range(0, len(path) - 1):
            if path[i + 1].x - path[i].x > 0 and path[i + 1].y - path[i].y > 0:
                pass
            elif path[i + 1].x - path[i].x > 0:
                gaps += 1
            elif path[i + 1].y - path[i].y > 0:
                gaps += 1

        return gaps