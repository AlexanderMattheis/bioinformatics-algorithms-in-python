from data.prediction_output_data import PredictionOutputData as Output
from maths.vector import Vector

import algorithms.prediction.base_pairs as pairs
import copy

class NussinovBacktracking:
    """
    Performs the traceback in the Nussinov algorithm.
    """

    global base_pairs

    def __init__(self):
        """Initializes Nussinov backtracking variables."""
        self.base_pairs = []

    def traceback_one(self, path, input_data):
        """
        Returns one traceback for a given input.
        It is based on the pseudo-code of Prof. Dr. Backofen
        http://www.bioinf.uni-freiburg.de//Lehre/Courses/2013_SS/V_RNA/slides/nussinov.pdf (slide 13)
        The code is just reversing the Nussinov recursion function.

        Args:
            path: list in which you store the traceback
            input_data: the input for which you create the traceback
        """
        current = path[-1]

        if current.x <= current.y+1:  # case 1: j<=i -> checking if already at inner diagonal or further
            return
        elif Output.table_values[current.y][current.x] == Output.table_values[current.y][current.x-1]:  # case 2
            path.append(Vector(current.x - 1, current.y))  # going one step to the left
            self.traceback_one(path, input_data)
            return
        else:  # case 3: checks for complementary bases on the left
            current_base = input_data.sequence[current.x-1]  # -1 because width of table is: len(sequence) + 1

            # for all k: i<=k<j
            for k in range(current.y, current.x-1-input_data.loop_length):
                base = input_data.sequence[k]

                # S_k and S_j complementary
                if pairs.complementary(base, current_base):
                    current_value = Output.table_values[current.y][current.x]
                    value_sum = Output.table_values[current.y][k] + Output.table_values[k+1][current.x-1] + 1

                    # if N_{i,j} = N_{i,k-1} + N_{k+1,j-1} + 1
                    if current_value == value_sum:
                        self.base_pairs.append((k+1, current.x))

                        # create vectors for new positions in the matrix from which to start
                        left = Vector(k, current.y)
                        bottom = Vector(current.x-1, k+1)

                        left_path = copy.copy(path)
                        bottom_path = copy.copy(path)

                        left_path.append(left)
                        bottom_path.append(bottom)

                        # jump to new positions
                        self.traceback_one(left_path, input_data)
                        self.traceback_one(bottom_path, input_data)
                        return