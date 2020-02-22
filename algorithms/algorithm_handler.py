from algorithms import available_algorithms
from gotoh import Gotoh
from needleman_wunsch import NeedlemanWunsch
from needleman_wunsch3 import NeedlemanWunsch3
from nussinov_alg import Nussinov

class AlgorithmHandler:
    """
    Selects an algorithm.
    """

    def let_compute(self, input):
        """
        Starts the algorithm with the right input parameters.

        Args:
            input: console input
        """
        algorithm = input[0]

        if algorithm == available_algorithms.GOTOH:
            procedure = Gotoh()
        elif algorithm == available_algorithms.NEEDLEMAN_WUNSCH:
            procedure = NeedlemanWunsch()
        elif algorithm == available_algorithms.NEEDLEMAN_WUNSCH_3D:
            procedure = NeedlemanWunsch3()
        elif algorithm == available_algorithms.NUSSINOV:
            procedure = Nussinov()

        procedure.compute(input)
