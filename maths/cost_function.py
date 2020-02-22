from maths.matrix_reader import MatrixReader

import maths.matrix_types as matrix_types
import system.string_symbols as strings
import system.words as words

class CostFunction:
    """
    Returns the cost for the alignment of two amino acids.
    """

    A = "A"
    R = "R"
    N = "N"
    D = "D"
    C = "C"
    Q = "Q"
    E = "E"
    G = "G"
    H = "H"
    I = "I"
    L = "L"
    K = "K"
    M = "M"
    F = "F"
    P = "P"
    S = "S"
    T = "T"
    W = "W"
    Y = "Y"
    V = "V"

    global __matrix_type
    global __matrix

    class AminoAcids:
        """
        Returns the amino acid position in the matrix.
        """

        A = 0
        R = 1
        N = 2
        D = 3
        C = 4
        Q = 5
        E = 6
        G = 7
        H = 8
        I = 9
        L = 10
        K = 11
        M = 12
        F = 13
        P = 14
        S = 15
        T = 16
        W = 17
        Y = 18
        V = 19
        Star = 20

    class Bases:
        """
        Returns the base position in the matrix.
        """
        A = 0
        C = 1
        G = 2
        T = 3
        Star = 4

    def __init__(self, type=words.UNDEFINED):
        """
        Initializes a cost function.

        Args:
            type: the matrix you want to use as a cost function
        """
        self.__matrix_type = type

        if type != words.UNDEFINED:
            self.__matrix = MatrixReader().get_matrix(type)

    def get_value(self, a, b):
        """
        Returns the matrix value for two amino acids.

        Args:
            a: first amino acid
            b: scond amino acid
        """
        if a != strings.EMPTY and b != strings.EMPTY:
            x = self.__get_position(a)
            y = self.__get_position(b)

            return int(self.__matrix[x][y])
        elif a == strings.EMPTY and b == strings.EMPTY:
            return 0
        return strings.NEGATIVE_INFINITY

    def __get_position(self, value):
        """
        Returns the position of an amino acid in the matrix.

        Args:
            value: amino acid or base from which you want to know the position
        """
        if len(self.__matrix) > 5:
            number = self.AminoAcids()
        else:
            number = self.Bases()

        if value.upper() == self.A:
            return number.A

        elif value.upper() == self.R:
            return number.R

        elif value.upper() == self.N:
            return number.N

        elif value.upper() == self.D:
            return number.D

        elif value.upper() == self.C:
            return number.C

        elif value.upper() == self.Q:
            return number.Q

        elif value.upper() == self.E:
            return number.E

        elif value.upper() == self.G:
            return number.G

        elif value.upper() == self.H:
            return number.H

        elif value.upper() == self.I:
            return number.I

        elif value.upper() == self.L:
            return number.L

        elif value.upper() == self.K:
            return number.K

        elif value.upper() == self.M:
            return number.M

        elif value.upper() == self.F:
            return number.F

        elif value.upper() == self.P:
            return number.P

        elif value.upper() == self.S:
            return number.S

        elif value.upper() == self.T:
            return number.T

        elif value.upper() == self.W:
            return number.W

        elif value.upper() == self.Y:
            return number.Y

        elif value.upper() == self.V:
            return number.V

        else:
            return number.Star

    def get_path(self):
        """Returns matrix path if a default matrix is used."""
        if self.__matrix_type == matrix_types.PAM250:
            return matrix_types.PAM250_PATH

        elif self.__matrix_type == matrix_types.BLOSUM62:
            return matrix_types.BLOSUM62_PATH

    def set_matrix(self, matrix_path):
        """Sets the matrix which is used for evaluation."""
        self.__matrix = MatrixReader().read_in(matrix_path)
