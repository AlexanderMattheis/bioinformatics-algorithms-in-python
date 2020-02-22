import maths.matrix_types as matrix_types
import system.string_symbols as strings

class MatrixReader:
    """
    Reads in the substitution matrices.
    """
    __READ = 'r'

    def get_matrix(self, matrix_type):
        """
        Returns a substitution matrix.

        Args:
            matrix_type: the matrix_type that should be returned
        """
        if matrix_type == matrix_types.PAM250:
            return self.read_in(matrix_types.PAM250_PATH)

        elif matrix_type == matrix_types.BLOSUM62:
            return self.read_in(matrix_types.BLOSUM62_PATH)

    def read_in(self, file_path):
        """
        Reads in all lines which are not comments or empty.

        Args:
            file_path: path to the file you want read in

        Returns:
            list of matrix line values
        """
        file = open(file_path, self.__READ)
        data = []

        for line in file.readlines():
            if not line.startswith(strings.COMMENT) \
                    and not line.startswith(strings.SPACE) \
                    and len(line) > 1:  # to skip commented/empty lines

                splitted_line = line.split(strings.SPACE)

                values = []
                for value in splitted_line:  # removing all non numbers
                    if len(value) > 0 and not value.isalpha() and value.find(strings.STAR) == -1:
                        value = value.replace(strings.NEW_LINE_2, strings.EMPTY)  # remove new line symbols
                        values.append(value)

                data.append(values)

        file.close()

        return data