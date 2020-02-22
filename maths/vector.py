class Vector:
    """
    Creates a vector for traceback computation.
    """
    __DEFAULT_VALUE = "S"

    global x
    global y
    global z
    global matrix_label

    def __init__(self, pos_x, pos_y, pos_z = 0):
        """
        Constructor to create a 2D or 3D vector.

        Args:
            pos_x: first coordinate
            pos_y: second coordinate
            pos_z: third coordinate (default 0)
        """

        self.x = pos_x
        self.y = pos_y
        self.z = pos_z
        self.matrix_label = self.__DEFAULT_VALUE

    def create(self, matrix_label):
        """
        Assigns a vector a matrix label.

        Args:
            matrix_label: label you want to assign
        """
        self.matrix_label = matrix_label
        return self

    def __eq__(self, other):
        """
        Defines when two vectors are equal.

        Args:
            other: second object
        """
        if other is Vector:
            return self.x == other.x and self.y == other.y and self.z == other.z \
                   and self.matrix_label == other.matrix_label
        return False
