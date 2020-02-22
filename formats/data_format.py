from abc import ABCMeta, abstractmethod

import system.string_symbols as strings

class DataFormat(metaclass = ABCMeta):
    """Abstract class for implementation of different data formats."""
    LINE_LENGTH = 80

    @abstractmethod
    def lines_output(self, data_objects, sequences_names):
        """
        Outputs formatted text for the console.

        Args:
            data_objects: list containing the alignments

        Returns:
            formatted text
            sequences_names
        """
        pass

    def _empty(self, num_positions):
        """
        Create empty spaces.

        Args:
            num_positions: num of empty spaces
        """
        space = strings.EMPTY

        for i in range(0, num_positions):
            space += strings.SPACE

        return space
