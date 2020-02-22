from formats.data_format import DataFormat

import system.string_symbols as strings

class CharacterFormat(DataFormat):
    """
    Returns the dot bracket notation.
    """

    def lines_output(self, data_objects, sequence):
        """
        Outputs formatted text in dot-bracket notation for the console.

        Args:
            data_objects: list containing the alignments
            sequence: the sequence for which a dot-bracket notation is created

        Returns:
            Fasta formatted text
        """
        final_text = []
        num_of_chars = 0

        for structure in data_objects:
            structure_length = len(structure)
            dot_bracket_notation = strings.EMPTY
            sequence_part = strings.EMPTY

            for i in range(0, structure_length):
                num_of_chars += 1

                if num_of_chars > self.LINE_LENGTH:
                    num_of_chars = 1

                    final_text.append(sequence_part)
                    final_text.append(dot_bracket_notation)
                    final_text.append(strings.EMPTY)

                    sequence_part = strings.EMPTY
                    dot_bracket_notation = strings.EMPTY

                sequence_part += sequence[i]
                dot_bracket_notation += structure[i]

            final_text.append(sequence_part)
            final_text.append(dot_bracket_notation)
            final_text.append(strings.EMPTY)

        return final_text