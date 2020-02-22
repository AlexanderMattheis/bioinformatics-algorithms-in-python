from formats.data_format import DataFormat

import system.string_symbols as strings

class ThreeBlocksFormat(DataFormat):
    """
    Returns a format
    where three sequences placed one below the other.
    """
    NUM_OF_SEQUENCES = 3
    DISTANCE_BEFORE_SEQUENCE = 3

    def lines_output(self, data_objects, sequences_names):
        """
        Outputs formatted text with three sequences for the console.

        Args:
            data_objects: list containing the alignments
            sequences_names: the ids of the sequences

        Returns:
            Fasta formatted text
        """
        final_text = []

        text_sequence_1 = strings.EMPTY
        text_sequence_2 = strings.EMPTY
        text_sequence_3 = strings.EMPTY

        for alignment in data_objects:
            sequence_1 = alignment[0]
            sequence_2 = alignment[1]
            sequence_3 = alignment[2]

            alignment_length = max(len(sequence_1), len(sequence_2), len(sequence_3))

            if len(sequences_names[0]) > 0 or len(sequences_names[1]) > 0 or len(sequences_names[2]) > 0:
                num_of_chars = max(len(sequences_names[0]), len(sequences_names[1]), len(sequences_names[2])) \
                               + self.DISTANCE_BEFORE_SEQUENCE
                text_sequence_1 += sequences_names[0] + self._empty(num_of_chars - len(sequences_names[0]))
                text_sequence_2 += sequences_names[1] + self._empty(num_of_chars - len(sequences_names[1]))
                text_sequence_3 += sequences_names[2] + self._empty(num_of_chars - len(sequences_names[2]))

            for i in range(0, alignment_length):
                num_of_chars += 1

                if num_of_chars > self.LINE_LENGTH:
                    num_of_chars = 1

                    final_text.append(text_sequence_1)
                    final_text.append(text_sequence_2)
                    final_text.append(text_sequence_3)
                    final_text.append(strings.EMPTY)

                    text_sequence_1 = strings.EMPTY
                    text_sequence_2 = strings.EMPTY
                    text_sequence_3 = strings.EMPTY

                text_sequence_1 += sequence_1[i]
                text_sequence_2 += sequence_2[i]
                text_sequence_3 += sequence_3[i]

            final_text.append(text_sequence_1)
            final_text.append(text_sequence_2)
            final_text.append(text_sequence_3)
            final_text.append(strings.EMPTY)

        return final_text
