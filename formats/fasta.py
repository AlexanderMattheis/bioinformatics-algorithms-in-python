from formats.data_format import DataFormat

import re
import system.messages as messages
import system.string_symbols as strings
import sys

class Fasta(DataFormat):
    """Implements classes to check for Fasta format and create strings in Fasta format."""

    # integers
    ALIGNMENT_COUNT = 3  # tells how many lines in a fasta alignment string
    DISTANCE_BEFORE_SEQUENCE = 3

    # strings
    ONLY_ENGLISH_LETTERS = "^[A-Za-z]*$"

    def get_id(self, sequence_text):
        """
        Returns the sequence id.

        Args:
            sequence_path: path of the sequence from which you want the id

        Returns:
            the id of a Fasta sequence
        """
        return self.__get_id(sequence_text)

    def get_sequence_id(self, sequence_path):
        """
        Returns the sequence id.

        Args:
            sequence_path: path of the sequence from which you want the id

        Returns:
            the id of a Fasta sequence
        """
        sequence_text = self._get_text(sequence_path)
        sequence_id = self.__get_id(sequence_text)

        return sequence_id

    def __get_id(self, sequence_text):
        """
        Returns the ID from some sequence text.

        Args:
            sequence_text: text sequence from which you want get the id

        Returns:
            sequence from raw text
        """
        lines = sequence_text.splitlines()
        id = strings.EMPTY

        if len(lines) > 0:
            line_number = 0
            for i in range(0, len(lines)):
                if lines[i].find(strings.HEADER_START) != -1:
                    line_number = i
                    break

            line = lines[line_number]
            id = line.replace(strings.HEADER_START, strings.EMPTY).lstrip().rstrip()

        return id

    def get_sequence(self, sequence_path):
        """
        Returns text from the path of some sequence and checks if it is in Fasta format.

        Args:
            sequence_path: path of the sequence from which you want the text

        Returns:
            checked sequence from raw text
        """
        sequence_text = self._get_text(sequence_path)

        if self.__is_fasta_format(sequence_text):
            sequence = self.get_sequence_from_string(sequence_text)
        else:
            print(messages.WRONG_FASTA_FORMAT)
            sys.exit(0)

        return sequence

    def __is_fasta_format(self, single_sequence):
        """
        Checks if a text sequence has Fasta format.

        Args:
            single_sequence: text sequence you want to check

        Returns:
            True if the text sequence has Fasta format.
        """
        single_sequence = str(single_sequence)

        if len(single_sequence) == 0:
            return True

        if not single_sequence.startswith(strings.HEADER_START):
            return False

        seperated_lines = single_sequence.splitlines()
        is_fasta_format = True
        sequence_started = False

        final_text = strings.EMPTY

        for i in range(0, len(seperated_lines)):
            if i > 0 and not seperated_lines[i].startswith(strings.COMMENT) and not sequence_started:
                sequence_started = True

            # so comments in between or after a sequence are not allowed
            if sequence_started and seperated_lines[i].startswith(strings.COMMENT):
                is_fasta_format = False

            if i > 0 and not seperated_lines[i].startswith(strings.COMMENT):
                line = seperated_lines[i].replace(strings.NEW_LINE_1, strings.SPACE) \
                    .replace(strings.NEW_LINE_2, strings.SPACE) \
                    .replace(strings.NEW_LINE_3, strings.SPACE)

                final_text = final_text + line

        only_english_letters = re.match(self.ONLY_ENGLISH_LETTERS, final_text, 0)
        return is_fasta_format and only_english_letters

    def get_sequence_from_string(self, sequence_text):
        """
        Returns the sequence from some text.

        Args:
            sequence_text: text sequence you want to check

        Returns:
            sequence from raw text
        """
        lines = sequence_text.splitlines()
        final_text = strings.EMPTY

        for i in range(0, len(lines)):
            if i > 0 and lines[i].find(strings.COMMENT) == -1:
                final_text = final_text + lines[i]

        sequence = str(final_text).\
            replace(strings.NEW_LINE_1, strings.SPACE).\
            replace(strings.NEW_LINE_2, strings.SPACE).\
            replace(strings.NEW_LINE_3, strings.SPACE).\
            upper()

        return sequence

    def _get_text(self, sequence_path):
        """
        Returns text from the path of some sequence.

        Args:
            sequence_path: path of the sequence from which you want the text

        Returns:
            text of the file which was read
        """
        file = open(sequence_path, 'r')
        sequence_text = file.read()
        file.close()

        return sequence_text

    def lines_output(self, data_objects, sequences_names):
        """
        Outputs formatted text for the console.

        Args:
            data_objects: list containing the alignments
            sequences_names: the ids of the sequences

        Returns:
            Fasta formatted text
        """
        normalized_text = []

        text_sequence_1 = strings.EMPTY
        text_match_mismatch = strings.EMPTY
        text_sequence_2 = strings.EMPTY

        for alignment in data_objects:
            sequence_1 = alignment[0]
            sequence_2 = alignment[1]
            match_mismatch = self.__match_mismatch_string(alignment)
            alignment_length = max(len(sequence_1), len(sequence_2))

            if len(sequences_names[0]) > 0 or len(sequences_names[1]) > 0:
                num_of_chars = max(len(sequences_names[0]), len(sequences_names[1])) \
                          + self.DISTANCE_BEFORE_SEQUENCE
                text_sequence_1 += sequences_names[0] + self._empty(num_of_chars - len(sequences_names[0]))
                text_sequence_2 += sequences_names[1] + self._empty(num_of_chars - len(sequences_names[1]))
                text_match_mismatch += self._empty(num_of_chars)

            for i in range(0,alignment_length):
                num_of_chars += 1

                if num_of_chars > self.LINE_LENGTH:
                    num_of_chars = 1

                    normalized_text.append(text_sequence_1)
                    normalized_text.append(text_match_mismatch)
                    normalized_text.append(text_sequence_2)

                    text_sequence_1 = strings.EMPTY
                    text_match_mismatch = strings.EMPTY
                    text_sequence_2 = strings.EMPTY

                text_sequence_1 += sequence_1[i]
                text_match_mismatch += match_mismatch[i]
                text_sequence_2 += sequence_2[i]

            normalized_text.append(text_sequence_1)
            normalized_text.append(text_match_mismatch)
            normalized_text.append(text_sequence_2)
            normalized_text.append(strings.EMPTY)

            text_sequence_1 = strings.EMPTY
            text_match_mismatch = strings.EMPTY
            text_sequence_2 = strings.EMPTY

        return normalized_text

    def __match_mismatch_string(self, alignment):
        """
        Outputs a string with stars and colons for an alignment.

        Args:
            alignment: alignment of two strings

        Returns:
            string with stars and colons
        """
        match_mismatch_string = strings.EMPTY

        for i in range(0, len(alignment[0])):
            if (alignment[0][i] == strings.GAP) or (alignment[1][i] == strings.GAP):
                match_mismatch_string += strings.SPACE
            elif alignment[0][i] == alignment[1][i]:
                match_mismatch_string += strings.STAR
            else:
                match_mismatch_string += strings.COLON

        return match_mismatch_string
