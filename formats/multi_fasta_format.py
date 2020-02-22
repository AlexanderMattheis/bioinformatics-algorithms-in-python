from formats.fasta import Fasta

import os
import re
import system.messages as messages
import system.string_symbols as strings
import sys

class MultiFasta(Fasta):
    def get_3_sequences(self, sequences_path):
        """
        Returns three strings from the file of a given path.

        Args:
            sequences_path: the file path from which you want read in sequences
        """
        sequences_text = self._get_text(sequences_path)
        sequences = []

        if self.__is_3d_fasta_format(sequences_text):
            lines = sequences_text.splitlines()

            sequence = strings.EMPTY

            for i in range(0, len(lines)):
                if i > 0 and lines[i].startswith(strings.HEADER_START):
                    sequences.append(sequence)
                    sequence = strings.EMPTY

                sequence = sequence + lines[i].strip() + os.linesep

            sequences.append(sequence)
        else:
            print(messages.WRONG_FASTA_FORMAT)
            sys.exit(0)

        return sequences

    def __is_3d_fasta_format(self, sequences):
        """
        Checks if a text sequence with three seqeunces has Fasta format.

        Args:
            single_sequence: text sequence you want to check

        Returns:
            True if the text sequence has Fasta format.
        """
        return self.__is_multi_fasta_format(sequences, 3)

    def __is_multi_fasta_format(self, sequences, dimension=-1):
        """
       Checks if a text sequence has Fasta format.

       Args:
           single_sequence: text sequence you want to check

       Returns:
           True if the text sequence has Fasta format.
       """
        sequences = str(sequences)

        if len(sequences) == 0:
            return True

        if not sequences.startswith(strings.HEADER_START):
            return False

        seperated_lines = sequences.splitlines()
        is_fasta_format = True
        sequence_started = False
        only_english_letters = True

        english_letters_text = strings.EMPTY
        sequences_number = 0

        for i in range(0, len(seperated_lines)):
            if seperated_lines[i].startswith(strings.HEADER_START):
                sequence_started = False
                sequences_number += 1
                english_letters_text = strings.EMPTY

            if len(seperated_lines[i]) > 80 and not seperated_lines[i].startswith(strings.HEADER_START):
                is_fasta_format = False

            if i > 0 and not seperated_lines[i].startswith(strings.COMMENT) and not sequence_started:
                sequence_started = True

            # so comments in between or after a sequence are not allowed
            if sequence_started and seperated_lines[i].startswith(strings.COMMENT):
                is_fasta_format = False

            if i > 0 and not seperated_lines[i].startswith(strings.COMMENT):
                line = seperated_lines[i].replace(strings.NEW_LINE_1, strings.SPACE) \
                    .replace(strings.NEW_LINE_2, strings.SPACE) \
                    .replace(strings.NEW_LINE_3, strings.SPACE)

                if not seperated_lines[i].startswith(strings.HEADER_START):
                    english_letters_text += line

            if not re.match(self.ONLY_ENGLISH_LETTERS, english_letters_text, 0):
                only_english_letters = False

        return is_fasta_format and only_english_letters and (sequences_number == dimension or dimension == -1)

    def get_all_sequences(self, sequences_path):
        """
        Returns all strings from the file of a given path.

        Args:
            sequences_path: the file path from which you want read in sequences
        """
        sequences_text = self._get_text(sequences_path)
        sequences = []

        if self.__is_multi_fasta_format(sequences_text):
            lines = re.split(strings.NEW_LINE_ALL_TYPES, sequences_text)

            sequence = strings.EMPTY

            for i in range(0, len(lines)):
                if i > 0 and lines[i].startswith(strings.HEADER_START):
                    sequences.append(sequence)
                    sequence = strings.EMPTY

                sequence = sequence + lines[i].strip() + os.linesep

            sequences.append(sequence)
        else:
            print(messages.WRONG_FASTA_FORMAT)
            sys.exit(0)

        return sequences