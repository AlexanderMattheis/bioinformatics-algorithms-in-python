import unittest

import needleman_wunsch3
from prakt.nw3 import NeedlemanWunsch3Base

class TestMethodsNeedleman3(unittest.TestCase):
    """Test class to test Needleman-Wunsch algorithm with three sequences."""

    def test_instance(self):
        """Check inheritance."""
        assert issubclass(needleman_wunsch3.NeedlemanWunsch3, NeedlemanWunsch3Base)
        assert isinstance(needleman_wunsch3.NeedlemanWunsch3(), NeedlemanWunsch3Base)

    def test_7(self):
        """Short sequence test. Test taken from Bioinformatics exercise sheet."""
        nw3 = needleman_wunsch3.NeedlemanWunsch3()
        result = nw3.run("../T_INPUT/7test_seq.fasta", "../T_INPUT/7eva.txt", -1, False)
        (id_seq1, seq1, id_seq2, seq2, id_seq3, seq3, score, alignments) = result

        assert id_seq1 == "Mus_musculus"
        assert id_seq2 == "Oryctolagus_cuniculus"
        assert id_seq3 == "Bos_taurus"
        assert seq1 == "CTCACA"
        assert seq2 == "CAC"
        assert seq3 == "GTAC"
        assert score == 2
        assert alignments == [["CTCACA", "C--AC-", "GT-AC-"]]  # order of elements is random!

    def test_empty(self):
        """Checks if everything is fine with an empty sequence."""
        nw3 = needleman_wunsch3.NeedlemanWunsch3()
        result = nw3.run("../T_INPUT/empty_test.fasta", "../T_INPUT/3eva.txt", -1, False)
        (id_seq1, seq1, id_seq2, seq2, id_seq3, seq3, score, alignments) = result

        assert id_seq1 == ""
        assert id_seq2 == ""
        assert id_seq2 == ""
        assert seq1 == ""
        assert seq2 == ""
        assert seq3 == ""
        assert score == 0
        assert alignments == [["", "", ""]]  # order of elements is random!

    def test_one_empty(self):
        """Checks if everything is fine with one empty sequence."""
        nw3 = needleman_wunsch3.NeedlemanWunsch3()
        result = nw3.run("../T_INPUT/one_empty_test.fasta", "../T_INPUT/3eva.txt", -1, True)
        (id_seq1, seq1, id_seq2, seq2, id_seq3, seq3, score, alignments) = result

        assert id_seq1 == "Mus_musculus"
        assert id_seq2 == "Oryctolagus_cuniculus"
        assert id_seq3 == "Bos_taurus"
        assert seq1 == "CTCACA"
        assert seq2 == "CAC"
        assert seq3 == ""
        assert score == -9
        assert alignments == [["CTCACA", "--CAC-", "------"]] or [["CTCACA", "C--AC-", "------"]]
