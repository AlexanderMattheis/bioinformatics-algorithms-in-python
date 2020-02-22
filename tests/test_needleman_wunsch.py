import unittest

import needleman_wunsch
from prakt.nw import NeedlemanWunschBase

class TestMethodsNeedleman(unittest.TestCase):
    """Test class to test Needleman-Wunsch algorithm."""

    def test_instance(self):
        """Check inheritance."""
        assert issubclass(needleman_wunsch.NeedlemanWunsch, NeedlemanWunschBase)
        assert isinstance(needleman_wunsch.NeedlemanWunsch(), NeedlemanWunschBase)

    def test_1(self):
        """Short sequence test."""
        nw = needleman_wunsch.NeedlemanWunsch()
        result = nw.run("../T_INPUT/1test_seq1.fasta", "../T_INPUT/1test_seq2.fasta", "../T_INPUT/1eva.txt", -2, True)
        (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

        assert id_seq1 == "Seq1"
        assert id_seq2 == "Seq2"
        assert seq1 == "AGTC"
        assert seq2 == "ATC"
        assert score == 1
        assert alignments == [["AGTC", "A-TC"]]  # order of elements is random!

    def test_2(self):
        """Long sequence test. Test taken from Bioinformatics exercise sheet."""
        nw = needleman_wunsch.NeedlemanWunsch()
        result = nw.run("../T_INPUT/2test_seq1.fasta", "../T_INPUT/2test_seq2.fasta", "../T_INPUT/2eva.txt", -2, True)
        (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

        assert id_seq1 == "Seq1"
        assert id_seq2 == "Seq2"
        assert seq1 == "CCCCGCGACUCGGGUUCAAGGG"
        assert seq2 == "GGGUGAGACCCCAGUUCAACCC"
        assert score == 33
        assert alignments == [["CCCCGCGACUCGGGUUCAAGGG", "GGGUGAGACCCCAGUUCAACCC"]]  # order of elements is random!

    def test_3(self):
        """Another sequence test. Test taken from Bioinformatics exercise sheet."""
        nw = needleman_wunsch.NeedlemanWunsch()
        result = nw.run("../T_INPUT/3test_seq1.fasta", "../T_INPUT/3test_seq2.fasta", "../T_INPUT/3eva.txt", -1, True)
        (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

        assert id_seq1 == "animal| seq1"
        assert id_seq2 == "plant| seq2"
        assert seq1 == "TCCGA"
        assert seq2 == "TACGCGC"
        assert score == 2
        assert alignments == [["T-C-CGA", "TACGCGC"]]  # order of elements is random!

    def test_empty(self):
        """Checks if everything is fine with empty sequences."""
        nw = needleman_wunsch.NeedlemanWunsch()
        result = nw.run("../T_INPUT/empty_test.fasta", "../T_INPUT/empty_test.fasta", "../T_INPUT/3eva.txt", -1, True)
        (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

        assert id_seq1 == ""
        assert id_seq2 == ""
        assert seq1 == ""
        assert seq2 == ""
        assert score == 0
        assert alignments == []  # order of elements is random!

    def test_one_empty(self):
        """Checks if everything is fine with one empty sequence."""
        nw = needleman_wunsch.NeedlemanWunsch()
        result = nw.run("../T_INPUT/empty_test.fasta", "../T_INPUT/5test_seq2.fasta", "../T_INPUT/3eva.txt", -1, True)
        (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

        assert id_seq1 == ""
        assert id_seq2 == "Seq2"
        assert seq1 == ""
        assert seq2 == "TCCGA"
        assert score == -5
        assert alignments == [["-----", "TCCGA"]]  # order of elements is random!
