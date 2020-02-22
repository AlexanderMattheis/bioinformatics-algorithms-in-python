import unittest

import gotoh
from prakt.gt import GotohBase

class TestMethodsGotoh(unittest.TestCase):
    """Test class to test Gotoh algorithm."""

    def test_instance(self):
        """Check inheritance."""
        assert issubclass(gotoh.Gotoh, GotohBase)
        assert isinstance(gotoh.Gotoh(), GotohBase)

    def test_4(self):
        """Test if run function can be called."""
        gt = gotoh.Gotoh()
        result = gt.run("../T_INPUT/4test_seq1.fasta", "../T_INPUT/4test_seq2.fasta", "../T_INPUT/4eva.txt", -3, -1,
                        True)
        (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

        assert id_seq1 == "Seq1"
        assert id_seq2 == "Seq2"
        assert seq1 == "TGGA"
        assert seq2 == "GG"
        assert score == -6
        assert alignments == [['TGGA', '--GG'], ['TGGA', 'GG--']]  # order of elements can be random!

    def test_5(self):
        """Test if run function can be called."""
        gt = gotoh.Gotoh()
        result = gt.run("../T_INPUT/5test_seq1.fasta", "../T_INPUT/5test_seq2.fasta", "../T_INPUT/5eva.txt", -4, -1,
                        True)
        (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

        assert id_seq1 == "Seq1"
        assert id_seq2 == "Seq2"
        assert seq1 == "TACGCAGA"
        assert seq2 == "TCCGA"
        assert score == -3
        assert alignments == [['TACGCAGA', 'T---CCGA'], ['TACGCAGA', 'TCC---GA'],
                              ['TACGCAGA', 'TCCG---A']]  # order of elements can be random!

    def test_6(self):
        """Test if run function can be called."""
        gt = gotoh.Gotoh()
        result = gt.run("../T_INPUT/6test_seq1.fasta", "../T_INPUT/6test_seq2.fasta", "../T_INPUT/6eva.txt", -4, -1,
                        True)
        (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

        assert id_seq1 == "Seq1"
        assert id_seq2 == "Seq2"
        assert seq1 == "CC"
        assert seq2 == "ACCT"
        assert score == -7
        assert alignments == [['--CC', 'ACCT'], ['CC--', 'ACCT']]  # order of elements can be random!

    def test_empty(self):
        """Checks if everything is fine with empty sequences."""
        gt = gotoh.Gotoh()
        result = gt.run("../T_INPUT/empty_test.fasta", "../T_INPUT/empty_test.fasta", "../T_INPUT/6eva.txt", -4, -1,
                        True)
        (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

        assert id_seq1 == ""
        assert id_seq2 == ""
        assert seq1 == ""
        assert seq2 == ""
        assert score == 0
        assert alignments == []  # order of elements can be random!

    def test_one_empty(self):
        """Checks if everything is fine with one empty sequence."""
        gt = gotoh.Gotoh()
        result = gt.run("../T_INPUT/empty_test.fasta", "../T_INPUT/5test_seq2.fasta", "../T_INPUT/6eva.txt", -4, -1,
                        True)
        (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

        assert id_seq1 == ""
        assert id_seq2 == "Seq2"
        assert seq1 == ""
        assert seq2 == "TCCGA"
        assert score == -9
        assert alignments == [["-----", "TCCGA"]]  # order of elements is random!
