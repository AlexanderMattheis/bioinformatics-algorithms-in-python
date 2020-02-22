import unittest

import nussinov_alg
from prakt.nussinov import NussinovBase

class TestMethodsNussinov(unittest.TestCase):
    """Test class to test Nussinov algorithm."""

    def test_instance(self):
        """Check inheritance."""
        assert issubclass(nussinov_alg.Nussinov, NussinovBase)
        assert isinstance(nussinov_alg.Nussinov(), NussinovBase)

    def test_8(self):
        """Short sequence test. Test taken from Bioinformatics exercise sheet."""
        nussinov = nussinov_alg.Nussinov()
        result = nussinov.run("../T_INPUT/8test_seq.fasta")
        (id_seq, seq, dot_bracket_notation) = result

        assert id_seq == "RNA1"
        assert seq == "CCUGUAAG"
        assert dot_bracket_notation == "((.)().)"

    def test_9(self):
        """Short sequence test. Test taken from my report."""
        nussinov = nussinov_alg.Nussinov()
        result = nussinov.run("../T_INPUT/9test_seq.fasta")
        (id_seq, seq, dot_bracket_notation) = result

        assert id_seq == "RNA_Seq"
        assert seq == "CGAAU"
        assert dot_bracket_notation == "()(.)"

    def test_10(self):
        """Short sequence test. Test taken from my presentation."""
        nussinov = nussinov_alg.Nussinov()
        result = nussinov.run("../T_INPUT/10test_seq.fasta")
        (id_seq, seq, dot_bracket_notation) = result

        assert id_seq == "RNA_Seq"
        assert seq == "GAUC"
        assert dot_bracket_notation == "(())"

    def test_empty(self):
        """Checks if everything is fine with an empty sequence."""
        nussinov = nussinov_alg.Nussinov()
        result = nussinov.run("../T_INPUT/empty_test.fasta")
        (id_seq, seq, dot_bracket_notation) = result

        assert id_seq == ""
        assert seq == ""
        assert dot_bracket_notation == ""

    def test_very_short(self):
        """Checks if everything is fine with a sequence of lenght one."""
        nussinov = nussinov_alg.Nussinov()
        result = nussinov.run("../T_INPUT/very_short_test.fasta")
        (id_seq, seq, dot_bracket_notation) = result

        assert id_seq == "RNA_one"
        assert seq == "A"
        assert dot_bracket_notation == "."

    def test_many(self):
        """Checks if everything is fine with a sequence of not matching characters."""
        nussinov = nussinov_alg.Nussinov()
        result = nussinov.run("../T_INPUT/multi_not_matching_test.fasta")
        (id_seq, seq, dot_bracket_notation) = result

        assert id_seq == "RNA_many"
        assert seq == "ACACACCAC"
        assert dot_bracket_notation == "........."
