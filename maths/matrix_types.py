# stores paths to used substitution matrices
__TXT = ".txt"
__INPUT = "INPUT/"

MATRICES = ["--blosum62", "--pam250"]
MATRICES_SHORT = ["-b", "-p"]

BLOSUM62 = MATRICES[0];
PAM250 = MATRICES[1];

BLOSUM62_SHORT = MATRICES_SHORT[0];
PAM250_SHORT = MATRICES_SHORT[1];

BLOSUM62_PATH = __INPUT + BLOSUM62.replace("--", "") + __TXT;
PAM250_PATH = __INPUT + PAM250.replace("--", "") + __TXT;
