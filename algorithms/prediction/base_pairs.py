ADENIN = "A"
CYTOSIN = "C"
GUANIN = "G"
URACIL = "U"

def complementary(base_a, base_b):
    """
    Defines which bases are complementary and which not.

    Args:
        base_a: first base
        base_b: second base

    Returns:
        True if the bases are complementary and False if not
    """
    if (base_a == ADENIN and base_b == URACIL) or (base_a == URACIL and base_b == ADENIN):
        return True
    elif (base_a == CYTOSIN and base_b == GUANIN) or (base_a == GUANIN and base_b == CYTOSIN):
        return True
    #elif (base_a == URACIL and base_b == GUANIN) or (base_a == GUANIN and base_b == URACIL):
        #return True

    return False