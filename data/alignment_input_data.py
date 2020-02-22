import system.string_symbols as strings
import system.defaults as d

class AlignmentInputData:
    '''
    Stores the input data of an alignment algorithm to give easy access on it.
    '''
    global cost_function
    global gap_cost
    global gap_alpha
    global gap_beta
    global gap_opening
    global sequence_a
    global sequence_b
    global sequence_c
    global ids

    def init_needleman(self,
                       cost_function,
                       gap_cost,
                       sequence_a,
                       sequence_b,
                       sequence_c = strings.EMPTY,
                       ids = [d.SEQ_1, d.SEQ_2, d.SEQ_3]):
        '''
        Initializes the input data with values you need to use Needleman-Wunsch algorithm.

        Args:
            cost_function: PAM or BLOSUM substitution matrix you want to use
            gap_cost: value for gap costs
            sequence_a: first string you want to align with another string
            sequence_b: second string you want to align with another string
            sequence_c: third string you want to align with another string

        Returns:
            object of that class
        '''
        self.cost_function = cost_function
        self.gap_cost = gap_cost

        self.sequence_a = sequence_a
        self.sequence_b = sequence_b
        self.sequence_c = sequence_c

        self.ids = ids
        return self

    def init_gotoh(self,
                   cost_function,
                   gap_alpha,
                   gap_beta,
                   sequence_a,
                   sequence_b,
                   ids = [d.SEQ_1, d.SEQ_2]):
        '''
        Initializes the input data with values you need to use Gotoh algorithm.

        Args:
            cost_function: PAM or BLOSUM substitution matrix you want to use
            gap_alpha: basic costs for a gap
            gap_beta: costs for gap extension
            sequence_a: first string you want to align with another string
            sequence_b: second string you want to align with another string

        Returns:
            object of that class
        '''
        self.cost_function = cost_function
        self.gap_alpha = gap_alpha
        self.gap_beta = gap_beta
        self.gap_opening = gap_alpha + gap_beta
        self.sequence_a = sequence_a
        self.sequence_b = sequence_b
        self.ids = ids
        return self