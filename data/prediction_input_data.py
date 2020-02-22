class PredictionInputData:
    '''
    Stores the input data of a prediction algorithm to give easy access on it.
    '''

    global loop_length
    global sequence
    global id

    def init_nussinov(self, loop_length, sequence, id):
        '''
        Initializes the input data with values you need to use Needleman-Wunsch algorithm.

        Args:
            loop_length: minimum distance between two bases
            sequence: first string you want to align with another string

        Returns:
            object of that class
        '''
        self.loop_length = loop_length
        self.sequence = sequence
        self.id = id
        return self