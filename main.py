from algorithms.algorithm_handler import AlgorithmHandler
from system import cmd_line_parsing

# "starting point of the program"

# parses and checks if an algorithm is available
# -> if 'yes':  it starts the algorithm handler with that algorithm
# -> else:      it prints an error message

input = cmd_line_parsing.start()

if len(input) > 1:
    handler = AlgorithmHandler()
    handler.let_compute(input)
