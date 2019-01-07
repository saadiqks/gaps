# gaps
This program runs the GAPS method to verify an upper bound for an f(m,s) problem.

## Instructions
This program can be run in two ways. In the first, the inputs are m, s, and the numerator and denominator of the upper bound. The proof will be in a file called `output.txt`. For example, one can run `python gaps.py 31 19 54 133` to verify f(31,19) <= 54/133. In the second, the input is a text file where each line is four space delimited numbers which each represent a f(m,s) problem. If `python gaps.py try.txt` is run and the first line of `try.txt` is `31 19 54 133`, the proof for f(31,19) will be contained in the file `31.19.ERIK.txt`, since this is an ERIK case. This follows for each line of `try.txt`.
