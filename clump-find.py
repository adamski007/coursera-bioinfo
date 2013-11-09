import os;
import sys;

# Exercice 4 from bio-informatics of coursera.

def find_clump():
    # 1. We got a function which does a find of a k-mers in a string -> use that one.
    # 2. Modify our algo for the find-most-k-mers, to return a list with the k-mers it-self
    #       and also the count of k-mers. [ For the ex, we need at least x k-mers, and not only the
    #       most ]
    # It should return us a list with all k-mers.
    #
    # The logic of the function :
    # 1. Check in the sub-string the most k-mers.
    # 2. Check if the k-mers is at lest present t times.
    # 3. If yes, add this k-mers to the output.
    # 4. Move on the next sub-string of the string input.