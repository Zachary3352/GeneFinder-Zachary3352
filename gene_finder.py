# -*- coding: utf-8 -*-
"""

@author: Zachary Sherman (Zachary3352)

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###

def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """

    if nucleotide == "A":
        return "T"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "C":
        return "G"
    elif nucleotide == "G":
        return "C"
    else:
        return "Failed, check nucleotide letters!"

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    reversed = []
    for nucleotide in dna:
        backwards = dna[::-1] # Found out how to do this on https://www.educative.io/edpresso/how-do-you-reverse-a-string-in-python
    for nucleotide in range(len(dna)):
        reversed_nucleotide = get_complement(backwards[nucleotide])
        reversed.append(reversed_nucleotide)
    final_reversed_string = ''.join(reversed)
    return final_reversed_string

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGTGCCC")
    'ATGTGCCC'
    """

    i = 0
    while i < len(dna):
        if dna[i:i+3] == codons[10][0] or dna[i:i+3] == codons[10][1] or dna[i:i+3] == codons[10][2]:
            return dna[:i]
            break
        i = i+3
    return dna

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    i = 0
    orf1 = []
    while i < len(dna):
        if dna[i:i+3] == codons[3][0]:
            orf1.append(rest_of_ORF(dna[i:]))
            i = i + len(orf1) - 3
        i = i+3
    return orf1

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    i = 0
    all_orfs = []
    while i < 3:
        all_orfs = all_orfs + find_all_ORFs_oneframe(dna[i:])
        i = i+1
    return all_orfs

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    dna_reversed = get_reverse_complement(dna)
    all_ORFs_both_strands = find_all_ORFs(dna) + find_all_ORFs(dna_reversed)
    return all_ORFs_both_strands

# WEEK 2 BEGINS HERE

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """

    all_orfs = find_all_ORFs_both_strands(dna)
    if len(all_orfs) > 0:
        longest_ORF = max(all_orfs, key=len)
    else:
        longest_ORF = []
    #print(longest_ORF)
    return longest_ORF

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    longest_found_ORF = 0

    for scramble in range(num_trials):
        shuffled_string = shuffle_string(dna)
        found_ORF = len(longest_ORF(shuffled_string))
        if found_ORF > longest_found_ORF:
            longest_found_ORF = found_ORF
    return longest_found_ORF

print(longest_ORF_noncoding("ATATGAAGCCCTAGCAGATCCAGTAAGCATAAGTGTGTTAGACATGTAATAGCACAC", 200)) # This works yay!

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

#if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    #doctest.run_docstring_examples(longest_ORF, globals(), verbose=True)
