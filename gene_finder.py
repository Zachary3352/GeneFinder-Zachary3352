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

    Didn't feel the need to add another unit test here because these two adequately test the function.
    I could add unit tests to get the complements for T and G, but this feels unnecessary with such a simple function.
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

    This function's ability is clearly visible using only the provided doctests, so I decided not to add more.
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

    Added this doctest:
    >>> rest_of_ORF("ATGGTTTCTAATTAA")
    'ATGGTTTCTAAT'

    Added the above doctest to demonstrate the third possible stop codon.
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

    Added this doctest:
    >>> find_all_ORFs_oneframe("ATGCATTATTGCCGTGGCTGAGTTATGATTCTTCTCTAA")
    ['ATGCATTATTGCCGTGGC', 'ATGATTCTTCTC']

    Added the above doctest as an extra example.
    """
    i = 0
    current_orf = []
    all_found_orfs = []
    while i < len(dna):
        if dna[i:i+3] == codons[3][0]:
            current_orf = rest_of_ORF(dna[i:])
            all_found_orfs.append(current_orf)
            i = i + len(current_orf) - 3
        i = i+3
    return all_found_orfs

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

    Didn't feel the need to add another doctest here because the above test already adequately tests the function.
    It shows non-nested ORFs, but includes ORFs that end at the same point.
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

    Added this doctest:
    >>> find_all_ORFs_both_strands("ATGCATTATTGCCGTGGCTGAGTTATGATTCTTCTCTAA")
    ['ATGCATTATTGCCGTGGC', 'ATGATTCTTCTC', 'ATGCAT']

    Added the above doctest as an additional, more complicated test.
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

    Didn't feel the need to add an additional doctest because the provided doctest includes shorter ORFs than the
    ORFs the doctest finds -- the doctest finds the longest ORF.
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
        returns: the maximum length longest ORF

        I don't believe it's possible to add a doctest that always works to this function because the maximum
        length longest ORF could change based on how many trials are used.
        """

    longest_found_ORF = 0

    for scramble in range(num_trials):
        shuffled_string = shuffle_string(dna)
        found_ORF = len(longest_ORF(shuffled_string))
        if found_ORF > longest_found_ORF:
            longest_found_ORF = found_ORF
    return longest_found_ORF

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

        Added this doctest:
        >>> coding_strand_to_AA("ATGAAAGACTGTCCCATTTGA")
        'MKDCPI|'

        Added the above doctest to test using stop codons in coding_strand_to_AA.
    """


    DNA_list = []
    AA_list = []

    for i in range(len(dna)//3):
        DNA_list.append(dna[(i+1)*3-3:(i+1)*3])

    for i in range(len(DNA_list)):
        AA_list.append(aa_table[DNA_list[i]])

    return ''.join(AA_list)

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

        Added this doctest:
        >>> gene_finder("ATAGATCCAGTAAGATAGATCCAGTAAGCATAAGTGTGTTAGACAGATCCAGTAAGCATAAGTGTGTTAGACATGAAGCCCTAGCAGATCCAGTAAGCATAAGTGTGTTAGACATAGATCCAGTAAGCATAAGTGTGTTAGACGTAATAGCACACATAGATCCAGTAAGCATAAGTGTGTTAGACAGATCCAGTAAGCATAAGTGTGTTAGACATGAAGCCCTAGCAGATCCAGTAAGCATAAGTGTGTTAGACATAGATCCAGTAAGCATAAGTGTGTTAGACGTAATAGCACACATAGATCCAGTAAGCATAAGTGTGTTAGACAGATCCAGTAAGCATAAGTGTGTTAGACATGAAGCCCTAGCAGATCCAGTAAGCATAAGTGTGTTAGACATAGATCCAGTAAGCATAAGTGTGTTAGACGTAATAGCACACATAGATCCAGTAAGCATAAGTGTGTTAGACAGATCCAGTAAGCATAAGTGTGTTAGACATGAAGCCCTAGCAGATCCAGTAAGCATAAGTGTGTTAGACATAGATCCAGTAAGCATAAGTGTGTTAGACGTAATAGCACACCATAAGTGTGTTAGACAGATCCAGTAAGCATAAGTGTGTTAGACATGAAGCCCTAGCAGATCCAGTAAGCATAAGTGTGTTAGACATAGATCCAGTAAGCATAAGTGTGTTAGACGTAATAGCACAC")
        ['MLTGSMSNTLMLTGSARASCLTHLCLLDLSNTLMLTGSMCAITSNTLMLTGSMSNTLMLTGSARASCLTHLCLLDLSNTLMLTGSMCAITSNTLMLTGSMSNTLMLTGSARASCLTHLCLLDLSNTLMLTGSMCAITSNTLMLTGSMSNTLMLTGSARASCLTHLCLLDLSNTLMLTGSILLDL']

        Added the above doctest to "test" my final gene_finder. However, because so much DNA is needed to "get past" the
        threshold, I needed to use my function to add the doctest. So I'm not really sure that it's a proper doctest.
    """

    orfs_to_code = []
    AA_sequences = []
    threshold = longest_ORF_noncoding(dna, 1500)
    #print(threshold)
    all_ORFs_both_strands = find_all_ORFs_both_strands(dna)
    #print(all_ORFs_both_strands)
    for orf in all_ORFs_both_strands:
        if len(orf) > threshold:
            orfs_to_code.append(orf)
    #print(orfs_to_code)
    for ORF in orfs_to_code:
        AA_sequences.append(coding_strand_to_AA(ORF))
    return AA_sequences

from load import load_seq
dna = load_seq("./data/X73525.fa")

print(gene_finder(dna))

# if __name__ == "__main__":
#     import doctest
#     doctest.testmod()
#     doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(), verbose=True)
