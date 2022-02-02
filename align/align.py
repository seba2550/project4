# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        # TODO: Fill in the Needleman-Wunsch Algorithm below
        to perform global sequence alignment of seqA and seqB
        and return a tuple with the following format
        (alignment score, seqA alignment, seqB alignment)
        Also, write up a docstring for this function using the
        _read_sub_matrix as an example.
        Don't forget to comment your code!
        """
        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for pointers used in backtrace procedure
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        ###### TODO Implement the global sequence alignment here
        # Get the lengths of both sequences
        n = len(self._seqA) 
        m = len(self._seqB)

        # Fill in the upper left corner of our alignment matrix
        self._align_matrix[0, 0] = 0 
        
        # Fill in the upper left corner of the gap A matrix
        self._gapA_matrix[0, 1:] = self.gap_open

        # Fill in the upper left corner of the gap B matrix
        self._gapB_matrix[1:, 0] = self.gap_open

        # Plug in initial pointer values for our backtracing matrices. We'll use these later on to trace back through the alignment. A 0 means we put a gap in A, whereas a 2 means we put a gap in B.
        self._back_A[0, 1:] = 0  
        self._back_B[1:, 0] = 2  


        # Now fill in the rest of the three matrices. 
        for i in range(1, n + 1):
            for j in range(1, m + 1):

                # For each of the three matrices calculate all possible moves
                alignment_match = self.sub_dict[(seqA[i - 1], seqB[j - 1])] + np.array([self._gapA_matrix[i - 1, j - 1], self._align_matrix[i - 1, j - 1], self._gapB_matrix[i - 1, j - 1]])
                insert_gapA = self.gap_extend + np.array([self._gapA_matrix[i, j - 1], self.gap_open + self._align_matrix[i, j - 1], self.gap_open + self._gapB_matrix[i, j - 1]])
                insert_gapB = self.gap_extend + np.array([self.gap_open + self._gapA_matrix[i - 1, j], self.gap_open + self._align_matrix[i - 1, j], self._gapB_matrix[i - 1, j]])

                # Now pick the best move (max value, essentially) and add it to its respective matrix at the corresponding index
                self._align_matrix[i, j] = np.max(alignment_match)
                self._gapA_matrix[i, j] = np.max(insert_gapA)
                self._gapB_matrix[i, j] = np.max(insert_gapB)

                # We then use numpy's argmax function to return the indices of the best moves we calculated for each matrix. Those indices are added to our pointer matrices so that we can backtrace through them
                self._back[i, j] = np.argmax(alignment_match)
                self._back_A[i, j] = np.argmax(insert_gapA)
                self._back_B[i, j] = np.argmax(insert_gapB)
                
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        # TODO Implement the traceback procedure method below
        based on the heuristic you implement in the align method.
        The traceback method should return a tuple of the alignment
        score, the seqA alignment and the seqB alignment respectively.
        """
        # Implement this method based upon the heuristic chosen in the align method above.
        
        # This gets us the index for the bottom right corner of the "best" matrix (i.e.: the one with the highest score at the bottom, or the optimal solution)
        i = len(self._seqA)
        j = len(self._seqB)

        # Get the starting point for our traceback. It'll be the highest scoring bottom corner of any one of our three matrices.
        backtrace_start = (self._gapA_matrix[i, j], self._align_matrix[i, j], self._gapB_matrix[i, j])
        self.alignment_score = np.max(backtrace_start) # This bottom corner is also, by definition, our alignment score. Save it so that we can return it at the end.
        path = np.argmax(backtrace_start) # Get the index of this corner to access the pointer matrices

        while i > 0 and j > 0: # This loop continues until we reach 0 on any of the indices.
            if path == 0: # If the value is 0, we inserted a gap in A. Add a dash to the seq A align object and update the path variable
                self.seqA_align += "-"
                self.seqB_align += self._seqB[j - 1]
                path = self._back_A[i, j]
                j -= 1

            elif path == 1: # If the value is 1, we matched the two amino acids. Add both the amino acids to the seq align objects and update the path variable
                self.seqA_align = self._seqA[i - 1] + self.seqA_align
                self.seqB_align = self._seqB[j - 1] + self.seqB_align
                path = self._back[i, j]
                i -= 1
                j -= 1

            if path == 2: # Finally, if the value is 2, we inserted a gap in B. Add a dash to the seq B align object and again update the path variable
                self.seqA_align = self._seqA[i - 1] + self.seqA_align
                self.seqB_align = "-" + self.seqB_align
                path = self._back_B[i, j]
                i -= 1
            
        
        return(self.alignment_score, self.seqA_align, self.seqB_align) # We've finished backtracing, return the alignment score and the two aligned sequences


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
