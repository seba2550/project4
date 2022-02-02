# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/Tursiops_truncatus_BRD2.fa")

    sub_mat = "./substitution_matrices/BLOSUM62.mat"
    nw = NeedlemanWunsch(sub_mat, -10, -1)

    seqs_dict = {"Gallus gallus": gg_seq, "Mus musculus": mm_seq, "Balaeniceps rex": br_seq, "Tursiops truncatus": tt_seq}
    scores_dict = {"Gallus gallus": None, "Mus musculus": None, "Balaeniceps rex": None, "Tursiops truncatus": None}
    for species, seq in seqs_dict.items():
        scores_dict[species] = nw.align(hs_seq, seq)[0]

    # This was a really cool implementation but then I remembered dictionaries cannot be "traditionally" sorted. Anyway here's a lambda function to actually sort the dictionary by values, because I'm in too deep already to turn back on this approach
    scores_dict = {k: v for k, v in sorted(scores_dict.items(), key=lambda x: x[1], reverse = True)} 
    print("These are the species in descending order of most similar to the BRD human sequence: " , scores_dict.keys())

    print("And now the scores for the alignments with the BRD human sequence (again in descending order: " , scores_dict.values())

if __name__ == "__main__":
    main()
