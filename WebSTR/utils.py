################# Motif complement ###################
def motif_complement(motif):
    rev_comp = ""
    for nuc in reversed(motif):
        if nuc == "A":
            rev_comp +="T"
        elif nuc == "T":
            rev_comp +="A"
        elif nuc == "C":
            rev_comp +="G"
        elif nuc == "G":
            rev_comp +="C"
        else:
            rev_comp +="-"

    result = motif +"/"+ rev_comp
   
    return result