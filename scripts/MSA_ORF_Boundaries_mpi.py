#/usr/bin/python

#identify open reading frames, start/stop codons

import glob
import os
import sys
#import numpy
import multiprocessing

#os.chdir("C:/Users/mshpak/Desktop/Covid_19")

#input_filename = sys.argv[1]
# use input filename toyaligned.fa

Start = ['ATG', 'atg']
Stop = ['TAA','taa', 'TAG', 'tag', 'TGA', 'tga']

#take column of array of arrays (matrix)
def column(my_matrix, i):
    return [row[i] for row in my_matrix]


#what fraction of sequence is N or X
def frac_unknown(seq, letter):
    seq_len = len(seq)
    n_letter = seq.count(letter)
    return(n_letter/seq_len)


#get rank based on ith element in row of a matrix
def reorder_by_index(my_matrix, i): 
    get_col = column(my_matrix, i)
    sorted_seq = sorted(get_col)
    sorted_indices = [get_col.index(v) for v in sorted_seq]
    return(sorted_indices)


#find all instances of mystring in mylist
def find_locs(mylist,mystring):
    v1 = [i for i,x in enumerate(mylist) if x == mystring]
    return(v1)


def find_all_locs(mylist, list_of_strings):
    v1 = [i for i,x in enumerate(mylist) if x in list_of_strings]
    return(v1)


def first_larger_element(my_list, to_compare):
    if(max(my_list) > to_compare):
        return(next(x[0] for x in enumerate(my_list) if x[1] > to_compare))
    else:
        pass

#create list of possible start and stop codons, as well as the length
#identify the shortest using list of stops
#output list of [start, stop], [length]
def orf_boundaries(codon_list):
    orf_starts = find_locs(codon_list, 'ATG')
    orf_stops = find_all_locs(codon_list, Stop)
    start_stop = []
    stop_list = []
    for posit in orf_starts:
        temp_ind = first_larger_element(orf_stops,posit)
        if type(temp_ind) is int:
            temp_stop = orf_stops[temp_ind]
            stop_list = stop_list + [temp_stop]
            start_stop = start_stop + [[posit, temp_stop, temp_stop-posit]]
        else:
            pass
    stop_list = sorted(list(set(stop_list)))
    final_orf_boundaries = []
    for stops in stop_list:
        pos_stops_temp = [i for i,x in enumerate(start_stop) if x[1]== stops][0]
        final_orf_boundaries = final_orf_boundaries + [start_stop[pos_stops_temp]]
    return(final_orf_boundaries)

#check if codon contain N
def translate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            #if 'N' in codon or 'K' in codon:
            if any(lets not in ['A','C','G','T'] for lets in list(set(codon))):
                protein += 'X'
            else:
                protein+= table[codon] 
    return protein 

def separate(name_temp, seq_ind, seqs) :
    name_temp = name_temp.replace('>','')
    name_temp = name_temp.replace(' ','')
    seq_temp = seqs.replace(name_temp, '') #remove name of sequences
    seq_temp = seq_temp.replace('-','')
    seq_temp = seq_temp.replace(' ','')
    reverse_seq_temp = seq_temp[::-1]
    seq_length = len(seq_temp)
    orf_1 = seq_temp
    orf_2 = seq_temp[1:seq_length]
    orf_3 = seq_temp[2:seq_length]
    orf_4 = reverse_seq_temp
    orf_5 = reverse_seq_temp[1:seq_length]
    orf_6 = reverse_seq_temp[2:seq_length]
    codons_1 = []
    codons_2 = []
    codons_3 = []
    for i in range(0, seq_length + 3, 3):
        codon_1 = [seq_temp[i:i+3]]
        codon_2 = [seq_temp[i+1:i+4]]
        codon_3 = [seq_temp[i+2:i+5]]
        codons_1 = codons_1 + codon_1
        codons_2 = codons_2 + codon_2
        codons_3 = codons_3 + codon_3
    codons_4 = []
    codons_5 = []
    codons_6 = []
    for i in range(0, seq_length + 3, 3):
        codon_4 = [reverse_seq_temp[i:i+3]]
        codon_5 = [reverse_seq_temp[i+1:i+4]]
        codon_6 = [reverse_seq_temp[i+2:i+5]]
        codons_4 = codons_4 + codon_4
        codons_5 = codons_5 + codon_5
        codons_6 = codons_6 + codon_6

    orf_bound_1 = orf_boundaries(codons_1)
    orf_keep_1_100 = [vec for vec in orf_bound_1 if vec[2]>100]
    orf_bound_2 = orf_boundaries(codons_2)
    orf_keep_2_100 = [vec for vec in orf_bound_2 if vec[2]>100]
    orf_bound_3 = orf_boundaries(codons_3)
    orf_keep_3_100 = [vec for vec in orf_bound_3 if vec[2]>100]
    orf_bound_4 = orf_boundaries(codons_4)
    orf_keep_4_100 = [vec for vec in orf_bound_4 if vec[2]>100]
    orf_bound_5 = orf_boundaries(codons_5)
    orf_keep_5_100 = [vec for vec in orf_bound_5 if vec[2]>100]
    orf_bound_6 = orf_boundaries(codons_6)
    orf_keep_6_100 = [vec for vec in orf_bound_6 if vec[2]>100]
    All_Frames = [codons_1,codons_2,codons_3]
    All_ORF = [orf_keep_1_100,orf_keep_2_100,orf_keep_3_100]
    Reverse_Frames = [codons_4,codons_5,codons_6]
    Reverse_ORF = [orf_keep_4_100,orf_keep_5_100,orf_keep_6_100]
    DNA_Seqs = []
    AA_Seqs = []
    frame_index = 0

    ORF_To_Order = orf_keep_1_100 + orf_keep_2_100 + orf_keep_3_100
    Reverse_ORF_To_Order = orf_keep_4_100 + orf_keep_5_100 + orf_keep_6_100
    for frame in All_ORF:
        for orf in frame:
            codons = All_Frames[frame_index]
            temp_codons = codons[orf[0]:orf[1]+1]
            temp_dna_seq = ''.join(temp_codons)
            frac_N = frac_unknown(temp_dna_seq, 'N')  #proportion of unknown nucleotides
            if frac_N < 0.25:
                DNA_Seqs = DNA_Seqs + [temp_dna_seq]
                temp_aa_seq = translate(temp_dna_seq)
                AA_Seqs = AA_Seqs + [temp_aa_seq]
                #print(ORF_To_Order)
                #print(orf)
            else:
                ORF_To_Order = [arr for arr in ORF_To_Order if orf != arr]
        frame_index = frame_index + 1

    DNA_Seqs = [DNA_Seqs[index] for index in reorder_by_index(ORF_To_Order,0)]
    AA_Seqs = [AA_Seqs[index] for index in reorder_by_index(ORF_To_Order,0)]
    fname_orfs = 'Corona_DNA_ORF' + name_temp + '.fa'
    fname_orfs = fname_orfs.replace('_','-')
    fname_orfs = fname_orfs.replace('|','_')
    fname_orfs = fname_orfs.replace('/','_')
    fname_peptides = 'Corona_Peptide_' + name_temp + '.fa'
    fname_peptides = fname_peptides.replace('_','-')
    fname_peptides = fname_peptides.replace('|','_')
    fname_peptides = fname_peptides.replace('/','_')
    corona_orfs = open(fname_orfs, "w")
    corona_peptides = open(fname_peptides, "w")
    pep_ind = 0
    for peptides in AA_Seqs:
        corona_peptides.write('>' + 'ORF_Peptide_' + str(pep_ind + 1))
        corona_peptides.write('\n')
        corona_peptides.write(peptides)
        corona_peptides.write('\n')
        corona_peptides.write('\n')
        pep_ind = pep_ind + 1
    corona_peptides.close()
    dna_ind = 0
    for dna_seq in DNA_Seqs:
        corona_orfs.write('>' + 'ORF_' + str(dna_ind + 1))
        corona_orfs.write('\n')
        corona_orfs.write(dna_seq)
        corona_orfs.write('\n')
        corona_orfs.write('\n')
        dna_ind = dna_ind + 1
    corona_orfs.close()
    seq_ind = seq_ind + 1


def main(seq): 
    seq_file = open(seq, "r")


    allseq = ''
    seq_identify = []
    for line in seq_file:
        line = line.upper()
        line = line.rstrip('\n')
        if 'HCOV' in line:
            seq_identify = seq_identify + [line]
        else:
            pass
        allseq = allseq + line

    separate_seqs = allseq.split('>')
    separate_seqs = separate_seqs[1:len(separate_seqs)]  #remove '' at beginning

    if not os.path.isdir("Outfiles"):
        os.mkdir("Outfiles")
    os.chdir("Outfiles")

    seq_ind = 0
    jobs = []
    for seqs in separate_seqs:
        name_temp = seq_identify[seq_ind]
        s = multiprocessing.Process(target=separate, args=(name_temp, seq_ind, seqs, ))
        jobs.append(s)
        s.start()
        seq_ind = seq_ind + 1
    os.chdir("..")
    [x.join() for x in jobs]
  
if __name__ == "__main__":
	main()
