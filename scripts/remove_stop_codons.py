#/usr/bin/python

#identify open reading frames, start/stop codons

import glob
import os
import sys
#import numpy


#input_filename = sys.argv[1]
# use input filename toyaligned.fa

#break string into triplets
def string_to_codon(s1):
    lens = len(s1)
    triplets = [(s1[i:i+3]) for i in range(0, lens, 3)]
    return(triplets) 

def replace_stops_gaps(s1, Stop):
    codons = string_to_codon(s1)
    for scod in Stop:
        codons = [w.replace(scod, '---') for w in codons]
        newseq = ''.join(codons)
    return(newseq)

def string_to_pieces(s1, l_size):
    lens = len(s1)
    in_pieces = [(s1[i:i+l_size]) for i in range(0, lens, l_size)]
    return(in_pieces)

def main(orf_no) :
    os.chdir("orf_"+str(orf_no))
    input_filename = "orf_"+str(orf_no)+"_aligned.phy"
    Start = ['ATG', 'atg']
    Stop = ['TAA','taa', 'TAG', 'tag', 'TGA', 'tga']

    al_file = open(input_filename, "r")
    new_al_file = open(input_filename.split(".phy")[0]+"_rmstop.phy", "w")

    header = al_file.readline()

    seq_concat = {} #dictionary for contatenated sequences
    seq_tags = {} #dictioanry for sequence names
    pos = 0
    for line in al_file:
        line = line.rstrip('\n')
        if len(line) > 0:
            pos = pos + 1
            temp_vec = line.split()
            if len(temp_vec) > 1:
                seq = temp_vec[1]
                tags = temp_vec[0]
                seq_concat[str(pos)] = seq
                seq_tags[str(pos)] = tags
            else:
                #pos = pos + 1
                seq = temp_vec[0]
                seq_concat[str(pos)] = seq_concat[str(pos)] + seq
        else:
            pos = 0

    for seqno in seq_concat:
        stemp = seq_concat[seqno]
        cleaned_up = replace_stops_gaps(stemp, Stop)
        seq_concat[seqno] = string_to_pieces(cleaned_up, 50)
        
    #number of pieces per sequence and number of sequences, respectively
    npieces = len(seq_concat['1'])
    nseqs = len(seq_concat)

    #write in blocks to mimic phylip file
    new_al_file.write(header)

    for i in range(0,npieces):
        countseqs = 0 #to include blank line between nseqs
        for j in range(0,nseqs):
            if i == 0:
                temp_tag = seq_tags[str(j+1)]
                temp_seq = seq_concat[str(j+1)]
                lshift = len(temp_tag)
                new_al_file.write(temp_tag + ' '*(10-lshift) + temp_seq[i] + '\n' )
            else:
                temp_seq = seq_concat[str(j+1)]
                new_al_file.write(temp_seq[i] + '\n')
        new_al_file.write('\n')

    new_al_file.close()
    al_file.close()
    os.chdir("..")

if __name__ == "__main__":
	main()

