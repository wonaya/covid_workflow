#/usr/bin/python

#identify open reading frames, start/stop codons

#import glob
import os,sys
#import numpy

#os.chdir("C:/Users/mshpak/Desktop/Covid_19/April_15_MSA")

def main(file): 
    input_filename = file

    seq_file = open(input_filename, "r")
    new_seq_file = open(input_filename.split(".")[0]+"_new_name_file", "w")  #create file with new sequence names
    table_for_dict = open(input_filename.split(".")[0]+"_name_dict_file", "w")  #tab-separated table with temporary name and sequence name
    seq_name_dict = {}

    allseq = ''
    index_for_dict = 0
    seq_identify = []
    for line in seq_file:
        #line = line.upper()
        #line = line.rstrip('\n')
        if 'HCOV' in line or '>' in line:
            index_for_dict = index_for_dict + 1
            ind_name = '>seq_' + str(index_for_dict)
            seq_name_dict.update({ind_name: line.rstrip('\n')})
            new_seq_file.write(ind_name + '\n')
            table_for_dict.write(ind_name + '	' + line)
        else:
            new_seq_file.write(line)
        allseq = allseq + line

    new_seq_file.close()

if __name__ == "__main__":
	main()
