#/usr/bin/python

#identify open reading frames, start/stop codons

import glob
import os
import sys
#import numpy
import os

def main():
    filenames = os.listdir(".")
    orf_filenames = [sname for sname in filenames if "Corona-DNA" in sname]
    pep_filenames = [sname for sname in filenames if "Corona-Peptide" in sname]


    #remove instances of >8 or <8 orfs
    for orf_file in orf_filenames:
        orf_data = open(orf_file, "r")
        all_dat = orf_data.read()
        occurrences = all_dat.count('>ORF')
        orf_data.close()
        if(occurrences != 8):
            os.remove(orf_file)
            orf_filenames = [ele for ele in orf_filenames if orf_file != ele]

    for pep_file in pep_filenames:
        pep_data = open(pep_file, "r")
        all_dat = pep_data.read()
        occurrences = all_dat.count('>ORF')
        pep_data.close()
        if(occurrences != 8):
            os.remove(pep_file)
            pep_filenames = [ele for ele in pep_filenames if pep_file != ele]

    for orf_file in orf_filenames:
        orf_dat = open(orf_file, "r")
        orf_index = 0
        for line in orf_dat:
            line = line.rstrip('\n')
            if '>ORF' in line:
                line_new = next(orf_dat)
                if not os.path.isdir('orf_'+str(orf_index + 1)) :
                    os.mkdir('orf_'+str(orf_index + 1))
                os.chdir('orf_'+str(orf_index + 1))
                orf_temp = open('orf_' + str(orf_index + 1), "a+")  #a+ indicates append or create new file
                orf_temp.write('>' + orf_file)
                orf_temp.write('\n')
                orf_temp.write(line_new[:-4])
                orf_temp.write('\n')
                orf_temp.write('\n')
                orf_index = orf_index + 1
                orf_temp.close()
                os.chdir("..")
        orf_dat.close()


    for pep_file in pep_filenames:
        pep_dat = open(pep_file, "r")
        pep_index = 0
        for line in pep_dat:
            line = line.rstrip('\n')
            if '>ORF' in line:
                line_new = next(pep_dat)
                if not os.path.isdir('pep_'+str(pep_index + 1)) :   
                    os.mkdir('pep_'+str(pep_index + 1))
                os.chdir('pep_'+str(pep_index + 1))
                pep_temp = open('pep_' + str(pep_index + 1), "a+")  #a+ indicates append or create new file
                pep_temp.write('>' + pep_file)
                pep_temp.write('\n')
                #print(line_new[:-1])
                pep_temp.write(line_new[:-2])
                pep_temp.write('\n')
                pep_temp.write('\n')
                pep_index = pep_index + 1
                pep_temp.close()
                os.chdir("..")
        pep_dat.close()

if __name__ == "__main__":
	main()
