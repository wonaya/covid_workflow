print "module load tacc-singularity biocontainers clustalo Rstats phylip paml"
import os,sys, argparse
from datetime import datetime

import multiprocessing
from glob import glob
    
class argChecker():
	def __init__(self, options, afterValid):
		self.options = options
		self.av = afterValid
	def check(self, x):
		if x in self.options:
			return x
		else:
			raise argparse.ArgumentTypeError("%s not a valid %s"%(x, self.av))

def clustalo(orf_no) :
    os.chdir("Outfiles/orf_"+orf_no)
    os.system("singularity run -B $PWD:/data docker://biocontainers/clustal-omega:v1.2.1_cv5 clustalo -i orf_"+str(orf_no)+"_new_name_file -o orf_"+str(orf_no)+"_aligned.phy --outfmt=phylip --threads=272 --force")
    print "orf_"+str(orf_no), "DONE"
    os.chdir("../..")

def formatR(orf_no) :
    from Bio import SeqIO
    
    os.chdir("Outfiles/orf_"+str(orf_no))
    records = SeqIO.parse("orf_"+str(orf_no)+"_aligned_rmstop.phy", "phylip")
    count = SeqIO.write(records, "orf_"+str(orf_no)+"_aligned_rmstop.fasta", "fasta")
    
    Rfile = open("convert.R", 'w')
    Rfile.write("library(ape)"+"\n")
    Rfile.write("library(seqinr)"+"\n")
    Rfile.write('f<-read.fasta("orf_'+str(orf_no)+'_aligned_rmstop.fasta")'+"\n")
    Rfile.write('write.dna(f,"orf_'+str(orf_no)+'_aligned_rmstop.phy", nbcol=1,colsep="", colw=1000000)'+"\n")
    Rfile.close()
    os.system("Rscript convert.R")
    os.system("rm -Rf convert.R")
    os.chdir("../..")

def space_out(orf_no) :
    os.chdir("Outfiles/orf_"+str(orf_no))
    outfile= open("orf_"+str(orf_no)+"_aligned_rmstop2.phy",'w')
    file =  open("orf_"+str(orf_no)+"_aligned_rmstop.phy",'r')
    lines = file.readlines()
    outfile.write(lines[0])
    for a in lines[1:]:
        if len(a.split(" ")[0]) == 5 :
            outfile.write(a.split(" ")[0])
            outfile.write("      \t")
            outfile.write(a.split(" ")[-1])
        elif len(a.split(" ")[0]) == 6 :
            outfile.write(a.split(" ")[0])
            outfile.write("     \t")
            outfile.write(a.split(" ")[-1])
        elif len(a.split(" ")[0]) == 7 :
            outfile.write(a.split(" ")[0])
            outfile.write("    \t")
            outfile.write(a.split(" ")[-1])
    outfile.close() 
    os.system("mv orf_"+str(orf_no)+"_aligned_rmstop2.phy orf_"+str(orf_no)+"_aligned_rmstop.phy")
    os.chdir("../..")

def dnapars(orf_no) :
    os.chdir("Outfiles/orf_"+str(orf_no))
    os.system("rm -Rf outfile outtree")
    outfile = open("bash.sh", 'w') 
    #outfile.write("singularity run -B $PWD:/data docker://jcuhpc/phylip:latest dnapenny << EOF\n")
    outfile.write("singularity run -B $PWD:/data docker://jcuhpc/phylip:latest dnapars << EOF\n")
    outfile.write("orf_"+str(orf_no)+"_aligned_rmstop.phy\n")
    outfile.write("5\n")
    outfile.write("Y\n")
    outfile.write("EOF\n")
    outfile.close()
    os.system("bash bash.sh > log.txt")
    print "orf_"+str(orf_no), "DONE"
    os.system("rm -Rf bash.sh")
    os.chdir("../..")

def codeml(orf_no, model) :
    os.chdir("Outfiles/orf_"+str(orf_no))
    cfile = open("codeml.ctl", 'w')
    cfile.write("seqfile = orf_"+str(orf_no)+"_aligned_rmstop.phy"+"\n")
    cfile.write("outfile = orf_"+str(orf_no)+"_dnds_output"+"\n")
    cfile.write("treefile = outtree"+"\n")
    cfile.write("noisy = 9\n")
    cfile.write("verbose = 0\n")
    cfile.write("runmode = 0\n\n")
    cfile.write("seqtype = 1\n")
    cfile.write("CodonFreq = 3\n\n")
    cfile.write("aaDist = 0\n")
    cfile.write("aaRatefile = wag.dat\n\n")
    #cfile.write("model = 0\n\n")
    cfile.write("model = "+str(model)+"\n\n")
    cfile.write("NSsites = 0\n\n")
    cfile.write("icode = 0\n")
    cfile.write("Mgene = 0\n\n")
    cfile.write("fix_kappa = 0\n")
    cfile.write("kappa = 2\n")
    cfile.write("fix_omega = 0\n")
    cfile.write("omega = .4\n")
    cfile.write("fix_alpha = 1\n")
    cfile.write("alpha = 0.\n")
    cfile.write("Malpha = 0\n")
    cfile.write("ncatG = 3\n")
    cfile.write("fix_rho = 1\n")
    cfile.write("rho = 0.\n")
    cfile.write("clock = 0\n")
    cfile.write("getSE = 0\n")
    cfile.write("RateAncestor = 0\n")
    cfile.write("Small_Diff = .5e-6\n")
    cfile.close()
    print os.getcwd() 
    os.system("rm -Rf tmp.txt")
    os.system("/work/02114/wonaya/stampede2/software/paml4.9j/src/codeml > tmp.txt")
    #os.system("singularity run -B $PWD:/data docker://biocontainers/paml:v4.9hdfsg-1-deb_cv1 codeml")
#    os.system("rm -Rf codeml.ctl")
    os.chdir("..")

def main():
    parser = argparse.ArgumentParser(description="Run through pipeline",formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-s", "--seq", metavar='FASTA',help="Aligned fasta", required=True)
    parser.add_argument("-p", "--process", help="process to run, default to run all", default="1-9", required=False)
    parser.add_argument("-m", "--model", help="PAML model for codeml, default:0", default=0, required=False)
    args = parser.parse_args()

    sys.path.append("/scratch/02114/wonaya/COVID-19_genomes/pipeline/scripts")

    process = []
    if args.process == None :
        process = range(1,10)
    else :
        if "-" in args.process : 
            process = range(int(args.process.split("-")[0]), int(args.process.split("-")[1])+1)
        if "," in args.process :
            for x in args.process.split(",") :
                process.append(int(x))
        elif len(args.process) == 1 :
            process.append(int(args.process))
    #os.system("rm -Rf Outfiles")
    if 1 in process :
        print "Step 1. Clustalo"
        #os.system("singularity run -B $PWD:/data docker://biocontainers/clustal-omega:v1.2.1_cv5 clustalo -i "+gisaid+" -o "+gisaid_out+" --outfmt=fa --threads=272 -v")
        print datetime.now()
   
    if 2 in process :
        print datetime.now()
        print "Step 2. Identify ORF boundaries"
        import MSA_ORF_Boundaries_mpi
        #MSA_ORF_Boundaries.main(args.seq)
        MSA_ORF_Boundaries_mpi.main(args.seq)
        print datetime.now()
    
    if 3 in process :
        print datetime.now()
        print "Step 3. Separate out ORFs"
        import Create_ORF_Files_NoStop
        os.chdir("Outfiles")
        Create_ORF_Files_NoStop.main()
        print datetime.now()
        os.chdir("..")
    
    if 4 in process :
        print datetime.now()
        print "Step 4. Shorten sequence name"
        import Shorten_Sequence_name
        jobs = []
        for dirs in glob("Outfiles/orf*/") :
            no = dirs.strip("Outfiles/orf_")
            s = multiprocessing.Process(target=Shorten_Sequence_name.main, args=("Outfiles/orf_"+no+"/orf_"+no, ))
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]      
        print datetime.now()
    
    if 5 in process :
        print datetime.now()
        print "Step 5. Align individual ORF sequences to MSA format"
        jobs = []
        for dirs in glob("Outfiles/orf*/") :
            no = dirs.strip("Outfiles/orf_")
            s = multiprocessing.Process(target=clustalo, args=(no, ))        
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]
        print datetime.now()
    
    if 6 in process :    
        print datetime.now()
        print "Step 6. Remove stop codon from ORF MSA"
        import remove_stop_codons
        jobs = []
        for dirs in glob("Outfiles/orf*/") :
            print dirs
            s = multiprocessing.Process(target=remove_stop_codons.main, args=(dirs.strip("Outfiles/orf_/"), ))
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]
        print datetime.now()
    
    if 7 in process :
        print datetime.now()
        print "Step 7. Reformat Phylip"
        jobs = []
        for dirs in glob("Outfiles/orf*/") :
            print dirs
            s = multiprocessing.Process(target=formatR, args=(dirs.strip("Outfiles/orf_/"), ))
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]
    
        jobs = []
        for dirs in glob("Outfiles/orf*/") :
            print dirs
            s = multiprocessing.Process(target=space_out, args=(dirs.strip("Outfiles/orf_/"), ))
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]
    
        print datetime.now()
    
    if 8 in process :
        print datetime.now()
        print "Step 8. DNApars"
        jobs = []
        for dirs in glob("Outfiles/orf*/") :
            print dirs
            s = multiprocessing.Process(target=dnapars, args=(dirs.strip("Outfiles/orf_/"), ))
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]
        print datetime.now()
    
    if 9 in process :
        print datetime.now()
        print "Step 9. CODEML"
        jobs = []
        for dirs in glob("Outfiles/orf*/") :
            s = multiprocessing.Process(target=codeml, args=(dirs.strip("Outfiles/orf_/"), args.model, ))
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]
        print datetime.now()
    
if __name__ == "__main__":
    main()
