import os,sys

outfile= open(str(sys.argv[1].split(".phy")[0])+"_2.phy",'w')
file =  open(sys.argv[1],'r')
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
    
