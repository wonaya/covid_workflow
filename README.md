# covid_workflow

Make sure to do the following on Stampede2

`module load tacc-singularity Rstats`

Usage :

`python main.py -s gisaid_aligned.fa -p 1-9 -m 1`

## Arguments 

`-s` : Aligned sequence in Fasta format 

`-p` : Processes to run, can use `,` or `-`, for example 1-9 runs all, 8,9 runs DNApars and CODEML

`-m` : Model to use in Codeml run, default is 0. Can be 0, 1 or 2

## Processes 
Step 1. Clustalo alignment (deprecated)

Step 2. Identify ORF boundaries

Step 3. Separate out ORFs

Step 4. Shorten sequence name (parallel)

Step 5. Align individual ORF sequences to MSA format using Clustalo (parallel)

Step 6. Remove stop codon from ORF MSA (parallel)

Step 7. Reformat Phylip (parallel)

Step 8. DNApars (parallel)

Step 9. CODEML (parallel)
