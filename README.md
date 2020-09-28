# covid_workflow

Make sure to do the following on Stampede2

`module load tacc-singularity biocontainers clustalo Rstats phylip paml`

Usage :

`python main.py -s gisaid_aligned.fa`

Step 1. Clustalo alignment (deprecated)
Step 2. Identify ORF boundaries
Step 3. Separate out ORFs
Step 4. Shorten sequence name (parallel)
Step 5. Align individual ORF sequences to MSA format using Clustalo (parallel)
Step 6. Remove stop codon from ORF MSA (parallel)
Step 7. Reformat Phylip (parallel)
Step 8. DNApars (parallel)
Step 9. CODEML (parallel)
