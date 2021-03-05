Samit Basnet
align.py 
README file

1)

Enter this on the terminal:
python3 align.py dnaseq1.txt dnaseq2.txt output1.txt F G
    - align.py is the program we are running and it takes 5 command line arguments
    - the first two are files with either DNA or protein sequence
    - the third command line argument is the output file
    - forth argument tells if we should use the protein matrix or dna matrix depending on our sequence
    -the last argument tells if we are doing a global or a semi-global sequence alignment

Output file: output1.txt
 Here, we are using dnaMatrix as the matrix for our scoring and doing a global alignment
 the table will be printed on the command line so it is easy to follow the scores

2) python3 align.py proteinseq1.txt proteinseq2.txt output2.txt T G 

Output file: output2.txt

Here, we are doing a global sequence alignment between two protein sequences using BLOSUM45 as the scoring matrix. 

3) python3 align.py proteinseq1.txt proteinseq2.txt output3.txt T S 

Output file: output3.txt

Semiglobal alignment of two protein sequences using BLOSUM45 scoring matrix 

4) python3 align.py dnaseq1.txt dnaseq2.txt output4.txt F S

Output file: output4.txt

Semiglobal alignment of dna sequences using dnaMatrix for scoring