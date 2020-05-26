#! /bin/sh
#$ -S /bin/bash
#$ -pe smp 2
#$ -cwd
#$ -N neoantigen_frameshift
#$ -e log/
#$ -o log/

. ~/.bashrc
. ~/.bash_profile
. ./frameshift.config

exec 1>log_neoantigen_frameshift.out
exec 2>log_neoantigen_frameshift.err

# make TCGA-ID file: one column with TCGAID 
# make tcgaid.genelist: run type-make-genelist.sh script
sh frameshift.make-genelist.sh $SOURCE $INPUT

# Extract coordinates of frameshift mutations.
python frameshift.get-coordinate.py $SOURCE

# Get CDS and UTR sequences.
Rscript $BIOMART frameshift.coordinate.txt coding frameshift.coding.fasta
Rscript $BIOMART frameshift.coordinate.txt 3utr frameshift.3UTR.fasta

# Mutate cds, translate, and get frameshift peptide: 8 AA before the mutation site till the stop codon.
python frameshift.cdna-to-peptide.py

#rm frameshift.coding.fasta
#rm frameshift.3UTR.fasta
#rm frameshift.strand.txt

# pre-process peptide list.
peptide=""
while read line
do
if [[ $line == ">"* ]]; then
echo $peptide >> frameshift.peptide.fasta.tmp
echo $line >> frameshift.peptide.fasta.tmp
peptide=""
else
peptide="$peptide$line"
fi
done < frameshift.peptide.fasta
mv frameshift.peptide.fasta.tmp frameshift.peptide.fasta

# split whole peptide file on peptide list per patient.
mkdir -p frameshift-peptides-per-patient
while read line
do
ID=$(echo $line | awk '{print $1}')
SeqID=$(echo $line | awk '{print $9}')
header=$(grep -w "$SeqID" frameshift.peptide.fasta)
peptide=$(grep -A1 -w "$SeqID" frameshift.peptide.fasta | grep -v "$SeqID")
echo $header >> "frameshift-peptides-per-patient/$ID.peptide.fasta"
echo $peptide >> "frameshift-peptides-per-patient/$ID.peptide.fasta"
done < frameshift.coordinate.txt

# put patients on smaller folders:
#mkdir frameshift-peptides-per-patient-501-1000
#ls frameshift-peptides-per-patient | head -500 | while read line; do mv frameshift-peptides-per-patient/$line frameshift-peptides-per-patient-501-1000/; done


# Predict neoantigens:
python frameshift.peptide-to-neoantigen.py "frameshift-peptides-per-patient" "$HLA" "frameshift.$PROJECT.antigen.txt"





