from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import sys
import subprocess
import os

def AddPeptideRecordToOutput(mtaa, sequence, header, output):
        flag = mtaa - 9 
        if flag > 0:
                record = SeqRecord(sequence.translate(to_stop=True)[flag : ], id=header, description="")
        else:
                record = SeqRecord(sequence.translate(to_stop=True)[0 : ], id=header, description="")
        output.append(record)

def grep_seq(filename, arg):
	cmd="grep -A1 "+arg+" "+filename+" | grep -v "+arg
	try:
		result = subprocess.check_output(cmd, shell=True).rstrip()
		result = Seq(result, alphabet=IUPAC.ambiguous_dna)
		return result
	except subprocess.CalledProcessError as exc:
		result = ""
		return result

def grep_strand(filename, arg):
	cmd = "grep "+arg+" "+filename
	result = subprocess.check_output(cmd, shell=True).rstrip()
	result = result.split("\t")[2]
	return result

my_OUTPUT = []
my_CountBadTranscript = 0
for my_record in open("frameshift.coordinate.txt"):
	my_record = my_record.rstrip()
	my_ALT = my_record.split(" ")[1].split("_")[3]
	my_REF = my_record.split(" ")[1].split("_")[2]
	my_MTPOS = int(my_record.split(" ")[6].split("/")[0].split("-")[0])
	my_MTEND = 0
	try:
		my_MTEND = int(my_record.split(" ")[6].split("/")[0].split("-")[1])
	except IndexError:
		my_MTEND = 0
	my_MTAA = int(my_record.split(" ")[7].split("/")[0].split("-")[0])
	my_PROTLEN = int(my_record.split(" ")[7].split("/")[1])
	my_SEQID = my_record.split(" ")[8]
	my_GENE = my_record.split(" ")[9].rstrip()
	my_SEQID = my_SEQID + "-" + my_GENE
	my_CDSLEN = int(my_record.split(" ")[6].split("/")[1])
	my_ENSTID = my_record.split(" ")[5]
	my_CDS = grep_seq("frameshift.coding.fasta", my_ENSTID)
	if my_CDS.startswith("Sequence"):
		continue
	if my_CDS == "":
		continue
	if not my_CDS.startswith("ATG"):
		my_CountBadTranscript += 1
	my_UTR = grep_seq("frameshift.3UTR.fasta", my_ENSTID)
	if my_UTR.startswith("Sequence"):
		my_UTR=""
	my_TAIL = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	my_STRAND = int(grep_strand("frameshift.strand.txt", my_ENSTID))
	my_MTLEN = my_MTEND - my_MTPOS
	my_MTPOS = my_CDSLEN - my_MTPOS					# position of mutation from the end of CDS
	my_MTPOS = len(my_CDS) - my_MTPOS				# real position of mutation from the begining of CDS
	my_MTAA = my_PROTLEN - my_MTAA
	my_MTAA = int(len(my_CDS) / 3) - my_MTAA

        if len(my_ALT) > len(my_REF):
                if my_STRAND > 0:
                        mutable_seq = my_CDS.tomutable()
                        my_ALT = list(my_ALT)
                        for i in my_ALT[1:]:
                                mutable_seq.insert(my_MTPOS, i)
                                my_MTPOS = my_MTPOS + 1 
                        my_IN_CDS = mutable_seq.toseq()
			my_IN_CDNA = my_IN_CDS + my_UTR + my_TAIL
			my_IN_CDNA = my_IN_CDNA.rstrip()
			AddPeptideRecordToOutput(my_MTAA, my_IN_CDNA, my_SEQID, my_OUTPUT)
                else:
			mutable_seq = my_CDS.tomutable()
                        mutable_seq.reverse_complement()
                        my_ALT = list(my_ALT)
                        for i in my_ALT[1:]:
                                mutable_seq.insert(len(my_CDS) - my_MTPOS, i)
                                my_MTPOS = my_MTPOS - 1 
                        mutable_seq.reverse_complement()
                        my_IN_CDS = mutable_seq.toseq()
			my_IN_CDNA = my_IN_CDS + my_UTR + my_TAIL
			my_IN_CDNA = my_IN_CDNA.rstrip()
			AddPeptideRecordToOutput(my_MTAA, my_IN_CDNA, my_SEQID, my_OUTPUT)
        else:
                mutable_seq = my_CDS.tomutable()
                if my_MTLEN > len(my_REF) - 1:
                        print("deletion goes into intron:", my_record)
                if my_MTLEN > 0:
                        mutable_seq[my_MTPOS - 1 : my_MTPOS + my_MTLEN] = 'M' 
                else:
                        mutable_seq[my_MTPOS - 1] = 'M' 
                mutable_seq.remove('M')
                my_DEL_CDS = mutable_seq.toseq()
		my_DEL_CDNA = my_DEL_CDS + my_UTR + my_TAIL
		my_DEL_CDNA = my_DEL_CDNA.rstrip()
		AddPeptideRecordToOutput(my_MTAA, my_DEL_CDNA, my_SEQID, my_OUTPUT)
SeqIO.write(my_OUTPUT, "frameshift.peptide.fasta", "fasta")    
print("Bad transcripts: ", my_CountBadTranscript)
