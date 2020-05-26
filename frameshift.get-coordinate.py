from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import sys
import subprocess
import os
import numpy as np

for my_record in open("frameshift.genelist"):
        my_record = my_record.rstrip()
	my_ID = my_record.split("\t")[0]
	my_GENE = my_record.split("\t")[1]
	my_COUNTPASS = 0
	my_DATA = ""
	my_SOURCE = sys.argv[1]
	cmd='grep "variant" '+my_SOURCE+'/'+my_ID+'.vcf | grep "|frameshift_variant|" | grep "|'+my_GENE+'|" | sed "s/|frameshift_variant|/ /g"'
	try:
		my_DATA = subprocess.check_output(cmd, shell=True)
	except subprocess.CalledProcessError as exc:
		my_DATA=""
	if my_DATA != "":
		for my_line in my_DATA.splitlines():
			MUTATION=[]
			GENE=[]
			ENSGID=[]
			ENSTID=[]
			CDSPOS=[]
			GENEPOS=[]
			GENELEN=[]
			PROTPOS=[]
			my_COUNTPASS = 0
			if "PASS" in my_line:
				my_COUNTPASS = 1
			for annotation in my_line.split(" ")[1:]:
				MUTATION.append(my_line.split("\t")[0])
				GENE.append(annotation.split("|")[1])
				ENSGID.append(annotation.split("|")[2])
				ENSTID.append(annotation.split("|")[4])
				GENEPOS.append(annotation.split("|")[10])
				GENELEN.append(int(annotation.split("|")[10].split("/")[1]))
				CDSPOS.append(annotation.split("|")[11])
				PROTPOS.append(annotation.split("|")[12])
			cmd="cat frameshift.coordinate.txt | wc -l"
			my_COUNTID = int(subprocess.check_output(cmd, shell=True).rstrip())
			index_max = GENELEN.index(max(GENELEN))
			my_MUTATION=MUTATION[index_max]
			my_GENE=GENE[index_max]
			my_ENSGID=ENSGID[index_max]
			my_ENSTID=ENSTID[index_max]
			my_CDSPOS=CDSPOS[index_max]
			my_PROTPOS=PROTPOS[index_max]
			cmd='echo "'+my_ID+' '+my_MUTATION+' '+str(my_COUNTPASS)+' '+str(my_COUNTPASS)+' '+my_ENSGID+' '+my_ENSTID+' '+my_CDSPOS+' '+my_PROTPOS+' seqID'+str(my_COUNTID)+' '+my_GENE+'" >> frameshift.coordinate.txt'
			return_code = subprocess.call(cmd, shell=True) 



