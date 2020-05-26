from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import sys
import subprocess
import os
from collections import defaultdict

def GrepCoordinateRecord(record, filename):
	ID = record.id.split("-")[0]
	GENE = record.id.split("-")[1]
	cmd='grep "'+ID+' '+GENE+'" '+filename
	try:
		result = subprocess.check_output(cmd, shell=True).rstrip()
		return result
	except subprocess.CalledProcessError:
		result = ""
		return result

def GrepHlaRecord(tcga_id, filename):
	cmd = 'grep "'+tcga_id+'" '+filename
	try:
		result = subprocess.check_output(cmd, shell=True).rstrip()
		return result
	except subprocess.CalledProcessError:
		result = ""
		return result

os.system('mkdir -p netMHCv4')
os.system('mkdir -p netMHCPANv4')
filelist = subprocess.check_output('ls frameshift-peptides-per-patient/', shell=True)
output = open("frameshift.neoantigen.txt","w")

for file_record in filelist.split():
	my_TCGAID = file_record.split(".")[0]
	my_hla_record = GrepHlaRecord(my_TCGAID, "hla.msi-h.tcga.txt")
	if my_hla_record == "":							# iterate if can not find HLA information
		continue
	my_hla_string = my_hla_record.split("\t")[1:7]
	my_hla_string = ' '.join(my_hla_string)
	my_HLAv4 = my_hla_string.replace('*','').replace(':','')
	my_HLAPANv4 = my_hla_string.replace('*','')
	my_list_HLAv4 = subprocess.check_output('cat frequent.hla.txt', shell=True)
	my_list_HLAPANv4 = subprocess.check_output('cat frequent.hla.txt', shell=True)
	my_list_HLAv4 = my_list_HLAv4.replace('*','').replace(':','')
	my_list_HLAPANv4 = my_list_HLAPANv4.replace('*','')
	file_record = "frameshift-peptides-per-patient/"+file_record

	# Run NetMHCv4.0 with Kd=500 threshold
	os.system('mkdir -p netMHCv4_all/'+my_TCGAID)
	for allele in my_list_HLAv4.split('\n'):
		os.system('netMHCv4 -a "'+allele+'" -f "'+file_record+'" -l 9 -s 1 -xls 1 -xlsfile "netMHCv4_all/'+my_TCGAID+'/'+allele+'.neoantigens.xls"')
	my_excellist = subprocess.check_output('ls netMHCv4_all/'+my_TCGAID+'/', shell=True)
	for excel_file in my_excellist.split():
		my_input = open('netMHCv4_all/'+my_TCGAID+'/'+excel_file,"r")
		my_allele = excel_file.split(".")[0]
		for line in my_input:
			if len(line.split()) > 3:
				try:
					MKD = line.split()[3]					# $(echo $line | awk '{print $4}')
					MSCORE = line.split()[6]
					MID = line.split()[2]					# $(echo $line | awk '{print $3}')
					MPEPTIDE = line.split()[1]				# $(echo $line | awk '{print $2}')
					if float(MKD) <= 500:
						string = my_TCGAID+" "+MID+" netMHC IC50 "+my_allele+" "+MPEPTIDE+" "+MKD+"\n"
						output.write(string)
					if float(MSCORE) <= 2:
						string = my_TCGAID+" "+MID+" netMHC SCORE "+my_allele+" "+MPEPTIDE+" "+MSCORE+"\n"
						output.write(string) 
				except ValueError:
					continue
		my_input.close()

	# Run NetMHCPANv4.0 with Kd500 threshold
	os.system('mkdir -p netMHCPANv4_all/'+my_TCGAID+'')
	for allele in my_list_HLAPANv4.split():
		os.system('netMHCPANv4 -a "'+allele+'" -f "'+file_record+'" -l 9 -s 1 -BA -xls 1 -xlsfile "netMHCPANv4_all/'+my_TCGAID+'/'+allele+'.neoantigens.xls"')
	my_excellist = subprocess.check_output('ls netMHCPANv4_all/'+my_TCGAID+'/', shell=True)
	for excel_file in my_excellist.split():
		my_input = open('netMHCPANv4_all/'+my_TCGAID+'/'+excel_file,"r")
		my_allele = excel_file.split(".")[0]
		for line in my_input:
			if len(line.split()) > 3:
				try:
					MKD = line.split()[7]				# $(echo $line | awk '{print $7}')
					MSCORE = line.split()[8]
					MID = line.split()[2]					# $(echo $line | awk '{print $3}')
					MPEPTIDE = line.split()[1]                              # $(echo $line | awk '{print $2}')
					if float(MKD) <= 500:
						string = my_TCGAID+" "+MID+" netMHCpan IC50 "+my_allele+" "+MPEPTIDE+" "+MKD+"\n"
						output.write(string)
					if float(MSCORE) <= 2:
						string = my_TCGAID+" "+MID+" netMHCpan SCORE "+my_allele+" "+MPEPTIDE+" "+MSCORE+"\n"
						output.write(string)
				except ValueError:
					continue
		my_input.close()
output.close()




