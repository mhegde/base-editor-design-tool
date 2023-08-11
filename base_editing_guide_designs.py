import pandas as pd
import csv, argparse, re
import requests, sys, os
from datetime import datetime
from Bio import SeqIO


def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--input-file',
		type=str,
		help='File with Ensembl transcript IDs or fasta file with nucleotide sequences')
	parser.add_argument('--variant-file',
		type=str,
		default='variant_summary.txt',
		help='File with ClinVar SNPs (variant_summary.txt)')
	parser.add_argument('--input-type',
		type=str,
		help='tid for transcript IDs and nuc for nucleotide sequence')
	parser.add_argument('--be-type',
		type=str,
		default='',
		help='Type of base editor')
	parser.add_argument('--pam',
		type=str,
		default='NGG',
		help='PAM sequence for guide design')
	parser.add_argument('--edit-window',
		type=str,
		default='4-8',
		help='Editing window')
	parser.add_argument('--sg-len',
		type=int,
		default=20,
		help='Length of sgRNA')
	parser.add_argument('--edit',
		type=str,
		default='all',
		help='Edit')
	parser.add_argument('--intron-buffer',
		type=int,
		default=30,
		help='How far to tile into introns (bp)')
	parser.add_argument('--filter-gc',
		type=str,
		choices=['True','False'],
		default='False',
		help='Whether to filter out edits in a GC motif')
	parser.add_argument('--output-name',
		type=str,
		help='Output name')
	return parser


def revcom(s):
	basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N','K':'M','M':'K','R':'Y','Y':'R','S':'S','W':'W','B':'V','V':'B','H':'D','D':'H'}
	letters = list(s[::-1])
	letters = [basecomp[base] for base in letters]
	return ''.join(letters)


def get_aa_map():
	aa_map = {'Phe': 'F', 'Leu': 'L', 'Ile': 'I', 'Met': 'M', 'Val': 'V', 'Ser': 'S', 'Pro': 'P', 'Thr': 'T',
			  'Ala': 'A', 'Tyr': 'Y', 'Ter': 'Ter', 'His': 'H', 'Gln': 'Q', 'Asn': 'N', 'Lys': 'K', 'Asp': 'D',
			  'Glu': 'E', 'Cys': 'C', 'Trp': 'W', 'Arg': 'R', 'Gly': 'G'}
	return aa_map


def get_codon_map():
	codon_map = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I',
				 'ATA':'I', 'ATG':'M', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
				 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A',
				 'GCA':'A', 'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'Ter', 'TAG':'Ter', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
				 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'TGT':'C', 'TGC':'C',
				 'TGA':'Ter', 'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
				 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
	return codon_map


def get_pam_window_len(be):
	be_types = {'BE1':'NGG_4-8_20_C-T', 'BE2':'NGG_4-8_20_C-T', 'BE3':'NGG_4-8_20_C-T', 'HF-BE3':'NGG_4-8_20_C-T', 'BE4':'NGG_4-8_20_C-T', 'BE4max':'NGG_4-8_20_C-T',
			   'BE4-Gam': 'NGG_4-8_20_C-T', 'YE1-BE3':'NGG_4-7_20_C-T', 'EE-BE3':'NGG_5-6_20_C-T', 'YE2-BE3':'NGG_5-6_20_C-T', 'YEE-BE3':'NGG_5-6_20_C-T',
			   'VQR-BE3': 'NGAN_4-11_20_C-T', 'VRER-BE3': 'NGCG_3-10_20_C-T', 'SaBE3': 'NNGRRT_3-12_21_C-T', 'SaBE4': 'NNGRRT_3-12_21_C-T',
			   'SaBE4-Gam': 'NNGRRT_3-12_21_C-T', 'Sa(KKH)-BE3': 'NNNRRT_3-12_21_C-T', 'Target-AID': 'NGG_2-4_20_C-T', 'Target-AID-NG': 'NG_2-4_20_C-T',
			   'xBE3': 'NG_4-8_20_C-T', 'eA3A-BE3': 'NG_4-8_20_C-T', 'A3A-BE3': 'NG_4-8_20_C-T', 'BE-PLUS':'NGG_4-14_20_C-T', 'ABE7.9':'NGG_5-8_20_A-G',
			   'ABE7.10': 'NGG_4-7_20_A-G', 'xABE':'NG_4-7_20_A-G', 'ABESa':'NNGRRT_6-12_21_A-G', 'VQR-ABE':'NGA_4-6_20_A-G', 'VRER-ABE':'NGCG_4-6_20_A-G',
			   'Sa(KKH)-ABE':'NNNRRT_6-12_21_A-G'}
	if be in be_types.keys():
		pam, window, sg_len, edit = be_types[be].split('_')
	else:
		print('Please enter ONE of the following: BE1,BE2,BE3,HF-BE3,BE4,BE4max,BE4-Gam,YE1-BE3,EE-BE3,YE2-BE3,'
				 'YEE-BE3,VQR-BE3,VRER-BE3,SaBE3,SaBE4,SaBE4-Gam,Sa(KKH)-BE3,Target-AID,Target-AID-NG,xBE3,eA3A-BE3,'
				 'A3A-BE3,BE-PLUS,ABE7.9,ABE7.10,xABE,ABESa,VQR-ABE,VRER-ABE,Sa(KKH)-ABE')
		sys.exit()
	return pam, window, int(sg_len), edit


def get_pam_pattern(pam):
	code = {'N':'ACTG', 'R':'AG', 'Y':'CT', 'S':'GC', 'W':'AT', 'K':'GT', 'M':'AC', 'B':'CGT', 'D':'AGT', 'H':'ACT', 'V':'ACG'}
	pattern = ''
	for p in pam:
		if p in code.keys():
			pattern = pattern +'['+ code[p] + ']'
		else:
			pattern+=p
	return pattern


def check_ressite_4t(sg):
	res_flag, t4_flag = '', ''
	if ('CGTCTC' in sg or 'GAGACG' in sg or sg.startswith('TCTC') or sg.startswith('AGACG') or sg.endswith('GAGAC')):
		res_flag = 'yes'
	if 'TTTT' in sg:
		t4_flag = 'yes'
	return res_flag, t4_flag


'''
Returns information about gene, assembly, chromosome of specified Ensembl transcript; Also returns absolute values for gene with respect to
genomic locations; Flags genes for absence of UTRs;
'''
def get_tr_info(tr, input_type):
	server = "https://rest.ensembl.org"
	ext = "/lookup/id/"+tr+"?expand=1"
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	if not r.ok:
		r.raise_for_status()
		sys.exit("Transcript not found")
	tr_info = r.json()
	gene_name = tr_info['display_name'].rsplit('-', 1)[0]
	assembly = tr_info['assembly_name']
	gene_strand = tr_info['strand']
	chromosome = tr_info['seq_region_name']
	gene_id = tr_info['Parent']
	if gene_strand == 1:
		gene_start = tr_info['start']
		gene_end = tr_info['end']
		length = gene_end - gene_start
	else:
		gene_start = tr_info['end']
		gene_end = tr_info['start']
		length = gene_start - gene_end
	abs_pos_map, fs = get_absolute_pos(gene_start,gene_end,gene_strand, input_type)
	exons, cds_map = get_exons(tr,length,gene_start,gene_strand)
	if exons != '':
		utr, cds_start_exon, utr5_flag, utr3_flag = get_utrs(tr_info,exons,gene_strand)
	else:
		utr = ''
		cds_start_exon = ''
		utr5_flag = ''
		utr3_flag = ''
	return gene_name, assembly, gene_strand, chromosome, gene_id, exons, cds_map, abs_pos_map, fs, utr, cds_start_exon, utr5_flag, utr3_flag

'''
Returns exon boundaries for transcript;Also returns absolute values for gene with respect to
genomic locations
'''
def get_exons(tr, length, gene_start, gene_strand):
	server = "https://rest.ensembl.org"
	ext = "/map/cds/"+tr+"/1.."+str(length)+"?"
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	if not r.ok:
		#sys.exit()
		return '',''
	decoded = r.json()
	exons_all = decoded['mappings']
	if gene_strand == 1:
		gene_start_pos = exons_all[0]['start'] - gene_start
	else:
		gene_start_pos = gene_start - exons_all[0]['end']
	exons = []
	for i,e in enumerate(exons_all):
		exons.append(str(e['start'])+':'+str(e['end']))
	#abs_pos_map = get_absolute_pos(exons,gene_start_pos,gene_strand)
	cds_map = get_cds_map(exons, gene_strand)
	return exons, cds_map

'''
Generates a hash of utr locations
'''
def get_utr_map(utr, start, end, strand,j):
	k=start
	if strand == 1:
		while k<=end:
			utr[k] = j
			k+=1
			j+=1
	elif strand == -1:
		while k>=end:
			utr[k] = j
			k-=1
			j+=1
	return utr,j

'''
Returns hash of coordinates of 5' and 3' UTR
'''
def get_utrs(tr_info, exons, gene_strand):
	utr_exons = tr_info['Exon']
	ce_1 = exons[0]
	ce_last = exons[-1]
	if gene_strand == 1:
		ce_start,ce_end= ce_1.split(':')
		ce_last_start, ce_last_end = ce_last.split(':')
	elif gene_strand == -1:
		ce_end,ce_start= ce_1.split(':')
		ce_last_end,ce_last_start = ce_last.split(':')
	utr = {}
	j=0
	utr5_flag = 0
	utr3_flag = 0
	for i,e in enumerate(utr_exons):
		if gene_strand == 1:
			if (e['start']<int(ce_start)):
				if (e['end']<int(ce_start)):
					utr,j = get_utr_map(utr,e['start'],e['end'],gene_strand,j)
				elif (e['end']>int(ce_start)):
					cds_start_exon = i+1
					utr,j = get_utr_map(utr,e['start'],int(ce_start)-1,gene_strand,j)
			if (e['end'] > int(ce_last_end)):
				if (e['start']>int(ce_last_end)):
					utr,j = get_utr_map(utr,e['start'],e['end'],gene_strand,j)
				elif (e['start']<int(ce_last_end)):
					utr,j = get_utr_map(utr,int(ce_last_end)+1,e['end'],gene_strand,j)
			if (e['start'] == int(ce_start)):
				if i == 0:
					utr5_flag = 1
					cds_start_exon = i+1
				else:
					cds_start_exon = i+1
			if (e['end'] == int(ce_last_end)):
				if i == len(utr_exons)-1:
					utr3_flag = 1
		elif gene_strand == -1:
			if (e['end']>int(ce_start)):
				if (e['start']>int(ce_start)):
					utr,j = get_utr_map(utr,e['end'],e['start'],gene_strand,j)
				elif (e['start']<int(ce_start)):
					cds_start_exon = i + 1
					utr,j = get_utr_map(utr,e['end'],int(ce_start)+1,gene_strand,j)
			if (e['start']<int(ce_last_end)):
				if (e['end']>int(ce_last_end)):
					utr,j = get_utr_map(utr,int(ce_last_end)-1,e['start'],gene_strand,j)
				elif (e['end']<int(ce_last_end)):
					utr,j = get_utr_map(utr,e['end'],e['start'],gene_strand,j)
			if (e['end'] == int(ce_start)):
				if i == 0:
					utr5_flag = 1
					cds_start_exon = i+1
				else:
					cds_start_exon = i+1
			if (e['start'] == int(ce_last_end)):
				if i == len(utr_exons) - 1:
					utr3_flag = 1
	return utr, cds_start_exon, utr5_flag, utr3_flag

'''
Generates absolute positions for genomic locations of transcript starting from 40 nucleotides flanking sequence preceding
the 5'UTR up to 40 nucleotides flanking sequence succeeding the end of 3'UTR 
'''
def get_absolute_pos(gene_start, gene_end, gene_strand, input_type):
	abs_pos_map = {}
	fs = {}                 #Genomic location of flanking sequences
	#i=gene_start
	j=0
	if gene_strand == 1:
		if input_type == 'tid':
			i = gene_start-40   #including 40 nucs of flanking sequence annotation to abs_pos_map
			end = gene_end+40
		else:
			i = gene_start
			end = gene_end
		while i<=end:
			abs_pos_map[i] = j
			if i < gene_start or i > gene_end:
				fs[i] = j
			i+=1
			j+=1

	else:
		if input_type == 'tid':
			i = gene_start+40   #including 40 nucs of flanking sequence annotation to abs_pos_map
			end = gene_end-40
		else:
			i = gene_start
			end = gene_end
		while i>=end:
			abs_pos_map[i] = j
			if i > gene_start or i < gene_end:
				fs[i] = j
			i-=1
			j+=1

	return abs_pos_map, fs

'''
Generates absolute cds positions for genomic locations of transcript
'''
def get_cds_map(exons, gene_strand):
	cds_map={}
	j=0
	for e in exons:
		if gene_strand == 1:
			e = e.split(':')
			i = int(e[0])
			while(i <= int(e[1])):
				cds_map[i] = j
				i=i+1
				j=j+1
		else:
			e = e.split(':')
			i = int(e[1])
			while(i >= int(e[0])):
				cds_map[i] = j
				i-=1
				j+=1
	#j=j-1
	return cds_map


'''
Returns sequence of transcript along with 5' and 3' UTR
'''
def get_tr_sequence(tr):
	server = "https://rest.ensembl.org"
	ext = "/sequence/id/"+tr+"?content-type=text/plain;expand_5prime=40;expand_3prime=40"
	r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
	if not r.ok:
		r.raise_for_status()
		return ''
	sequence = r.text
	return sequence

'''
Returns protein sequence of transcript
'''
def get_pro_sequence(tr):
	server = "https://rest.ensembl.org"
	ext = "/sequence/id/"+tr+"?content-type=text/plain;type=protein"
	r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
	if not r.ok:
		r.raise_for_status()
		return ''
	pro_sequence = r.text
	return pro_sequence

'''
Returns CDS sequence of transcript
'''
def get_cds_sequence(tr):
	server = "https://rest.ensembl.org"
	ext = "/sequence/id/"+tr+"?content-type=text/plain;type=cds"
	r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
	if not r.ok:
		r.raise_for_status()
		return ''
	cds_sequence = r.text
	return cds_sequence

'''
Translates sgRNA sequence and annotates frame
'''
def get_sgrna_translated_seq(sgrna, cds_map, abs_pos, fs, sgrna_start_pos, gene_strand, sgrna_strand, utr, e, label):
	map_key = abs_pos.keys()[abs_pos.values().index(sgrna_start_pos)]
	sgrna_trans = {}
	if sgrna_strand == 'sense':
		for i,n in enumerate(sgrna):
			if gene_strand == 1:
				if map_key not in cds_map.keys():
					if map_key in utr.keys():
						sgrna_trans[n+str(i+1)] = '0_U'
					elif map_key in fs.keys():
						sgrna_trans[n + str(i + 1)] = '0_FS'
					else:
						if map_key < int(e[0]):
							sgrna_trans[n+str(i+1)] = '0_'+label+':-'+str(int(e[0]) - map_key)
						elif map_key > int(e[1]):
							sgrna_trans[n+str(i+1)] = '0_'+label+':+'+str(map_key - int(e[1]))
					map_key +=1
				else:
					sgrna_trans[n+str(i+1)] = str((cds_map[map_key]/3)+1)+'_'+str(cds_map[map_key]%3)
					map_key +=1
			else:
				if map_key not in cds_map.keys():
					if map_key in utr.keys():
						sgrna_trans[n+str(i+1)] = '0_U'
					elif map_key in fs.keys():
						sgrna_trans[n + str(i + 1)] = '0_FS'
					else:
						if map_key < int(e[0]):
							sgrna_trans[n+str(i+1)] = '0_'+label+':+'+str(int(e[0]) - map_key)
						elif map_key > int(e[1]):
							sgrna_trans[n+str(i+1)] = '0_'+label+':-'+str(map_key - int(e[1]))
					map_key-=1
				else:
					sgrna_trans[n+str(i+1)] = str((cds_map[map_key]/3)+1)+'_'+str(cds_map[map_key]%3)
					map_key -=1
	elif sgrna_strand == 'antisense':
		for i,n in enumerate(sgrna):
			if gene_strand == 1:
				if map_key not in cds_map.keys():
					if map_key in utr.keys():
						sgrna_trans[n+str(len(sgrna)-i)] = '0_U'
					elif map_key in fs.keys():
						sgrna_trans[n + str(len(sgrna) - i)] = '0_FS'
					else:
						if map_key < int(e[0]):
							sgrna_trans[n+str(len(sgrna)-i)] = '0_'+label+':-'+str(int(e[0]) - map_key)
						elif map_key > int(e[1]):
							sgrna_trans[n+str(len(sgrna)-i)] = '0_'+label+':+'+str(map_key - int(e[1]))
					map_key+=1
				else:
					sgrna_trans[n+str(len(sgrna)-i)] = str((cds_map[map_key]/3)+1)+'_'+str(cds_map[map_key]%3)
					map_key+=1
			else:
				if map_key not in cds_map.keys():
					if map_key in utr.keys():
						sgrna_trans[n+str(len(sgrna)-i)] = '0_U'
					elif map_key in fs.keys():
						sgrna_trans[n + str(len(sgrna) - i)] = '0_FS'
					else:
						if map_key < int(e[0]):
							sgrna_trans[n+str(len(sgrna)-i)] = '0_'+label+':+'+str(int(e[0]) - map_key)
						elif map_key > int(e[1]):
							sgrna_trans[n+str(len(sgrna)-i)] = '0_'+label+':-'+str(map_key - int(e[1]))
					map_key-=1
				else:
					sgrna_trans[n+str(len(sgrna)-i)] = str((cds_map[map_key]/3)+1)+'_'+str(cds_map[map_key]%3)
					map_key -=1
	return sgrna_trans


'''
Cleans variant file from ClinVar database
'''
def parse_variant_df(variant_df):
	# Remove non-GRCh38 rows, non-SNPs, and chromosomes other than 1-22, X, and Y.
	temp_variant_df = variant_df[(variant_df.Assembly == 'GRCh38')
								& (variant_df.Type == 'single nucleotide variant')
								& (variant_df.Chromosome != 'MT')
								& (variant_df.Chromosome != 'na')
								& (variant_df.ReferenceAlleleVCF != 'na')]
	parsed_variant_df = temp_variant_df.copy()
	parsed_variant_df.index = range(0,len(parsed_variant_df))
	parsed_variant_df.rename(columns = {'Start':'ClinVar_SNP_Position'}, inplace=True)
	parsed_variant_df.loc[:,'RefSeqID'] = parsed_variant_df.loc[:,'Name'].str.split('(',1).str[0]
	parsed_variant_df = parsed_variant_df[['#AlleleID',
										   'RefSeqID',
										   'Name',
										   'GeneSymbol',
										   'ClinicalSignificance',
										   'PhenotypeList',
										   'ClinVar_SNP_Position',
										   'ReferenceAlleleVCF',
										   'AlternateAlleleVCF',
										   'ReviewStatus']]
	return parsed_variant_df

'''
Parses the ClinVar SNP name
'''
def parse_snp_name(snp_name, aa_map):
	snp_aa = ''
	snp_aa_to = ''
	snp_aa_from = ''
	snp_aa_num = ''
	# If the SNP is in an exon, pull amino acid change out of ClinVar SNP name
	if snp_name.find('(p.') != -1:
		snp_aa = snp_name[snp_name.find('(p.')+3:-1]
		three_letter_aa = re.split(r'(\D+)', snp_aa)
		snp_aa_from = three_letter_aa[1]
		snp_aa_to = three_letter_aa[3]
		snp_aa_num = three_letter_aa[2]
		# Silent mutations are annotated as '='
		if snp_aa_to == '=':
			snp_aa_to = snp_aa_from
		snp_aa = snp_aa_from + snp_aa_num + snp_aa_to
	return snp_aa, snp_aa_from, snp_aa_num, snp_aa_to


'''
Returns:
	codon_pos_list, a list of 3 genomic positions for each nucleotide in a codon
	edit_gen_pos_list, a list of the genomic positions for all edited nucleotides in the codon
'''
def get_genomic_pos_list(edit_indices, gene_strand, sgrna_strand, sg_gen_pos):
	edit_gen_pos_list = []
	codon_pos_list = []
	for edit_pos,frame in edit_indices.iteritems():
		# if sgRNA is in + strand
		if ((gene_strand == 1) and (sgrna_strand == 'sense')) or ((gene_strand == -1) and (sgrna_strand == 'antisense')):
			edit_gen_pos = (sg_gen_pos+int(edit_pos)-1)
		# elif sgRNA is in - strand
		elif ((gene_strand == 1) and (sgrna_strand == 'antisense')) or ((gene_strand == -1) and (sgrna_strand == 'sense')):
			edit_gen_pos = (sg_gen_pos-int(edit_pos)+1)
		edit_gen_pos_list.append(str(edit_gen_pos))
		# For edits in exon, add the other two positions in the codon to the codon_pos_list
		# Therefore we can annotate any SNPs that affect that codon
		if (frame != 'U') and ('Exon' not in frame) and (frame!='FS'):
			frame = int(frame)
			if gene_strand == 1:
				codon_pos_list = [edit_gen_pos-frame,edit_gen_pos-frame+1,edit_gen_pos-frame+2]
			elif gene_strand == -1:
				codon_pos_list = [edit_gen_pos+frame,edit_gen_pos+frame-1,edit_gen_pos+frame-2]
		elif (frame == 'U') or ('Exon' in frame) or (frame == 'FS'):
			# For utr and introns, only check for snps at the location of the edit
			codon_pos_list.extend([edit_gen_pos])
	return codon_pos_list, edit_gen_pos_list

'''
Returns dataframe containing information about the pathogenicity and position of SNPs created by edit
'''
def get_snps(edit_map, edit, sg_gen_pos, gene_strand, sgrna_strand, gene_variant_df, aa_map):
	snp_type_list = []
	snp_info = []
	edit_nuc, edit_to = edit.split('-')
	all_snps = gene_variant_df['ClinVar_SNP_Position'].tolist()
	# Iterate through all amino acid changes
	for k,v in edit_map.iteritems():
		temp_snp_type_list = []
		edit_indices = {}
		edits = k.split('_')
		sgrna_edit_nuc = edits[0]
		for i in edits[1:]:
			edit_index, frame = i.split('-',1)
			edit_indices[edit_index] = frame
		edit_cat = v.split('_')[1]
		# UTR or flanking sequence mutations
		if 'UTR' in v or 'Flanking' in v:
			aa_edit,edit_type = v.split('_')
			old_codon, new_codon = '', ''
			aa_from = ''
			aa_num = ''
			aa_to = ''
		# Coding mutations
		elif '_' in v and 'Exon' not in v:
			aa_edit, edit_type, old_codon, new_codon = v.split('_')
			aa_from = aa_edit[0:3]
			aa_num = aa_edit[3:-3]
			aa_to = aa_edit[-3:]
		# Intronic mutations
		elif 'Exon' in v:
			aa_edit = 'intron'
			edit_type, old_codon, new_codon = '', '', ''
			aa_from = ''
			aa_num = ''
			aa_to = ''
		codon_pos_list, edit_gen_pos_list = get_genomic_pos_list(edit_indices, gene_strand, sgrna_strand, sg_gen_pos)
		# Check if there is any overlap between codon_pos_list and all_snps
		if any(i in codon_pos_list for i in all_snps):
			clinvar_snps_df = gene_variant_df[gene_variant_df.ClinVar_SNP_Position.isin(codon_pos_list)].loc[:,['Name','ClinicalSignificance','ClinVar_SNP_Position','ReferenceAlleleVCF','AlternateAlleleVCF','ReviewStatus']]
			for index,row in clinvar_snps_df.iterrows():
				snp_aa, snp_aa_from, snp_aa_num, snp_aa_to = parse_snp_name(row.Name, aa_map)
				# First check for nucleotide position
				if str(row.ClinVar_SNP_Position) not in edit_gen_pos_list:
					same_nucleotide_pos = False
					same_nucleotide_change = False
				elif str(row.ClinVar_SNP_Position) in edit_gen_pos_list:
					same_nucleotide_pos = True
					# If sgRNA is in the forward strand, C>T or A>G SNPs will be created
					if ((gene_strand == 1) and (sgrna_strand == 'sense') or ((gene_strand == -1) and sgrna_strand == 'antisense')):
						if row.AlternateAlleleVCF == edit_to:
							same_nucleotide_change = True
						else:
							same_nucleotide_change = False
					# If sgRNA is in the reverse strand, G>A or T>C SNPs will be created
					elif ((gene_strand == 1) and (sgrna_strand == 'antisense') or ((gene_strand == -1) and sgrna_strand == 'sense')):
						if row.AlternateAlleleVCF == revcom(edit_to):
							same_nucleotide_change = True
						else:
							same_nucleotide_change = False
				if snp_aa != '':
					if str(aa_num) != str(snp_aa_num):
						same_aa_pos = False
						if (aa_from == snp_aa_from) and (aa_to == snp_aa_to):
							same_aa_change = True
						else:
							same_aa_change = False
					elif str(aa_num) == str(snp_aa_num):
						same_aa_pos = True
						if (aa_from == snp_aa_from) and (aa_to == snp_aa_to):
							same_aa_change = True
						else:
							same_aa_change = False
				else:
					same_aa_pos = 'N/A'
					same_aa_change = 'N/A'
				# Append the clinical significance of exact match SNPs to snp_type_list
				# Do not require that the aa number matches
				if same_nucleotide_pos and same_nucleotide_change:
					if same_aa_pos and same_aa_change:
						temp_snp_type_list.append(row.ClinicalSignificance)
					elif (same_aa_pos == 'N/A') and (same_aa_change == 'N/A'):
						temp_snp_type_list.append(row.ClinicalSignificance)
				snp_info.append([str(sgrna_edit_nuc+'_'+'_'.join(edit_indices.keys())),
									 ';'.join(edit_gen_pos_list),
									 aa_edit,
									 old_codon,
									 new_codon,
									 edit_cat,
									 snp_aa] +
									 row.tolist() +
									 [same_nucleotide_pos,
									 same_nucleotide_change,
									 same_aa_pos,
									 same_aa_change])
		else:
			snp_info.append([str(sgrna_edit_nuc+'_'+'_'.join(edit_indices.keys())), ';'.join(edit_gen_pos_list),
				aa_edit,
				old_codon,
				new_codon,
				edit_cat])
		if not temp_snp_type_list:
			temp_snp_type_list.append('None')
		snp_type_list.extend(temp_snp_type_list)		
	return snp_type_list, snp_info


'''
Returns clinical significances as ;-separated string
'''
def get_clinical_sig(snp_type_list):
	if not snp_type_list:
		clinical_sig = ''
	else:
		clinical_sig = ';'.join(snp_type_list)
	return clinical_sig



'''Filters out edits that are in a GC motif'''
def filter_gc_motifs(sgrna_context, context_index, sgrna_strand):
	if sgrna_strand == 'sense':
		motif = sgrna_context[context_index-1]
		if motif == 'G':
			return False
	elif sgrna_strand == 'antisense':
		motif = sgrna_context[len(sgrna_context) - context_index-2]
		if motif == 'G':
			return False			
	return True

def filter_gc_motifs_for_aa(sgrna_strand,sgrna_context,codon_start,k):
	if sgrna_strand == 'sense':
		if sgrna_context[codon_start+k-1] == 'G':
			return False
	elif sgrna_strand == 'antisense':
		if sgrna_context[len(sgrna_context) - (codon_start+k)-2] == 'G':
			return False
	return True


'''
Returns edits for sgRNA in specified window
Also returns number of silent edits in window
'''
def get_edits(edit_map, context, window, edit, sgrna_trans, codon_map, j, sgrna_strand, pam, window_start, window_end, aa_map, filter_gc, sgrna_context):
	error = ''
	num_silent = 0
	num_stop = 0
	edit_nuc, edit_to = edit.split('-')
	if sgrna_strand == 'antisense':
		edit_nuc, edit_to = revcom(edit_nuc), revcom(edit_to)
		context_index_track = len(sgrna_trans)-j+len(pam)+3
	for i,n in enumerate(window):
		if n == edit_nuc:
			if sgrna_strand == 'sense':
				nuc_edit_pos = edit_nuc+str(j+1)
				aa_num,frame = sgrna_trans[nuc_edit_pos].split('_')
				context_index = j+4
			elif sgrna_strand == 'antisense':
				nuc_edit_pos = edit_nuc+str(j)
				aa_num,frame = sgrna_trans[nuc_edit_pos].split('_')
				nuc_edit_pos = revcom(edit_nuc)+str(j)
				context_index = context_index_track

			# Check motif
			if filter_gc:
				proceed = filter_gc_motifs(sgrna_context, context_index, sgrna_strand)
			else:
				proceed = True
			if proceed:
				if frame == '0':
					codon_start = context_index
					codon_end = context_index+3
					old_codon = context[codon_start:codon_end]
				elif frame == '1':
					codon_start = context_index-1
					codon_end = context_index+2
					old_codon = context[codon_start:codon_end]
				elif frame == '2':
					codon_start = context_index-2
					codon_end = context_index+1
					old_codon = context[codon_start:codon_end]
				else:
					old_codon = ''
				if old_codon != '':
					if old_codon in codon_map.keys():
						old_aa = codon_map[old_codon]
						old_aa_3 = list(aa_map.keys())[list(aa_map.values()).index(old_aa)]
						new_codon = []
						edit_indices = []
						for k,x in enumerate(old_codon):
							if x == edit_nuc:
								if filter_gc:
									motif_check = filter_gc_motifs_for_aa(sgrna_strand,sgrna_context,codon_start,k)
								else:
									motif_check = True
								if motif_check:								
									if aa_num+'_'+str(k) in sgrna_trans.values():
										nuc_index = sgrna_trans.keys()[sgrna_trans.values().index(aa_num+'_'+str(k))][1:]
										if (int(nuc_index) >= window_start) and (int(nuc_index) <= window_end):
											new_codon.append(edit_to)
											nuc_index = nuc_index + '-' + str(k)
											edit_indices.append(nuc_index)
										else:
											new_codon.append(x)
									else:
										new_codon.append(x)
								else:
									new_codon.append(x)
							else:
								new_codon.append(x)
					else:
						error = 'Codon '+old_codon+'not standard codon'
						return '', '', error, '', '', '', '', ''
					new_codon = ''.join(new_codon)
					new_aa = codon_map[new_codon]
					new_aa_3 = list(aa_map.keys())[list(aa_map.values()).index(new_aa)]
					#aa_edit = old_aa+str(aa_num)+new_aa
					aa_edit = old_aa_3 + str(aa_num) + new_aa_3
					edit_indices.sort(key=lambda x: int(x.split('-')[0]))
					if sgrna_strand == 'antisense':
						nuc_edit_pos = revcom(edit_nuc)+'_'+'_'.join(edit_indices)
					else:
						nuc_edit_pos = edit_nuc+'_'+'_'.join(edit_indices)
					if nuc_edit_pos not in edit_map.keys():
						if old_aa == new_aa:
							edit_cat = 'Silent'
							num_silent+=1
						elif old_aa != 'Ter' and new_aa == 'Ter':
							edit_cat = 'Nonsense'
							num_stop+=1
						else:
							edit_cat = 'Missense'
						edit_map[nuc_edit_pos] = aa_edit+'_'+edit_cat+'_'+old_codon+'_'+new_codon
				else:
					nuc_edit_pos = nuc_edit_pos[0]+'_'+nuc_edit_pos[1:]+'-'+frame
					if nuc_edit_pos not in edit_map.keys():
						if 'Exon' in frame:
							if (frame.split(':')[1] == '-1') or (frame.split(':')[1] == '-2'):
								edit_cat = 'Splice-acceptor'
							elif (frame.split(':')[1] == '+1') or (frame.split(':')[1] == '+2'):
								edit_cat = 'Splice-donor'
							else:
								edit_cat = 'Intron'
							edit_map[nuc_edit_pos] = frame+'_'+edit_cat
						elif frame == 'U':
							edit_map[nuc_edit_pos] = 'utr_UTR'
						elif frame == 'FS':
							edit_map[nuc_edit_pos] = 'flankseq_Flanking'

		if sgrna_strand == 'sense':
			j+=1
		else:
			j-=1
			context_index_track+=1
		transcript_ref_allele = edit_nuc
		transcript_alt_allele = edit_to
		if gene_strand == 1:
			genome_ref_allele = edit_nuc
			genome_alt_allele = edit_to
		elif gene_strand == -1:
			genome_ref_allele = revcom(edit_nuc)
			genome_alt_allele = revcom(edit_to)
	return edit_map, num_silent, error, num_stop, transcript_ref_allele,transcript_alt_allele,genome_ref_allele,genome_alt_allele

'''
Returns edits for sgRNA, also returns number of silent edits
'''
def get_edit_info(context, sgrna, sgrna_strand, edit, window, pam, sgrna_trans, codon_map, sg_gen_pos, gene_strand, gene_variant_df, aa_map, filter_gc, sgrna_context):
	window_start, window_end = window.split('-')
	sgrna_window = sgrna[int(window_start)-1:int(window_end)]
	j_window = int(window_start)-1
	if sgrna_strand == 'antisense':
		sgrna_window = revcom(sgrna_window)
		j_window = int(window_end)

	edit_map = {}
	edit_map, window_silent, cds_error, num_stop, transcript_ref_allele, transcript_alt_allele, genome_ref_allele, genome_alt_allele = get_edits(edit_map,context,sgrna_window,edit,sgrna_trans,codon_map,j_window,sgrna_strand,pam,int(window_start),int(window_end), aa_map,filter_gc, sgrna_context)
	snp_type_list, snp_info = get_snps(edit_map, edit, sg_gen_pos, gene_strand, sgrna_strand, gene_variant_df, aa_map)
	clinical_sig = get_clinical_sig(snp_type_list)
	return edit_map, window_silent, cds_error, num_stop, clinical_sig, snp_info, transcript_ref_allele, transcript_alt_allele, genome_ref_allele, genome_alt_allele

'''
Returns edits in format suitable for writing to file
'''
def get_print_edits(edit_map):
	nuc_edits = ''
	aa_edits = ''
	cat = ''
	old_codon = ''
	new_codon = ''
	num_edits = 0
	for k in sorted(edit_map, key= lambda x: int(x.split('_')[1].split('-')[0])):
		v = edit_map[k]
		if '_' in v:
			vals = v.split('_')
			aa_edits = aa_edits + vals[0] + ';'
			cat = cat+vals[1]+';'
			# len(vals) > 2 for coding sequence, <= 2 for non-coding (intron, UTR, flanking)
			if len(vals) > 2:
				old_codon = old_codon + vals[2] + ';'
				new_codon = new_codon + vals[3] + ';'
		else:
			aa_edits = aa_edits + v + ';'
		ne_edits = k.split('_')
		nuc_edits = nuc_edits+ne_edits[0]
		for ne in ne_edits[1:]:
			nuc_edits = nuc_edits +'_' + ne.split('-')[0]
		nuc_edits = nuc_edits + ';'
		num_edits += 1
	return nuc_edits, aa_edits, old_codon, new_codon, cat, num_edits


def get_context_for_trans(ct_index, ct_index_check, abs_pos, cds_map, fs, sgrna_context):
	context_for_trans = ''
	flag = 0
	i_count = 0
	cds_pos = ''
	while ct_index < ct_index_check:
		gen_pos = abs_pos.keys()[abs_pos.values().index(ct_index)]
		if gen_pos in cds_map.keys():
			if flag == 0: #Check to see if ct_index has encountered CDS
				cds_pos = cds_map[gen_pos]
				flag = 1
			cds_str = cds_sequence[cds_map[gen_pos]:cds_map[gen_pos] + (len(sgrna_context) - len(context_for_trans))]
			context_for_trans += cds_str
			ct_index += len(cds_str)
		elif gen_pos in utr.keys():
			context_for_trans += 'U'
			ct_index += 1
		elif gen_pos in fs.keys():
			context_for_trans += 'F'
			ct_index += 1
		else:
			context_for_trans += 'I'
			if flag == 0:
				i_count += 1
			ct_index += 1

	return context_for_trans, cds_pos, i_count

'''
Designs sgRNAs for specified PAM sequence and writes to output    
'''


def design_sgrnas(gene_name, assembly, chromosome, gene_id,w, gene_seq, abs_pos, fs, cds_map, utr, t, pam, exons, gene_strand, edit, window, cds_sequence, pam_len, sg_len, cds_start_exon, w_error, w_clin, gene_variant_df, aa_map, input_type, intron_buffer,filter_gc):
	current_exon = cds_start_exon
	for i,e in enumerate(exons):
		e = e.split(':')
		label = 'Exon'+str(current_exon)
		if gene_strand == 1:
			start_pos = abs_pos[int(e[0])]
			pos_end = abs_pos[int(e[1])]
		else:
			start_pos = abs_pos[int(e[1])]
			pos_end = abs_pos[int(e[0])]

		if input_type == 'tid':
			pos = start_pos - intron_buffer
			pos_end = pos_end + intron_buffer
		else:
			pos = start_pos+4

		while pos < pos_end:
			sg_str_anti = 0
			sg_str_sense = 0
			target = gene_seq[pos:(pos+pam_len+sg_len)]
			context = gene_seq[(pos-4):(pos+pam_len+sg_len+4)]
			if len(context) == sg_len+pam_len+8:  #Check if full context is available for target sequence
				m = re.search('[^ATCG]', context)
				if m is None:
					start = target[0:pam_len]
					finish = target[sg_len:sg_len + pam_len]
					m_fwd = re.search(get_pam_pattern(pam), finish)
					m_rev = re.search(get_pam_pattern(revcom(pam)), start)
					if m_rev is not None:
						sg_str_anti = 1
					if m_fwd is not None:
						sg_str_sense = 1
			if sg_str_anti == 1:
				error = ''
				sgrna_for_trans = target[pam_len:pam_len+sg_len]
				sgrna = revcom(sgrna_for_trans)
				res_flag, t4_flag = check_ressite_4t(sgrna)
				sgrna_strand = 'antisense'
				sgrna_pam = revcom(target[0:pam_len])
				sgrna_context = revcom(gene_seq[(pos-3):(pos+pam_len+sg_len+4)])
				sgrna_start_pos = pos+pam_len
				sgrna_end_pos = pos+sg_len+pam_len-1

				ct_index = pos-3
				context_for_trans, cds_pos, i_count = get_context_for_trans(ct_index, pos+sg_len+pam_len+4, abs_pos, cds_map, fs, sgrna_context)
				m = re.search('[^F]',context_for_trans)
				if m is None:
					context_for_trans = ''
					error = 'Entire context in flanking sequence'


				if 'I' in context_for_trans:
					map_key_context_start = abs_pos.keys()[abs_pos.values().index(pos - 3)]
					map_key_context_end = abs_pos.keys()[abs_pos.values().index(pos + pam_len + sg_len + 3)]
					if map_key_context_start in cds_map.keys():
						# True if sgRNA is in first coding exon or inner exon, false if in last coding exon
						if (cds_map[map_key_context_start] + pam_len + sg_len + 7) <= len(cds_sequence):
							# Get CDS from following exons
							context_for_trans = cds_sequence[cds_map[map_key_context_start]:cds_map[map_key_context_start] + pam_len + sg_len + 7]
						else:
							# Get all following CDS, then fill in U's
							context_for_trans = cds_sequence[cds_map[map_key_context_start]:len(cds_sequence)]
							context_for_trans = context_for_trans + (sg_len + pam_len + 7 - len(context_for_trans))*'U'
					elif map_key_context_end in cds_map.keys():
						# True if sgRNA is in inner exon or last coding exon, false if in first coding exon
						if (cds_map[map_key_context_end] - sg_len - pam_len - 6) >= 0:
							# Get CDS from previous exons
							context_for_trans = cds_sequence[(cds_map[map_key_context_end] - sg_len - pam_len - 6):cds_map[map_key_context_end]+1]
						else:
							# Get all of preceding CDS, then fill in U's
							context_for_trans = cds_sequence[0:cds_map[map_key_context_end]+1] #+1 because ceiling is not included
							context_for_trans = (sg_len + pam_len + 7 - len(context_for_trans))*'U' + context_for_trans
					else:
						if cds_pos != '':
							context_for_trans = cds_sequence[cds_pos-i_count:cds_pos-i_count+len(sgrna_context)]
						else:
							context_for_trans = sgrna_context

				if context_for_trans == '':
					if error == '':
						error = 'No context_for_trans found'
					print(error)
					w_error.writerow([gene_name, tr, sgrna, sgrna_strand, error])

				if context_for_trans != '':
					sg_gen_pos = abs_pos.keys()[abs_pos.values().index(sgrna_end_pos)]
					sgrna_trans = get_sgrna_translated_seq(sgrna_for_trans, cds_map, abs_pos, fs, sgrna_start_pos, gene_strand, sgrna_strand, utr, e, label)
					edit_map, window_silent, cds_error, num_stop, clinical_sig, snp_info, transcript_ref_allele, transcript_alt_allele, genome_ref_allele, genome_alt_allele = get_edit_info(context_for_trans, sgrna, sgrna_strand, edit, window, pam, sgrna_trans, codon_map, sg_gen_pos, gene_strand, gene_variant_df, aa_map, filter_gc, sgrna_context)
					if cds_error != '':
						w_error.writerow([gene_name, tr, sgrna, sgrna_strand, cds_error])
						return 0
					nuc_edits, aa_edits, old_codon, new_codon, cat, num_edits = get_print_edits(edit_map)
					w.writerow([sgrna, sgrna_context, gene_name, gene_id,t, gene_strand, assembly, transcript_ref_allele, transcript_alt_allele,
								genome_ref_allele, genome_alt_allele, chromosome, sg_gen_pos, sgrna_strand,
								sgrna_pam, edit, num_edits, window_silent, nuc_edits, aa_edits, cat, clinical_sig, res_flag, t4_flag])
					for snp_row in snp_info:
						w_clin.writerow([sgrna, sgrna_strand, sgrna_context, chromosome, gene_name, gene_strand, edit, transcript_ref_allele,
										 transcript_alt_allele, genome_ref_allele, genome_alt_allele] + snp_row)

			if sg_str_sense == 1:
				error = ''
				sgrna = target[0:sg_len]
				res_flag, t4_flag = check_ressite_4t(sgrna)
				sgrna_strand = 'sense'
				sgrna_pam = target[sg_len:sg_len+pam_len]
				sgrna_context = gene_seq[pos-4:pos+sg_len+pam_len+3]
				sgrna_start_pos = pos
				sgrna_end_pos = pos+sg_len-1

				ct_index = pos-4
				context_for_trans, cds_pos, i_count = get_context_for_trans(ct_index, pos+sg_len+pam_len+3, abs_pos, cds_map, fs, sgrna_context)
				m = re.search('[^F]',context_for_trans)
				if m is None:
					context_for_trans = ''
					error = 'Entire context in flanking sequence'


				if 'I' in context_for_trans:
					map_key_context_start = abs_pos.keys()[abs_pos.values().index(pos-4)]
					map_key_context_end = abs_pos.keys()[abs_pos.values().index(pos + pam_len + sg_len + 3)]
					if map_key_context_start in cds_map.keys():
						if (cds_map[map_key_context_start]+sg_len+pam_len+7) <= len(cds_sequence):
							context_for_trans = cds_sequence[cds_map[map_key_context_start]:cds_map[map_key_context_start]+sg_len+pam_len+7]
						else:
							context_for_trans = cds_sequence[cds_map[map_key_context_start]:len(cds_sequence)]
							context_for_trans = context_for_trans + (sg_len + pam_len + 7 - len(context_for_trans))*'U'
					elif map_key_context_end in cds_map.keys():
						if (cds_map[map_key_context_end]-pam_len-sg_len-7) >= 0:
							# Get CDS from previous exons
							context_for_trans = cds_sequence[cds_map[map_key_context_end]-pam_len-sg_len-7:cds_map[map_key_context_end]]
						else:
							context_for_trans = cds_sequence[0:cds_map[map_key_context_end]]
							context_for_trans = (sg_len + pam_len + 7 - len(context_for_trans))*'U' + context_for_trans
					else:
						if cds_pos != '':
							context_for_trans = cds_sequence[cds_pos - i_count:cds_pos - i_count + len(sgrna_context)]
						else:
							context_for_trans = sgrna_context

				if context_for_trans == '':
					if error == '':
						error = 'No context_for_trans found'
					print(error)
					w_error.writerow([gene_name, tr, sgrna, sgrna_strand,error])

				if context_for_trans != '':
					sg_gen_pos = abs_pos.keys()[abs_pos.values().index(sgrna_start_pos)]
					sgrna_trans = get_sgrna_translated_seq(sgrna, cds_map, abs_pos, fs, sgrna_start_pos, gene_strand, sgrna_strand, utr, e, label)
					edit_map, window_silent, cds_error, num_stop, clinical_sig, snp_info, transcript_ref_allele, transcript_alt_allele, genome_ref_allele, genome_alt_allele = get_edit_info(context_for_trans, sgrna, sgrna_strand, edit, window, pam, sgrna_trans, codon_map, sg_gen_pos, gene_strand, gene_variant_df, aa_map, filter_gc, sgrna_context)
					if cds_error != '':
						w_error.writerow([gene_name, tr, sgrna, sgrna_strand, cds_error])
						return 0
					nuc_edits, aa_edits, old_codon, new_codon, cat, num_edits = get_print_edits(edit_map)
					w.writerow([sgrna, sgrna_context, gene_name, gene_id,t, gene_strand, assembly, transcript_ref_allele, transcript_alt_allele,
								genome_ref_allele, genome_alt_allele, chromosome, sg_gen_pos, sgrna_strand,
								sgrna_pam, edit, num_edits, window_silent, nuc_edits, aa_edits, cat, clinical_sig, res_flag, t4_flag])
					for snp_row in snp_info:
						w_clin.writerow([sgrna, sgrna_strand, sgrna_context, chromosome, gene_name, gene_strand, edit, transcript_ref_allele,
										 transcript_alt_allele, genome_ref_allele, genome_alt_allele] + snp_row)
			pos = pos+1
		current_exon += 1
	return 1  


def get_seq_info(seq):
	exons = []
	cds = ''
	cur_val = 0
	ex = '0'
	for i,p in enumerate(seq):
		if p.islower() and seq[i-1].isupper():
			ex = ex +':'+ str(i-1)
			exons.append(ex)
			cds = cds + seq[cur_val:i]
		elif p.islower() and seq[i+1].isupper():
			cur_val = i+1
			ex = str(cur_val)
		else:
			continue
	ex = ex +':' +str(i)
	cds += seq[cur_val:]
	exons.append(ex)
	return exons, cds


def get_seq_pro(cds, codon_map):
	i=0
	pro = ''
	while(i < len(cds)-3):
		codon = cds[i:(i+3)]
		aa = codon_map[codon]
		pro += aa
		i+=3
	return pro


def write_readme(output_folder, input_file, pam, window, edit, variant_file, intron_buffer, filter_gc):
	with open(output_folder+'/README.txt','w') as o:
		w = csv.writer(o,delimiter='\t')
		w.writerow((['Input file: '+input_file]))
		w.writerow((['PAM: '+pam]))
		w.writerow((['Edit window: '+ window]))
		w.writerow((['Edit: '+edit]))
		w.writerow((['Intron Buffer: ' + str(intron_buffer)]))
		w.writerow((['Filter out GC motifs: ' + str(filter_gc)]))
		w.writerow((['Variant file: ' + variant_file]))
		w.writerow((['Output folder: '+output_folder]))
	return


def read_args(args):
	input_type = args.input_type
	input_file = args.input_file
	intron_buffer = args.intron_buffer
	if args.filter_gc == 'True':
		filter_gc = True
	else:
		filter_gc = False
	if input_type == 'tid':
		input_df = pd.read_table(args.input_file)
	else:
		input_df = pd.DataFrame(columns=['Sequence', 'ID'])
		fasta_seqs = SeqIO.parse(open(args.input_file), 'fasta')
		i=0
		for fasta in fasta_seqs:
			input_df.loc[i, 'ID'] = fasta.id
			input_df.loc[i, 'Sequence'] = str(fasta.seq)
			i+=1
	if args.be_type != '':
		pam, window, sg_len, edit = get_pam_window_len(args.be_type)
	else:
		pam = args.pam
		window = args.edit_window
		sg_len = args.sg_len
		edit = args.edit
	variant_file = args.variant_file
	variant_df = pd.read_table(variant_file, dtype = {'#AlleleID':str, 'Name':str, 'ClinicalSignificance':str, 'Assembly':str,
													  'Chromosome':str,'Type':str,'Start':int,'ReferenceAlleleVCF':str,
													  'AlternateAlleleVCF':str,'ReviewStatus':str}, usecols = ['#AlleleID',
																											'GeneSymbol',
																											'Name',
																											'ClinicalSignificance',
																											'PhenotypeList',
																											'Assembly',
																											'Chromosome',
																											'Type',
																											'Start',
																											'ReferenceAlleleVCF',
																											'AlternateAlleleVCF',
																											'ReviewStatus'])
	parsed_variant_df = parse_variant_df(variant_df)
	output_name = args.output_name
	output_folder = output_name + '_' + str(datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)
	codon_map = get_codon_map()
	aa_map = get_aa_map()
	return input_type, input_file, input_df, pam, window, sg_len, edit, parsed_variant_df, variant_file, output_name, output_folder, codon_map, aa_map, intron_buffer, filter_gc


def check_sequences(seq):
	m = re.search('[^ACTGactg]',seq)
	if m is not None:
		error = 'Sequences contain non-ACTG characters'
	else:
		error = ''
	return error


if __name__ == '__main__':
	args = get_parser().parse_args()
	input_type, input_file, input_df, pam, window, sg_len, edit, parsed_variant_df, variant_file, output_name, output_folder, codon_map, aa_map, intron_buffer, filter_gc = read_args(args)
	write_readme(output_folder, input_file, pam, window, edit, variant_file, intron_buffer, filter_gc)

	with open(output_folder+'/sgrna_designs_'+output_name+'.txt','w') as o_sgrna, open(output_folder + '/error_report.txt', 'w') as o_error, open(output_folder+'/clinvar_annotations_'+output_name+'.txt','w') as o_clin:
		w_error = csv.writer(o_error,delimiter='\t')
		w_error.writerow(['Gene Symbol','Ensembl transcript ID','sgRNA','sgRNA Strand','Error'])
		w = csv.writer(o_sgrna,delimiter='\t')
		w.writerow((['sgRNA sequence','sgRNA context sequence','Gene Symbol','Ensembl Gene ID','Ensembl transcript ID','Gene strand','Genome assembly',
					 'Transcript reference allele','Transcript alternate allele','Genome reference allele',
					 'Genome alternate allele', 'Chromosome', 'sgrna genomic position', 'sgRNA Strand','PAM',
					 'Edit','# edits','#silent edits', 'Nucleotide edits','Amino acid edits',
					 'Mutation category', 'Clinical significance', 'BsmBI flag', '4T flag']))
		w_clin = csv.writer(o_clin, delimiter='\t')
		w_clin.writerow(['sgRNA sequence', 'sgRNA strand', 'sgRNA context sequence','Chromosome','Gene', 'Gene strand', 'Edit', 
						 'Transcript reference allele','Transcript alternate allele', 'Genome reference allele', 'Genome alternate allele',
						 'Edit nucleotide', 'Edit nucleotide position(s)', 'sgRNA amino acid change', 'Original codon',
						 'Edited codon', 'Mutation category','SNP amino acid change', 'SNP name', 'SNP clinical significance',
						 'SNP nucleotide position', 'SNP reference allele', 'SNP alternate allele', 'SNP review status',
						 'Same nucleotide position', 'Same nucleotide change', 'Same amino acid position',
						 'Same amino acid change'])
		for i,r in input_df.iterrows():
			print('Designing for '+r[1])
			tr = r[0]
			if input_type == 'tid':
				gene_name, assembly, gene_strand, chromosome, gene_id, exons, cds_map, abs_pos_map, fs, utr, cds_start_exon, utr5_flag, utr3_flag = get_tr_info(tr, input_type)
				if exons == '':
					error = 'Transcript not found'
					print(error)
					w_error.writerow([gene_name,tr,'N/A','N/A',error])
					continue
				tr_seq = get_tr_sequence(tr)
				seq_error = check_sequences(tr_seq)
				if seq_error != '':
					print(seq_error)
					w_error.writerow([gene_name, tr, 'N/A', 'N/A', seq_error])
					continue
				if tr_seq == '':
					error = 'Transcript sequence not found'
					print(error)
					w_error.writerow([gene_name, tr, 'N/A', 'N/A', error])
					continue

				pro_sequence = get_pro_sequence(tr)
				if pro_sequence == '':
					error = 'Protein sequence not found'
					print(error)
					w_error.writerow([gene_name, tr, 'N/A', 'N/A', error])
					continue
				cds_sequence = get_cds_sequence(tr)
				seq_error = check_sequences(cds_sequence)
				if seq_error != '':
					print(seq_error)
					w_error.writerow([gene_name, tr, 'N/A', 'N/A', seq_error])
					continue
				if cds_sequence == '':
					error = 'Coding sequence not found'
					print(error)
					w_error.writerow([gene_name, tr, 'N/A', error])
					continue
				gene_variant_df = parsed_variant_df[parsed_variant_df.GeneSymbol == gene_name]
			elif input_type == 'nuc':
				tr_seq = r[0]
				seq_error = check_sequences(tr_seq)
				if seq_error != '':
					print(seq_error)
					w_error.writerow([r[1], r[1], 'N/A', 'N/A', seq_error])
					continue
				gene_name, assembly, chromosome, gene_id,gene_start, gene_end, gene_strand, tr = r[1], '', '', '',0, len(tr_seq), 1, r[1]
				exons, cds_sequence = get_seq_info(r[0])
				pro_sequence = get_seq_pro(cds_sequence, codon_map)
				abs_pos_map, fs = get_absolute_pos(gene_start, gene_end, gene_strand, input_type)
				cds_map = get_cds_map(exons, gene_strand)
				utr = {}
				gene_variant_df = parsed_variant_df[parsed_variant_df.GeneSymbol == '']
				cds_start_exon = 1
			else:
				print('Please enter a valid input type; tid for list of transcripts or nuc for fasta sequences')
			if edit == 'all':
				edit_1 = 'C-T'
				val = design_sgrnas(gene_name, assembly, chromosome, gene_id, w, tr_seq, abs_pos_map, fs, cds_map, utr, tr, pam, exons, gene_strand, edit_1, window, cds_sequence, len(pam), sg_len, cds_start_exon, w_error, w_clin, gene_variant_df, aa_map, input_type, intron_buffer, filter_gc)
				edit_2 = 'A-G'
				val = design_sgrnas(gene_name, assembly, chromosome, gene_id, w, tr_seq, abs_pos_map, fs, cds_map, utr, tr, pam, exons, gene_strand, edit_2, window, cds_sequence, len(pam), sg_len, cds_start_exon, w_error, w_clin, gene_variant_df, aa_map, input_type, intron_buffer, filter_gc)
			else:
				val = design_sgrnas(gene_name, assembly, chromosome, gene_id, w, tr_seq, abs_pos_map, fs, cds_map, utr, tr, pam, exons, gene_strand, edit, window, cds_sequence, len(pam), sg_len, cds_start_exon, w_error, w_clin, gene_variant_df, aa_map, input_type, intron_buffer, filter_gc)

			if input_type !='tid':
				os.remove(output_folder+'/clinvar_annotations_'+output_name+'.txt')

