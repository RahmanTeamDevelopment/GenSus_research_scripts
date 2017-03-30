#!/env/bin/python


from __future__ import division
from collections import defaultdict
import os
import sys
import datetime
import optparse
import subprocess



path_samples = "/scratch/cancgene/syost/current_working_set_of_exome_calls_files/TXTs/"  # path will be changed in the future

#  Stable input files:
path_countICR1000 = "/scratch/cancgene/gst/databases/20170224_ICR1000_unionOfAllVars_counts.txt"
path_InHouse419 = "/scratch/cancgene/gst/databases/20170224_InHouse419_unionOfAllVars_counts.txt"
path_exac = "/scratch/cancgene/gst/databases/20170223_ExAC.r0.3.nonTCGA.sites_cavaOutput_v1.2.txt"  # NOTE: cavaOutput_v1.2
path_GeneSummaryPTV_ICR1000 = "/scratch/cancgene/gst/databases/20170302_GeneSummaryPTV_ICR1000.txt"
path_GeneSummaryPTV_InHouse419 = "/scratch/cancgene/gst/databases/20170302_GeneSummaryPTV_InHouse419.txt"


# Printing out welcome message
def printStartTime():
	start_time = datetime.datetime.now()
	print "\n--------------------------------------------------------------------------------------------"
	print 'variant_inheritance_set_of_genes.py is now running.'
	print 'Started: ', str(start_time)
	return start_time

# Printing out summary message
def printEndTime(start_time):
	end_time = datetime.datetime.now()
	print "\n--------------------------------------------------------------------------------------------"
	print 'Finished successfully.'
	print 'Ended: ', str(end_time)
	print 'Total runtime: ' + str(end_time - start_time)
	print "\n--------------------------------------------------------------------------------------------"

# ICR1000 variant count
def countICR1000():
	icr1000_dict = dict()
	for line in open(path_countICR1000):
		cols=line.split()
		if line[0] == '#': continue
		key = cols[5] + '$' + cols[9]  # TRANSCRIPT+CSN unique combination
		TOT = cols[19]  # TOT column in ICR1000 counts file
		HOM = cols[23]  # HOM column in ICR1000 counts file
		MAF = ((float(TOT) + float(HOM)) / (993 * 2)) * 100  # Minor Allele Frequency(%)
		value = (TOT,MAF)
		icr1000_dict[key] = value  # Variant count,MAF
	return icr1000_dict

# InHouse419 variant count
def countInHouse419():
	inHouse419_dict = dict()
	for line in open(path_InHouse419):
		cols = line.split()
		key = cols[5] + '$' + cols[9]  # TRANSCRIPT+CSN unique combination
		value = cols[19]  # TOT column in inHouse419 counts file.
		inHouse419_dict[key] = value
	return inHouse419_dict

# ExAC variant count
def countEXAC():
	exac_dict = dict()
	for line in open(path_exac):
		if line[0] == '#': continue
		cols = line.split()
		AC_Adj = cols[21]
		AC_Hom = cols[26]
		AC_NFE = cols[31]
		HOM_NFE = cols[52]
		AN_Adj = cols[23]
		tot_variantCount = int(AC_Adj) - int(AC_Hom)  # Total variant count
		nfe_variantCount = int(AC_NFE) - int(HOM_NFE) # Non Finish European variant count
		try:
			MAF = (float(AC_Adj) / float(AN_Adj)) * 100 # Minor Allele Frequency (%)
		except:
			MAF = 0
		key = cols[6] + '$' + cols[10] # TRANSCRIPT+CSN unique combination
		pair = (tot_variantCount, nfe_variantCount, MAF)
		exac_dict[key] = pair
	return exac_dict


# Read input with samples keys to dictionary
def read_key_file(input_path):
	input_dict = defaultdict(list)
	fin_samplesList = open(input_path)
	for line in fin_samplesList:
		cols = line.split('\t')
		key = cols[1] # Family ID
		value = (cols[0], cols[2], cols[3].split('\n')[0])  # Sample_ID, Relationship, Affected/Unaffected
		input_dict[key].append(value)
	fin_samplesList.close()
	return input_dict
	
# Writing variants that meet the requirements to VariantInheritanceSummary_setOfGenes.txt
def extract_variants(input_path, icr1000_dict, inHouse419_dict, exac_dict, input_dictm, genes_list):

	flag = 0
	first_line = 1
	fin_samplesList = open(input_path)
	out_pathfile = os.getcwd() + "/VariantInheritanceSummary_setOfGenes.txt"
	
	for point_line in fin_samplesList:
		if first_line == 1:
			if point_line[0] != '#':
				sys.exit("Error: input file " + input_path + " must have header starting with # ")
			else:
				rest = point_line.split('\t')[1:]
				header_list_of_samples = '\t'.join(rest)
				first_line = 0
				fout = open(out_pathfile, 'w+')
				continue
				
		sample_ID = point_line.split('\t')[0]
		family_ID = point_line.split('\t')[1]
		relationship = point_line.split('\t')[2]
		if 'Proband' not in relationship: continue
		
		rest = point_line.split('\t')[1:]
		additional_columns = '\t'.join(rest)
		path_to_txt = path_samples + sample_ID + "_annotated_calls.txt"
		if not os.path.isfile(path_to_txt):
			sys.exit("Error: Sample " + sample_ID + " has no corresponding " + sample_ID + "_annotated_calls.txt file in: " + path_samples) 
		fin_sample_txt = open(path_to_txt)
		for line in fin_sample_txt:
			if line[0] == '#' and flag == 0:  # First txt file in the list
				header = line.split()
				fout.write('\t'.join(header[0:20]) + '\t' + '\t'.join(header[23:26]) + '\t' + "CountICR1000\tMAF_ICR1000\tCountInHouse419\tCountExacNFE\tCountExacTotal\tMAF_ExAC\tTrio_analyzed\tNumSamplesPerFamily\tInheritance\t" + header_list_of_samples)
				flag = 1
				continue
			cols = line.split()
			chr = cols[0]
			pos = cols[1]
			ref = cols[2]
			alt = cols[3]
			qual = cols[4]
			qualflag = cols[5]
			filter = cols[6]
			TR = cols[7] 
			TC = cols[8]
			sample = sample_ID
			GT = cols[10]
			TYPE = cols[11]
			ENST = cols[12]
			GENE = cols[13]
			
			if GENE not in genes_list: continue
			
			TRINFO = cols[14]
			LOC = cols[15]
			CSN = cols[16]
			CLASS = cols[17]
			SO = cols[18]
			IMPACT = cols[19]
			
			if IMPACT == '3' and CLASS != 'SY': continue
			
			ALTANN = cols[23]
			ALTCLASS = cols[24]
			ALTSO = cols[25]
			
			try:
				count_ICR1000 = icr1000_dict[ENST + '$' + CSN][0]  # variant count
				MAF_ICR1000 = icr1000_dict[ENST + '$' + CSN][1]  # Minor Allele Frequency
				if MAF_ICR1000 > 0.5:continue
			except:
				count_ICR1000 = '0'
				MAF_ICR1000 = '0'
                                pass
								
			try:	
				count_InHouse419 = inHouse419_dict[ENST + '$' + CSN]  # Variant count
			except:
				count_InHouse419 = '0'
                                pass
								
			try:
				countExACNFE = exac_dict[ENST+'$'+CSN][1]  # Variant count
				countExACTotal = exac_dict[ENST+'$'+CSN][0]  # Variant count
				MAF_ExAC = exac_dict[ENST+'$'+CSN][2]  # Minor Allele Frequency
				if MAF_ExAC > 0.5:continue
			except:
				countExACNFE = '0'
				countExACTotal = '0'
				MAF_ExAC = '0'
				pass
				
			members_trio = 0
			affected_list = []
			genotype_list = []
			num_samples_per_family = len(input_dict[family_ID])
			found_match_in_parent = 0
			list_sex_of_match = [] # For X-linked genes
			
			for n in range(0, num_samples_per_family):
				if 'Father' in input_dict[family_ID][n][1] or 'Mother' in input_dict[family_ID][n][1] or 'Proband' in input_dict[family_ID][n][1]:
					members_trio = members_trio + 1	
					
				if 'Father' in input_dict[family_ID][n][1] or 'Mother' in input_dict[family_ID][n][1]:
					sampleID_parent = input_dict[family_ID][n][0]
					fin_txt_parent = open(path_samples + sampleID_parent + "_annotated_calls.txt")
					
					for line in fin_txt_parent:
						cols_parent = line.split()
						gene_parent = cols_parent[13]
						csn_parent = cols_parent[16]
						genotype_list.append(cols_parent[10])
						if GENE == gene_parent and CSN == csn_parent:
							found_match_in_parent = found_match_in_parent + 1
							list_sex_of_match.append(input_dict[family_ID][n][1])
							affected_list.append(input_dict[family_ID][n][2])
							continue
					fin_txt_parent.close()
					
			if members_trio == 1:
				inheritanceInfo_column = 'no_parents' 
				is_trio = 'no'
			elif members_trio == 2:  
				is_trio = 'no'
				if found_match_in_parent == 0:
					inheritanceInfo_column = 'not_in_only_parent'
				else:
					if 'Unaffected' in affected_list:
						inheritanceInfo_column = 'inherited_from_only_unaffected'
					elif 'Affected' in affected_list:
						inheritanceInfo_column = 'inherited_from_only_affected'
					elif 'Unknown' in affected_list:
						inheritanceInfo_column = 'inherited_from_only_unknown'
					else:
						print ("\nWarning < Sample_ID: "+sample_ID+" > inheritance column in input file must contain only one of the following: Unaffected/Affected/Unknown")
			elif members_trio >= 3:
				is_trio = 'yes'
				if found_match_in_parent == 0:
					inheritanceInfo_column = 'de_novo'
				else:
					if 'Unknown' in affected_list and 'affected' not in affected_list:
						inheritanceInfo_column = 'inherited_from_unknown'
					elif 'Affected' not in affected_list and 'Unknown' not in affected_list:
						inheritanceInfo_column = 'inherited_from_unaffected'
					else:
						if found_match_in_parent == 1 or (found_match_in_parent == 2 and 'unaffected' not in affected_list):
							inheritanceInfo_column = 'inherited_from_affected'
						if found_match_in_parent == 2:  # both parents have a match
							
							# Special case where both parents are het and only one is affected:
							if '1/1' not in genotype_list:
								inheritanceInfo_column = 'inherited_from_parent' 
							else:
								# When one of the parents is hom we know the affection status is his
								index_hom = genotype_list.index('1/1')
								affection_of_hom = affected_list[index_hom]
								if affection_of_hom == 'Affected':
									inheritanceInfo_column = 'inherited_from_affected'
								else:
									inheritanceInfo_column = 'inherited_from_unaffected'
			if chr == 'X' and 'inherited' in inheritanceInfo_column and 'parent' not in inheritanceInfo_column and found_match_in_parent == 1:
				inheritanceInfo_column = inheritanceInfo_column + '_' + list_sex_of_match[0]
				
			out=[chr, pos, ref, alt, qual, qualflag, filter, TR, TC, sample, GT, TYPE, ENST, GENE, TRINFO, LOC, CSN, CLASS, SO, IMPACT, ALTANN, ALTCLASS, ALTSO, count_ICR1000, str(MAF_ICR1000), count_InHouse419, str(countExACNFE), str(countExACTotal), str(MAF_ExAC), is_trio, str(num_samples_per_family), inheritanceInfo_column]
			# Writing variant info to output file
			fout.write('\t'.join(out) + '\t' + additional_columns)
		
		# Closing input file
		fin_sample_txt.close()
		
		affected_list = []
		genotype_list = []
		list_sex_of_match = []
		inheritanceInfo_column = ''
		
	# Closing output file
	fout.close()
		

#---------------------------------------------------------------------------------------------------
#----------------------------------------      MAIN      -------------------------------------------
#---------------------------------------------------------------------------------------------------
if __name__ == "__main__":
	
	# Parse command line arguments 
	descri = 'variant_inheritance_set_of_genes.py'
	parser = optparse.OptionParser(usage='python path/to/variant_inheritance_set_of_genes.py <options>', description=descri)
	parser.add_option('-I', dest='input', help='input TSV file of 4 columns: SampleID, FamilyID, Relationship (Father/Mother/Proband/...), Affected?(Affected,Unaffected,Unknown). File MUST have header starting with #')	
	parser.add_option('--genes', dest='path_genes_list', help='file with genes list, each gene in a separate line')
	
	# Get and check arguments 
	(opts, args) = parser.parse_args()
	input = opts.input
	path_genes_list = opts.path_genes_list

	start_time = printStartTime()
	
	if not input:
		sys.exit("\nError: please enter input file with a list of samples")
	if not path_genes_list:
		sys.exit("\nError: please enter a file with a list of candidate genes")
	if not os.path.isfile(input):
		sys.exit("\nError: file "+input+" doesn't exist")
	if not os.path.isfile(path_genes_list):
		sys.exit("\nError: file "+path_genes_list+" doesn't exist")

	
	# Creating dictionaries from stable input files
	icr1000_dict = countICR1000()
	inHouse419_dict = countInHouse419()
	exac_dict = countEXAC()
	
	# Creating a dictionary from input file with sample_IDs
	input_dict = read_key_file(input)
	
	# Reading genes list from genes input file
	genes_list = []
	for line in open(path_genes_list):
		genes_list.append(line.split('\n')[0])
	
	# Writing variants that meet the requirements to file
	extract_variants(input, icr1000_dict, inHouse419_dict, exac_dict, input_dict, genes_list)


	printEndTime(start_time)


