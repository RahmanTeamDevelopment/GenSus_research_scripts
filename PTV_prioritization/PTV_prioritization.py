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

default_out_path = os.getcwd()


# Printing out welcome message
def printStartTime():
	start_time = datetime.datetime.now()
	print "\n--------------------------------------------------------------------------------------------"
	print 'PTV prioritization.py is now running.'
	print 'Started: ', str(start_time)
	return start_time

# Printing out summary message
def printEndTime(start_time):
	end_time = datetime.datetime.now()
	print 'Finished successfully.'
	print 'Ended: ', str(end_time)
	print 'Total runtime: ' + str(end_time - start_time)
	print "\n--------------------------------------------------------------------------------------------"

# ICR1000 variant count
def countICR1000():
	icr1000_dict = dict()
	for line in open(path_countICR1000):
		cols = line.split()
		key = cols[5]+'$'+cols[9]  # TRANSCRIPT+CSN unique combination
		value = cols[19]  # TOT column in ICR1000 counts file.
		icr1000_dict[key] = value
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
		cols=line.split()
		AC_Adj = cols[21]
		AC_Hom = cols[26]
		AC_NFE = cols[31]
		HOM_NFE = cols[52]
		tot_variantCount = int(AC_Adj) - int(AC_Hom)  # Total variant count
		nfe_variantCount = int(AC_NFE) - int(HOM_NFE)  # Non Finish European variant count
		key = cols[6] + '$' + cols[10]  # TRANSCRIPT+CSN unique combination
		pair = (tot_variantCount, nfe_variantCount)
		exac_dict[key] = pair
	return exac_dict

	
# Outputs SampleSummaryPTV.txt
def run_step_1(input_path, out_path, name,icr1000_dict, inHouse419_dict, exac_dict):

	flag = 0
	first_line = 1
	fin_samplesList = open(input_path)
	out_pathfile = out_path+"/SampleSummaryPTV"+name+".txt"
	
	for point_line in fin_samplesList:
		if first_line == 1:
			if point_line[0] != '#':
				sys.exit("Error: input file " + input_path + " must have header starting with # ")
			else:
				rest = point_line.split('\t')[1:]
				header_list_of_samples = '\t'.join(rest)
				if header_list_of_samples: 
					end = ''
				else:
					end = '\n'
				first_line = 0
				fout = open(out_pathfile,'w+')
				continue
		if point_line in ['\n', '\r\n']:
			continue
			
		sample_ID = point_line.split('\t')[0]
		if sample_ID[-1] == '\n':
			sample_ID = sample_ID[:-1]
		rest = point_line.split('\t')[1:]
		additional_columns = '\t'.join(rest)
		path_to_txt = path_samples + sample_ID + "_annotated_calls.txt"
		if not os.path.isfile(path_to_txt):
			sys.exit("Error: Sample " + sample_ID + " has no corresponding " + sample_ID + "_annotated_calls.txt file in: " + path_samples) 
			
		fin_sample_txt = open(path_samples + sample_ID + "_annotated_calls.txt")
		for line in fin_sample_txt:
			if line[0] == '#' and flag == 0: #first txt file in the list
				header = line.split()
				fout.write('\t'.join(header[0:20]) + '\t' + '\t'.join(header[23:26]) + '\t' + "CountICR1000\tCountInHouse419\tCountExacNFE\tCountExacTotal\t" + end + header_list_of_samples)
				flag = 1
				continue
			cols = line.split()
			chr = cols[0]
			pos = cols[1]
			ref = cols[2]
			alt = cols[3]
			if len(ref) >= 11 or len(alt) >= 11: continue
			qual = cols[4]
			qualflag = cols[5]
			if qualflag != "high": continue
			filter = cols[6]
			TR = cols[7] 
			if int(TR) < 10: continue  # Number of reads supporting the alternate allele
			TC = cols[8]
			if int(TC) < 20: continue  # Number of reads covering the position in the genome 
			sample = sample_ID
			GT = cols[10]
			TYPE = cols[11]
			ENST = cols[12]
			GENE = cols[13]
			TRINFO = cols[14]
			LOC = cols[15]
			CSN = cols[16]
			CLASS = cols[17]
			SO = cols[18]
			IMPACT = cols[19]
			if IMPACT != '1': continue
			ALTANN = cols[23]
			ALTCLASS = cols[24]
			ALTSO = cols[25]
			try:
				count_ICR1000 = icr1000_dict[ENST + '$' + CSN]  # variant count
				if int(count_ICR1000) > 1:continue
			except:
				count_ICR1000 = '0'
                                pass
			try:	
				count_InHouse419 = inHouse419_dict[ENST + '$' + CSN]  # variant count
			except:
				count_InHouse419 = '0'
                                pass
			try:
				countExACNFE = exac_dict[ENST + '$' + CSN][1] #  variant count
				countExACTotal = exac_dict[ENST + '$' + CSN][0] #  variant count
				if countExACNFE > 10: continue
			except:
				countExACNFE = '0'
				countExACTotal = '0'
				pass
				
			out = [chr, pos, ref, alt, qual, qualflag, filter, TR, TC, sample, GT, TYPE, ENST, GENE, TRINFO, LOC, CSN, CLASS, SO, IMPACT, ALTANN, ALTCLASS, ALTSO, count_ICR1000, count_InHouse419, str(countExACNFE), str(countExACTotal)]
			to_write = '\t'.join(out) + '\t' + additional_columns
			if to_write[-1] != '\n': to_write = to_write + '\n'
			
			# Writing variant info to output file
			fout.write(to_write)
		
		# Closing input file
		fin_sample_txt.close()
		
	# Closing output file
	fout.close()
		
# SampleSummaryPTV.txt count
def countCase_for_step_2(path_in_file):
	sampleSummary_dict = defaultdict(list)
	if not os.path.isfile(path_in_file):
		sys.exit("Error: file " + path_in_file + " doesn't exist")
		
	in_path_to_SampleSummaryPTV = open(path_in_file)
	for line in in_path_to_SampleSummaryPTV:
		cols = line.split()
		key = cols[12] + '$' + cols[16]  # ENST$CSN
		value = cols[9]  # sample_id
		try:
			sampleSummary_dict[key].append(value)
		except:
			sampleSummary_dict[key] = value
	in_path_to_SampleSummaryPTV.close()
	return sampleSummary_dict

# Outputs VariantSummaryPTV.txt
def run_step_2(out_path, path_fin, name, sampleSummary_dict):
	fin_SampleSummary = open(path_fin)
	out_pathfile = out_path + "/VariantSummaryPTV" + name + ".txt"
	fout = open(out_pathfile, 'w+')
	header = ['#CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'ENST', 'GENE', 'TRINFO', 'LOC', 'CSN', 'CLASS', 'SO', 'IMPACT', 'ALTANN', 'ALTCLASS', 'ALTSO', 'CountCase', 'countICR1000', 'countInHoluse419', 'countExACNFE', 'countExACTotal']
	fout.write('\t'.join(header) + '\n')
	
	for line in fin_SampleSummary:
		cols = line.split()
		sample_id = cols[9]
		ENST = cols[12]
		CSN = cols[16]
		if sampleSummary_dict[ENST + '$' + CSN][0] != sample_id: continue 
		chr = cols[0]
		pos = cols[1]
		ref = cols[2]
		alt = cols[3]
		type = cols[11]
		ENST = cols[12]
		GENE = cols[13]
		TRINFO = cols[14]
		LOC = cols[15]
		CSN = cols[16]
		CLASS = cols[17]
		SO = cols[18]
		IMPACT = cols[19]	
		ALTANN = cols[20]
		ALTCLASS = cols[21]
		ALTSO = cols[22]
		countCase = str(len(sampleSummary_dict[ENST + '$' + CSN]))  # Number of times enst+csn is seen in SampleSummaryPTV.txt
		count_ICR1000 = cols[23]
		count_InHouse419 = cols[24]
		count_ExACNFE = cols[25]
		count_ExACTotal = cols[26]
		out = [chr, pos, ref, alt, type, ENST, GENE, TRINFO, LOC, CSN, CLASS, SO, IMPACT, ALTANN, ALTCLASS, ALTSO, countCase, count_ICR1000, count_InHouse419, count_ExACNFE, count_ExACTotal]
		if line[0] != '#': fout.write('\t'.join(out) + '\n')
		
	fin_SampleSummary.close()
	fout.close()

# VariantSummaryPTV.txt count
def countCase_for_step_3(path_in_file):
	variantSummary_dict = defaultdict(list)
	if not os.path.isfile(path_in_file):
		sys.exit("Error: file " + path_in_file + " doesn't exist")
	in_path_to_VariantSummaryPTV = open(path_in_file)
	for line in in_path_to_VariantSummaryPTV:
		cols = line.split()
		key = cols[5]
		value = (cols[9], cols[16])  # CSN, CountCase
		try:
			variantSummary_dict[key].append(value)
		except:
			variantSummary_dict[key] = value
	in_path_to_VariantSummaryPTV.close()
	return variantSummary_dict

# Create a dictionary for GeneSummaryPTV_ICR1000.txt
def gene_summary_ICR1000_for_step_3(path_GeneSummaryPTV_ICR1000):
	gene_summary_icr1000_dict = dict()
	if not os.path.isfile(path_GeneSummaryPTV_ICR1000):
		sys.exit("Error: file " + path_GeneSummaryPTV_ICR1000 + " doesn't exist")
	f_in = open(path_GeneSummaryPTV_ICR1000)
	for line in f_in:
		cols = line.split()
		if cols[0] == 'GENE': continue
		key = cols[2]
		value = (cols[3], cols[4], cols[5])  # CaseTotalRarePTV, CaseDifferentRarePTV, CaseSingletonPTV
		gene_summary_icr1000_dict[key] = value
	f_in.close()
	return gene_summary_icr1000_dict

# Create a dictionary for GeneSummaryPTV_InHouse419.txt
def gene_summary_InHouse419_for_step_3(path_GeneSummaryPTV_InHouse419):
	gene_summary_inHouse419_dict = dict()
	if not os.path.isfile(path_GeneSummaryPTV_InHouse419):
		sys.exit("Error: file " + path_GeneSummaryPTV_InHouse419 + " doesn't exist")
	f_in = open(path_GeneSummaryPTV_InHouse419)
	for line in f_in:
		cols = line.split()
		if cols[0] == 'GENE': continue
		key = cols[2]
		value = (cols[3], cols[4], cols[5])  # CaseTotalRarePTV, CaseDifferentRarePTV, CaseSingletonPTV
		gene_summary_inHouse419_dict[key] = value
	f_in.close()
	return gene_summary_inHouse419_dict

# Outputs GeneSummaryPTV.txt
def run_step_3(out_path, path_fin, name, variantSummary_dict, geneSummary_Icr1000_dict, geneSummary_InHouse419_dict):
	fin_VariantSummary = open(path_fin)
	out_pathfile = out_path + "/GeneSummaryPTV" + name + ".txt"
	fout = open(out_pathfile, 'w+')
	header = ['GENE', 'TRINFO', 'ENST', 'CaseTotalRarePTV', 'CaseDifferentRarePTV', 'CaseSingletonPTV', 'ICR1000_TotalRarePTV', 'ICR1000_DifferentRarePTV', 'ICR1000_SingletonPTV', 'InHouse419_TotalRarePTV', 'InHouse419_DifferentRarePTV', 'InHouse419_SingletonPTV']
	fout.write('\t'.join(header)+'\n')
	
	for line in fin_VariantSummary:
		if line[0] == '#': continue
		cols = line.split()
		GENE = cols[6]
		TRINFO = cols[7]
		ENST = cols[5]
		CSN = cols[9]
		if variantSummary_dict[ENST][0][0] != CSN: continue  
		sum_of_countCase = 0 #  CaseTotalRarePTV
		num_of_countCase_rows = len(variantSummary_dict[ENST])  # CaseDifferentRarePTV
		sum_of_singleton_countCase = 0  # CaseSingletonPTV
		
		for i in range(0,num_of_countCase_rows):
			countCase_enst_csn = variantSummary_dict[ENST][i][1]
			sum_of_countCase = sum_of_countCase + int(countCase_enst_csn)
			if countCase_enst_csn == '1':
				sum_of_singleton_countCase = sum_of_singleton_countCase + 1
		blank_cols = '\t'.join(('0') * 3)
		try:
			ICR1000_out = '\t'.join([geneSummary_Icr1000_dict[ENST][0], geneSummary_Icr1000_dict[ENST][1], geneSummary_Icr1000_dict[ENST][2]])
		except:
			ICR1000_out = blank_cols
		try:
			InHouse419_out = '\t'.join([geneSummary_InHouse419_dict[ENST][0], geneSummary_InHouse419_dict[ENST][1], geneSummary_InHouse419_dict[ENST][2]])
		except:
			InHouse419_out = blank_cols	
		out = [GENE, TRINFO, ENST, str(sum_of_countCase), str(num_of_countCase_rows), str(sum_of_singleton_countCase)]		
		fout.write('\t'.join(out) + '\t' + ICR1000_out + '\t' + InHouse419_out + '\n')
		
	fin_VariantSummary.close()
	# Closing output file
	fout.close()

def get_input_path(out_path, file_name_string):
		name_file = ''
		filenames = os.listdir(out_path)
		for f in filenames:
			if file_name_string in f:
				name_file = f
				break
		if name_file == '':
			sys.exit("Error: there is no " + file_name_string + " file in the current directory")
		return name_file


#---------------------------------------------------------------------------------------------------
#----------------------------------------      MAIN      -------------------------------------------
#---------------------------------------------------------------------------------------------------
if __name__ == "__main__":
	# Parse command line arguments 
	descri = 'PTV_prioritization.py'
	parser = optparse.OptionParser(usage='python path/to/PTV_prioritization.py <options>', description=descri)
	steps = {'0':'Run steps 1,2,3 ','1':'Sample_Summary_PTV,','2':'Variant_Summary_PTV ','3':'Gene_Summiary_PTV '}
	parser.add_option('-s', dest='step', help='0:Run steps 1,2,3 consecutively 1:Sample_Summary_PTV 2:Variant_Summary_PTV 3:Gene_Summary_PTV')
		
	step_0_1_option=optparse.OptionGroup(parser, "steps 0 and 1 argument")
	step_0_1_option.add_option('-I', dest='input', help='input file: Tab-separated list of samples with annotated_calls.txt files. More columns are optional. File MUST have header starting with #')
	parser.add_option_group(step_0_1_option)
	
	step0_options=optparse.OptionGroup(parser, "step 0 optional parameters")
	step0_options.add_option('-o', dest='out_path', help='directory where output files will be stored [cwd]')
	step0_options.add_option('--name', dest='name', help='add a representative name to out files'\
							'[1-SampleSummaryPTV.txt 2-VariantSummaryPTV.txt 3-GeneSummaryPTV.txt]')
	parser.add_option_group(step0_options)	
	

	# Get and check arguments 
	arguments = ['step'] 
	(opts, args) = parser.parse_args()
	for arg in arguments:
        	if getattr(opts,arg) is None:
            		print "Error: " + arg + " must be specified. Please use -h for details"
            		sys.exit(1)
	step = opts.step 
	input = opts.input
	out_path = opts.out_path 
	name = opts.name
	
	if step not in steps.keys():
		print "Error: step is illegal. Please run with -h for more details"
 	   	sys.exit(1)
	if not out_path: 
		out_path = default_out_path
	else:
		if len(out_path.split('/')) < 2: # A full path is not given 
			out_path = os.getcwd() + '/' + out_path
		if os.path.isdir(out_path) == False:
			os.makedirs(out_path)	
	if not name: 
		name = ''
	else:
		name = '_' + name

		

	start_time = printStartTime()
	if step =='1' or step == '0':
		if not input:
			sys.exit("Error: please enter input file with a list of samples")
		if not os.path.isfile(input):
			sys.exit("Error: file "+input+" doesn't exist")
		print "\n-----------------------------------------------------------------------"
		print "Step 1 running -> Sample Summary PTV"	
		print "-----------------------------------------------------------------------"

		icr1000_dict = countICR1000()
		inHouse419_dict = countInHouse419()
		exac_dict = countEXAC()
		run_step_1(input, out_path, name, icr1000_dict, inHouse419_dict, exac_dict)
	

	if step == '2' or step == '0':
		print "\n-----------------------------------------------------------------------"
		print "Step 2 running -> Variant Summary PTV"
		print "-----------------------------------------------------------------------"
		
		if step == 2:
			name_file = get_input_path(out_path, "SampleSummaryPTV")
			path_fin=out_path + '/' + name_file + name 
		else:
			path_fin = out_path + "/SampleSummaryPTV" + name + ".txt"
		
		sampleSummary_dict = countCase_for_step_2(path_fin)
		run_step_2(out_path, path_fin, name, sampleSummary_dict)

	if step=='3' or step =='0':
		print "\n-----------------------------------------------------------------------"
		print "Step 3 running -> Gene Summary PTV"
		print "-----------------------------------------------------------------------"
		
		if step == 3:
			name_file = get_input_path(out_path, "VariantSummaryPTV")
			path_fin=out_path + '/' + name_file + name
		else:
			path_fin = out_path + "/VariantSummaryPTV" + name + ".txt"
		
		variantSummary_dict = countCase_for_step_3(path_fin)
		geneSummary_Icr1000_dict = gene_summary_ICR1000_for_step_3(path_GeneSummaryPTV_ICR1000)
		geneSummary_InHouse419_dict = gene_summary_InHouse419_for_step_3(path_GeneSummaryPTV_InHouse419)
		run_step_3(out_path, path_fin, name, variantSummary_dict, geneSummary_Icr1000_dict, geneSummary_InHouse419_dict)

	printEndTime(start_time)


