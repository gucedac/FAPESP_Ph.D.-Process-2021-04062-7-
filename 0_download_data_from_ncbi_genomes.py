import os
import sys
import csv
import time
import shutil
import subprocess

def get_defensefinder_systems(files):
	defensefinder_info = []
	for l in csv.reader(open(files),delimiter=','):
		if l[0].startswith('###'):
			data = l[0]
		if data == '###defense_finder_systems':
			defensefinder_info.append(l)

	defensefinder_syst = {}
	for n,i in enumerate(defensefinder_info[2:]): # rows by immune systems
		syst_num  = n+1
		for l_info in i[5].split(','):
			syst_name = i[0]
			gene_syst_locus = l_info
			gene_syst_contig= dic_locations[gene_syst_locus][0]
			gene_syst_start = dic_locations[gene_syst_locus][1]
			gene_syst_end   = dic_locations[gene_syst_locus][2]
			gene_syst_strand= dic_locations[gene_syst_locus][3]
			sys_info = (syst_name,gene_syst_locus,gene_syst_contig,
		           gene_syst_start,gene_syst_end,gene_syst_strand)
			if syst_num in defensefinder_syst:
				defensefinder_syst[syst_num].append(sys_info)
			else:
				defensefinder_syst[syst_num] = [sys_info]
	defensefinder_systems = []
	for k,v in defensefinder_syst.items():
		g_names , g_locuss , g_contigs = [] , [] , []
		g_starts, g_ends   , g_strands = [] , [] , []
		for vi in v:
			g_names.append(vi[0]) ; g_locuss.append(vi[1]) ; g_contigs.append(vi[2])
			g_starts.append(vi[3]); g_ends.append(vi[4])   ; g_strands.append(vi[5])
		defensefinder_row = ['defensefinder',k]+['*'.join(set(g_names))]+\
						                          ['*'.join(set(g_locuss))]+\
						                          ['*'.join(set(g_contigs))]+\
						                          ['*'.join(g_starts)]+\
						                          ['*'.join(g_ends)]+\
						                          ['*'.join(g_strands)]
		defensefinder_systems.append(defensefinder_row)
	return defensefinder_systems

def get_padloc_systems(files):
	padloc_info = []
	for l in csv.reader(open(files),delimiter=','):
		if l[0].startswith('###'):
			data = l[0]
		if data == '###padloc.csv file ':
			padloc_info.append(l)

	padloc_syst = {}
	for i in padloc_info[2:]: # rows by immune gene
		syst_num = i[0]
		syst_name  = i[2]
		gene_syst_locus = i[3]
		gene_syst_contig= dic_locations[gene_syst_locus][0]
		gene_syst_start = dic_locations[gene_syst_locus][1]
		gene_syst_end   = dic_locations[gene_syst_locus][2]
		gene_syst_strand= dic_locations[gene_syst_locus][3]
		sys_info = (syst_name,gene_syst_locus,gene_syst_contig,\
		           gene_syst_start,gene_syst_end,gene_syst_strand)

		if syst_num in padloc_syst:
			padloc_syst[syst_num].append(sys_info)
		else:
			padloc_syst[syst_num] = [sys_info]

	padloc_systems = []
	for k,v in padloc_syst.items():
		g_names , g_locuss , g_contigs = [] , [] , []
		g_starts, g_ends   , g_strands = [] , [] , []
		for vi in v:
			g_names.append(vi[0]) ; g_locuss.append(vi[1]) ; g_contigs.append(vi[2])
			g_starts.append(vi[3]); g_ends.append(vi[4])   ; g_strands.append(vi[5])
		padloc_row = ['padloc',k]+['*'.join(set(g_names))]+\
		                          ['*'.join(set(g_locuss))]+\
		                          ['*'.join(set(g_contigs))]+\
		                          ['*'.join(g_starts)]+\
		                          ['*'.join(g_ends)]+\
		                          ['*'.join(g_strands)]
		padloc_systems.append(padloc_row)
	return padloc_systems


''' >> inputs to start the pipeline << '''
s_running_in = sys.argv[1] # server or 'local
f_genomes_ftp = sys.argv[2] # prokaryotic_assembly_summary_smpl100_genom90498.csv
s_start, s_end = list(map(int,sys.argv[3].split('-'))) # 1-10

''' >> Parent directories of the Immune system project << '''
if s_running_in == 'local':
	d_Work_proj_ImmuSyst = '/media/guillermo/Work/proj_ImmuSyst' ## (scripts)
	d_Data_proj_ImmuSyst = '/media/guillermo/Data/proj_ImmuSyst' ## (raw data)
elif s_running_in == 'server':
	d_Work_proj_ImmuSyst = '/work/xylella/gucedac/proj_ImmuSyst' ## (scripts)
	d_Data_proj_ImmuSyst = '/work/xylella/gucedac/proj_ImmuSyst' ## (raw data)
else:
	print('local or server?')

''' >> Sub-directories of the Immune system project << '''
d_analysis_from_ncbi_geno = os.path.join(d_Work_proj_ImmuSyst,'analysis/A1.Retrieve_defense_system_data/from_ncbi_genomes')
d_data_from_ncbi_geno = os.path.join(d_Data_proj_ImmuSyst,'data/retrieved_data/data_from_ncbi_genomes')

''' >> general directory where the ALL raw data will be obtained << '''
d_destination = d_data_from_ncbi_geno

os.chdir(d_analysis_from_ncbi_geno)

'''#### 1. SELECTION OF THE GENOMES TO WORK (ACCORDING s_start and s_end) ###'''
f_genomes_ftp_path = os.path.join(d_analysis_from_ncbi_geno,f_genomes_ftp)
open_csv = open(f_genomes_ftp_path)
reader_csv = csv.reader(open_csv,delimiter='\t')
rows_set = {}
for row in reader_csv:
	if int(row[0])>=s_start and int(row[0]) <= s_end:
		rows_set[row[0]] = row
open_csv.close()

for index , row in rows_set.items():
	sp = row[1][0]+'_'+row[2].replace(' ','_')
	genome = row[3]
	sp_genome = str(index)+'.'+sp+'_'+genome
	if row[7] == 'Y':
		'''### 2. MAKING THE SPECIES-GENOME FOLDER ###'''
		sp_genome_path = os.path.join(d_destination,sp_genome)
		print('\n***',sp_genome,'***')

		genome_path = sp_genome_path
		if not os.path.isdir(genome_path):
			os.mkdir(genome_path)

		'''### 3 CHECKIN PADLOC AND DEFENSEFINDER RESULTS IN GENOME FOLDER ###'''
		padloc_systems , defensefinder_systems = {},{}
		ispadloc_and_defensefinder_file = []
		l_predictor_totalpredic = []
		for files in os.listdir(genome_path):
			os.chdir(genome_path)
			if files.endswith('_genomic.Padloc'):
				ispadloc_and_defensefinder_file.append(files)
				padloc_info = []
				for l in csv.reader(open(files),delimiter=','):
					if l[0].startswith('###'):
						data = l[0]
					if data == '###padloc.csv file ':
						padloc_info.append(l)
				#padloc_syst = {}
				for i in padloc_info[2:]: # rows by immune gene
					syst_num = i[0]
					syst_name  = i[2]
					gene_syst_locus = i[3]
					padloc_systems['padloc_'+syst_num] = syst_name

			elif files.endswith('_genomic_prodigal.DefenseFinder'):
				ispadloc_and_defensefinder_file.append(files)
				defensefinder_info = []
				for l in csv.reader(open(files),delimiter=','):
					if l[0].startswith('###'):
						data = l[0]
					if data == '###defense_finder_systems':
						defensefinder_info.append(l)
				#defensefinder_systems = {}
				for n,i in enumerate(defensefinder_info[2:]): # rows by immune systems
					syst_num  = n+1
					syst_name = i[0]
					checking_contigs = list('_'.join(c.split('_')[:-1]) for c in i[5].split(','))
					key = 'defensefinder_'+str(syst_num)
					if len(set(checking_contigs)) == 1:
						defensefinder_systems[key] = syst_name
					else:
						print('systems in different contgis:',key,syst_name)
				
			elif files.endswith('fna'):
				predictor , immusyst_numbr = files.split('_')[-2] , files.split('_')[-1].split('.')[0]
				predictor_immusyst_numbr = predictor+'_'+immusyst_numbr
				l_predictor_totalpredic.append(predictor_immusyst_numbr)

		not_processed = []
		for syst in padloc_systems:
			if not syst in l_predictor_totalpredic:
				not_processed.append(syst)
		for syst in defensefinder_systems:
			if not syst in l_predictor_totalpredic:
				not_processed.append(syst)
		if len(not_processed) != 0:
			os.chdir(genome_path)
			for files in os.listdir(genome_path):
				if files.endswith('.Padloc'):
					os.remove(files)

		if (len(not_processed) == 0) and len(ispadloc_and_defensefinder_file)==2:
			print('genome processed!')

		else:
			for files in os.listdir(genome_path):
				os.chdir(genome_path)
				if files.endswith('_PADLOC'):
					pass#shutil.rmtree(files)

			'''### 4. DOWNLOADING THE GENOME FROM NCBI ###'''
			os.chdir(genome_path)
			Genbank_RefSeq_Link = row[6]
			genome_fna = Genbank_RefSeq_Link.split('/')[-1] + '_genomic.fna'
			genome_fna_zip = genome_fna + '.gz'
			genome_fna_link = Genbank_RefSeq_Link + '/' + genome_fna_zip
			if not genome_fna in os.listdir(genome_path):
				subprocess.call('wget -q ' + genome_fna_link , shell=True)
				if genome_fna_zip in os.listdir(genome_path):
					subprocess.call('gunzip '+ genome_fna_zip , shell=True)
			
			'''### 5. IMMUNE SYSTEMS PREDICTION ###'''
			input_folder = genome_path
			output_folder = genome_path
			if True:
				os.chdir(d_analysis_from_ncbi_geno)
				script_5 = '5_RUN_immune_systems.sh'
				if s_running_in == 'server':
					script_5 = '5_RUN_immune_systems_server.sh'
				comand_5 = str('bash '+script_5+' '+input_folder+' '+output_folder)
				subprocess.call(comand_5,shell=True)

			'''### 6. GETTING ANNOTATED FILE NAMES ###'''
			file_fna , file_gff , file_faa = '' , '' , ''
			for files in os.listdir(genome_path):
				dict_id_location = {}
				if files.endswith('_PADLOC'):
					padloc_folder_path = os.path.join(genome_path , files)
					os.chdir(padloc_folder_path)
					for f in os.listdir(padloc_folder_path):
						if f.endswith('_genomic_prodigal.gff'):
							file_gff = os.path.join(padloc_folder_path,f)
						if f.endswith('_genomic_prodigal.faa'):
							file_faa = os.path.join(padloc_folder_path,f)
				if files.endswith('_genomic.fna'):
					file_fna = os.path.join(genome_path,files)

			'''### 7. GETTING LOCATION OF EACH LOCUS_TAG ###'''
			if True:
				os.chdir(d_analysis_from_ncbi_geno)
				comand_7 = str('bash 7_RUN_immune_systems_location.sh '+file_gff+' '+file_faa+' '+file_fna+' '+genome_path)
				subprocess.call(comand_7,shell=True)

			'''### 8. MAKING DICTIONARY OF THE LOCATION OF EACH LOCUS_TAG ###'''
			os.chdir(genome_path)
			dic_locations = {} # locus_tag contig start end strand
			for l in open('locations.txt').readlines():
				l = l.rstrip().split('\t')
				dic_locations[l[0]] = tuple(l[1:])

			'''### 9. RETRIEVING THE INFORMATION OF THE IMMUNE SYSTEMS PREDICTED ###'''
			padloc_defensefinder_systems = []
			padloc_systems = []
			defensefinder_systems = []
			for files in os.listdir(genome_path):
				if files.endswith('_genomic.Padloc'):
					# PADLOC FILE
					padloc_systems = get_padloc_systems(files)

				if files.endswith('_genomic_prodigal.DefenseFinder'):
					# DEFENSE FINDER FILE
					defensefinder_systems = get_defensefinder_systems(files)

			os.chdir(genome_path)
			padloc_defensefinder_systems = padloc_systems + defensefinder_systems
			immune_systems_locations_txt = os.path.join(genome_path,'immune_systems_locations.txt')
			with open('immune_systems_locations.txt','w') as fh:
				wfh = csv.writer(fh,delimiter='\t')
				wfh.writerows(padloc_defensefinder_systems)

			'''### 10. RETRIEVING THE SEQUENCES OF EACH IMMUNE SYSTEM (fna and faa) ###'''
			if True:
				os.chdir(d_analysis_from_ncbi_geno)
				comand_10 = str('bash 10_RUN_immune_systems_sequences.sh '+immune_systems_locations_txt+' '+file_faa+' '+file_fna+' '+genome_path)
				subprocess.call(comand_10,shell=True)

			os.chdir(genome_path)
			for final_files in os.listdir(genome_path):
				if final_files.endswith('.zip'):# genome_fna_zip in os.listdir(genome_path):
					os.remove(final_files)
				if final_files.endswith('locations.txt'):
					pass
					#os.remove(final_files)
				if final_files.endswith('_genomic.fna'):
					pass
					#os.remove(final_files)
				if final_files.endswith('_PADLOC'):
					shutil.rmtree(final_files)
	else:
		print('\n***',sp_genome,'***')
		print('genome filtered!')