# PROST-SEQ
A tool for predicting the effects of missense mutations on protein stability changes upon missense mutation using protein sequence only.

Requirements: Listed separately as requirementsPy2.txt and requirementsPy3.txt for virtual environments.

Installation of Anaconda3 is prefered

 1) create python3 virtual environment and fulfil (install packages) requirements3.txt

 2) create python2 virtual environment and fulfil (install packages) requirementspy2.txt [Required for running run_list_spd33.sh on a new sequence file]. Activate python 2  in spd33_run_list.sh
[Required for running on a new sequence] 

Download the following databases and unzip

 i) unire50 (https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50)

 ii) uniclust30_2018_08 (http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz)

 iii) uniprot20_2016_02 (https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_02/uniref/uniref2016_02.tar.gz)

2.1 Check and rectify paths to DATABASES in Mutation_pred.py and spd33_run_list.sh

Activate python3 virtual environment and run the python script (Mutation_pred.py):

Command-line arguments:

	{-file,--file}	protein sequence (FASTA format)

	{-mutation, --mutation}	missence mutation (example: A 12 W or GLN 10 ALA)

	{-mutlist, --mutlist, --ml, --mutation-list}	list of mutations

	{-outdir, --outdir, --out-dir}	directory name for results

	{-out-file, --out-file} Name for the result output file

	{-h, --help}	command-line summary
single mutation

python Mutation_pred.py -file fasta.txt -mutation wild-residue position mutant-residue  -outdir(optional) Result -out-file (optional) mutation_result
list of mutations

python Mutation_pred.py -file fasta.txt -mutlist Mut_list.txt -outdir(optional) Result -out-file(optional) mut_list_Result
Example:

	1) python Mutation_pred.py -file Input/Frataxin.txt -mutlist Input/Frataxin_mut.txt -outdir Result
	
	2) python Mutation_pred.py -file Input/Frataxin.txt -mutation D 21 G  -outdir Result -out-file D21G_result
**Internal files will be stored inside Gen_Files folder.
