
#include "KmerHash.h"
#include "HomoKmer.h"
#include "GeneralSet.h"
#include "Consensus.h"
#include "ReadUtility.h"
#include <iostream>
#include <ctime>
#include <errno.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
#include <getopt.h>


using namespace std;


struct option opts[] = {
	{"kmer_length",   required_argument,   0,   'k'},
	{"out_dir",        required_argument,   0,   'o'},
	{"pair reads",          no_argument,         0,   'p'},
	{"file_type",          required_argument,         0,   't'},
	{"help",          no_argument,         0,   'h'},
	{"double_stranded_mode", no_argument,  0,   OPT_DOUBLE_STRANDED_MODE},
	{"fr",     required_argument,   0,   OPT_FR_STRAND},
	{"left",          required_argument,   0,   OPT_LEFT},
	{"right",         required_argument,   0,   OPT_RIGHT},
	{"singlefile",    required_argument,   0,   OPT_SINGLEFILE},
	{0,0,0,0}

};


string usage() {

	stringstream usage_info;
	usage_info
		<< endl
		<< "===============================================================================" << endl
		<< " FC-Virus Usage " << endl
		<< "===============================================================================" << endl
		<< " ** Options: **" <<endl
		<< "  -k <int>: length of kmer, default 21. " << endl
		<< "  -o <string>: output directory. " << endl
		<< "  -p : paired-end reads. " << endl
		<< "  -t <string>: type of file, fa or fq. " << endl
		<< "  -h : help information. " << endl
		<< " If pair end reads: " << endl
		<< "  --left <string>: left reads file name (.fasta or .fastq). " << endl
		<< "  --right <string>: right reads file name (.fasta or .fastq). " << endl
		<< " If single end reads: " << endl
		<< "  -singlefile <string>: reads file name (.fasta or .fastq). " << endl
		//<< "  --double_stranded_mode: indicate the pair-end read is double stranded mode" << endl
		//<< "  --fr <int>: only used for pair-end reads. 1: --1--> <--2--  2: <--1-- --2-->  3: --1--> --2--> or <--1-- <--2--, default 1. " << endl
		<< "===============================================================================" << endl
		<< endl;

	return usage_info.str();

}


int parse_options(int argc, char* argv[]) {

	int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, "k:o:p:t:h", opts, &option_index);
		switch (next_option) {
		case -1:
			break;
		case 'k':
			g_kmer_length = atoi(optarg);
			break;
		case 'o':
			out_dir = optarg;
			break;
		case 'p':
			g_is_paired_end = true;
			break;
		case 'h':
			g_help = true;
			break;
		case 't':
			g_file_type = optarg;
			break;
		case OPT_DOUBLE_STRANDED_MODE:
			g_double_stranded_mode = true;
			break;
		case OPT_FR_STRAND:
			g_fr_strand = atoi(optarg);
			break;
		case OPT_LEFT:
			g_left_file = optarg;
			break;
		case OPT_RIGHT:
			g_right_file = optarg;
			break;
		case OPT_SINGLEFILE:
			g_reads_file = optarg;
            break;
		default:
			exit(1);
		}

	} while (next_option != -1);

	if (g_help) {
		cout << usage();
		exit (1);
	}



	if (g_kmer_length > 32) {
		cout << "Error: the kmer length should shorter than 32." << endl;
		exit(1);
	}

	if (g_reads_file.length() > 1)
		g_is_paired_end = false;

	if (g_fr_strand != 1 && g_fr_strand != 2 && g_fr_strand != 3) {
		cout << "Error: --fr can only be 1, 2 or 3" << endl;
		exit(1);
	}

	return 0;

}



int main(int argc, char* argv[]){

	map<string, int> data;

	time_t s_time = time(NULL);
	data.clear();
	int parse_ret = parse_options(argc,argv);
	if (parse_ret)
		return parse_ret;
	if (!g_is_paired_end) {
		load_reads(g_reads_file, data, false);
	} else {
		if (g_double_stranded_mode) {								
			load_reads(g_left_file, data, false);								
			load_reads(g_right_file, data, false);               
		} else {              
			if (g_fr_strand == 2) {//--1-->  <--2--				
				load_reads(g_left_file, data, false);				
				load_reads(g_right_file, data, true);
            }                	
			if (g_fr_strand == 1) {//<--1-- --2-->				
				load_reads(g_left_file, data, true);				
				load_reads(g_right_file, data, false);               	 
			}
			if (g_fr_strand == 3) {//--1--> --2-->  or <--1-- <--2--				
				load_reads(g_left_file, data, false);				
				load_reads(g_right_file, data, false);              	
			}               
		}
	}

	KmerHash kmer_hash(g_kmer_length);
	kmer_hash.get_hash(data);
	vector<int> kmer_cov_vec;
	cout << "Building consensus..." << endl;
	kmer_hash.get_homologous_kmers(data, kmer_cov_vec);
	HomoKmer homo_hash(g_kmer_length);
	vector<string> homo_reads;
	vector<int> homo_reads_cov;
	int homo_seed = homo_hash.get_hash(kmer_hash, data, homo_reads, homo_reads_cov, kmer_cov_vec);

	construct_consensus(kmer_hash, homo_hash, homo_reads, homo_reads_cov, homo_seed);
	time_t e_time = time(NULL);
	cout << "Success! (elapsed time: " << (e_time - s_time) << " s)" << endl;
	return 1;

}

