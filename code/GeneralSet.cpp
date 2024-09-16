// This file was modified from a file named common.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpack Software LICENSE.



#include "GeneralSet.h"
#include <cstring>
#include <cstdlib>
#include <errno.h>

using namespace std;

//k-mer
int g_kmer_length = 25;
int g_min_kmer_coverage = 2;
float g_min_kmer_entropy = 0.0f;
int g_min_seed_coverage = 1;
float g_min_seed_entropy = 1.5f;
float g_min_ratio_non_error = 0.04f;
int g_ave_kmer_depth = 0;
int g_homo_start = 2;
int g_homo_end = 2;

//reads
int g_mid_read_id = 0;
bool g_is_paired_end = true; 
int g_fr_strand = 1;
bool g_double_stranded_mode = false; 



//trunk
int g_min_trunk_length = 200;

//file
string g_reads_file = "";
string g_left_file = "";
string g_right_file = "";
string out_dir = ""; 
string g_file_type = "";
//others
bool g_help = false;


