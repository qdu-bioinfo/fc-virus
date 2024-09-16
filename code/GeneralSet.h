// This file was modified from a file named common.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpack Software LICENSE.



#ifndef GENERALSET_H
#define GENERALSET_H



#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>


using namespace std;

//k-mer
extern int g_kmer_length;
extern int g_min_kmer_coverage;
extern float g_min_kmer_entropy;
extern int g_min_seed_coverage;
extern float g_min_seed_entropy;
extern float g_min_ratio_non_error;
extern int g_ave_kmer_depth;
extern int g_homo_start;
extern int g_homo_end;

//reads
extern int g_mid_read_id;
extern bool g_is_paired_end; //false;
extern int g_fr_strand;
extern bool g_double_stranded_mode; //false;


//file
extern string g_reads_file;
extern string g_left_file;
extern string g_right_file;
extern string out_dir; 
extern string g_file_type;

//others
extern bool g_help;


#define OPT_DOUBLE_STRANDED_MODE		303
#define OPT_FR_STRAND		304
#define OPT_LEFT		308
#define OPT_RIGHT	309
#define OPT_SINGLEFILE			311



#endif
