
#ifndef CONSENSUS_H
#define CONSENSUS_H

#include"KmerUtility.h"
#include"KmerHash.h"
#include"HomoKmer.h"
#include<algorithm>
#include<assert.h>
#include<fstream>
#include<map>
#include<numeric>
#include<iomanip>
#include<set>
#include<string>
#include<stdlib.h>
#include<sstream>
#include<vector>
#include<list>
#include <bits/stdc++.h>



using namespace std;


string forward_extend_by_homo_reads(KmerHash& kmer_hash, HomoKmer& homo_hash, vector<string>& homo_reads, vector<int>& homo_reads_cov, kmer_int_type seed);
string reverse_extend_by_homo_reads(KmerHash& kmer_hash, HomoKmer& homo_hash, vector<string>& homo_reads, vector<int>& homo_reads_cov, kmer_int_type seed, vector<kmer_int_type>& seed_kmers, vector<size_t>& reads1);
float coverage_estimation(KmerHash& kmer_hash, kmer_int_type seed);
string forward_extend_by_kmers(KmerHash& kmer_hash, kmer_int_type seed, map<kmer_int_type, bool>& used_kmers, float& ave_cov);
string reverse_extend_by_kmers(KmerHash& kmer_hash, kmer_int_type seed, map<kmer_int_type, bool>& used_kmers, float& ave_cov);	
void construct_consensus(KmerHash& kmer_hash, HomoKmer& homo_hash, vector<string>& homo_reads, vector<int>& homo_reads_cov, int seed);
float compute_consensus_depth(KmerHash& kmer_hash, string& contig);
void update_kmers_status(KmerHash& kmer_hash, string& contig, float ave_cov);

#endif
