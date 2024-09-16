#ifndef HOMOKMER_H
#define HOMOKMER_H


#include "KmerUtility.h"
#include "ReadUtility.h"
#include "KmerHash.h"
#include "GeneralSet.h"
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <boost/unordered_map.hpp>


using namespace std;


class HomoKmer {

private:

	typedef  boost::unordered_map<kmer_int_type, vector<size_t> > homo_kmer_type;
	typedef  boost::unordered_map<kmer_int_type, vector<size_t> >::iterator homo_kmer_type_iterator;
	typedef  boost::unordered_map<kmer_int_type, vector<size_t> >::const_iterator homo_kmer_type_const_iterator;

	homo_kmer_type homo_kmer_hash;
	int kmer_length;

public:

	HomoKmer () { }
	HomoKmer (int kmer_length);
	vector<size_t> & operator[](kmer_int_type kmer);

	homo_kmer_type_iterator find_kmer(kmer_int_type kmer){
		//if (g_double_stranded_mode)
			//kmer = get_DS_kmer_val(kmer, kmer_length);
		return homo_kmer_hash.find(kmer);
	}

	int get_hash(KmerHash& kmer_hash, map<string, int>& data, vector<string>& homo_reads, vector<int>& homo_reads_cov, vector<int>& kmer_cov_vec);

};



#endif
