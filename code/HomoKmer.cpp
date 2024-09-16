#include "HomoKmer.h"
#include <algorithm>
#include <iostream>
#include <math.h>
#include <time.h>
#include <numeric>


using namespace std;



HomoKmer::HomoKmer (int kmer_length) {

	this -> kmer_length = kmer_length;
	

}


vector<size_t>& HomoKmer::operator[](kmer_int_type kmer) {

	return homo_kmer_hash[kmer];

}


int HomoKmer::get_hash(KmerHash& kmer_hash, map<string, int>& data, vector<string>& homo_reads, vector<int>& homo_reads_cov, vector<int>& kmer_cov_vec) {

	time_t start_time = time(NULL);
		
	map<string, int>::iterator it;
    int ave_read_len = 0;
    int read_num = 0;
    for (it = data.begin(); it != data.end(); ++it) {
        read_num++;
        const string& read = it -> first;
        ave_read_len = ave_read_len + read.length();
        if (read_num > 500)
            break;
    }
    ave_read_len = ave_read_len / read_num;

    int min_homo = (ave_read_len / g_kmer_length) * 0.3;
    if (min_homo < 1)
        min_homo = 1;
    
    int max_homo_kmer = 100000 / ((ave_read_len - min_homo*g_kmer_length)/min_homo + g_kmer_length);
 


    int idx = g_homo_end;
    if (idx >= kmer_cov_vec.size()) {
        idx = kmer_cov_vec.size()-1;
        g_homo_end = idx;
    }

    int homo_sum = 0;
    while(idx >= g_homo_start) {
        homo_sum = homo_sum + kmer_cov_vec[idx];
        //*
        if (homo_sum > max_homo_kmer) {
            g_homo_start = idx;
            break;
        }
        //*/
        idx--;
    }

/*
    if (homo_sum < max_homo_kmer*0.3) {
       
        while(g_homo_end < kmer_cov_vec.size() && homo_sum < max_homo_kmer*0.2) {
            g_homo_end++;
            homo_sum = homo_sum + kmer_cov_vec[g_homo_end];
        }

         while(g_homo_start > 2 && homo_sum < max_homo_kmer*0.2) {
            g_homo_start--;
            homo_sum = homo_sum + kmer_cov_vec[g_homo_start];
        }


    }

*/

    int max_homo_seed = -1;
    int temp_max_homo = 0;
	for (it = data.begin(); it != data.end(); ++it) {
		const string& read = it -> first;
		if (read.length() < kmer_length)
			continue;

        vector<kmer_int_type> temp_homo;
		for (int j = 0; j <= read.length()-kmer_length;++j) {
			const string& kmer = read.substr(j, kmer_length);
            
            int cng = contains_non_gatc(kmer);
			if (cng != -1) {
                j = j + cng;
                continue;
            }
            kmer_int_type kmer_int = kmer_to_int(kmer);
			if (kmer_hash[kmer_int] >= g_homo_start && kmer_hash[kmer_int] <= g_homo_end) {
                temp_homo.push_back(kmer_int);
            }
		}
        if (temp_homo.size() > 1) {
            for (int i = 0; i < temp_homo.size(); ++i) {
                homo_kmer_hash[temp_homo[i]].push_back(homo_reads.size());
            }
            if (temp_max_homo < temp_homo.size()) {
                temp_max_homo = temp_homo.size();
                max_homo_seed = homo_reads.size();
            }
            homo_reads.push_back(read);
            homo_reads_cov.push_back(it->second);
            
        }
	}


    return max_homo_seed;

}

