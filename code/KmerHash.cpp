#include "KmerHash.h"
#include <algorithm>
#include <iostream>
#include <math.h>
#include <time.h>
#include <numeric>


using namespace std;


bool cmp_by_value(const pair<int,int>& a, const pair<int, int>& b) {
	return a.second > b.second;
}

//functions in class KmerHash
KmerHash::KmerHash (int kmer_length) {

	this -> kmer_length = kmer_length;
	

}


int& KmerHash::operator[](kmer_int_type kmer) {

	return kmer_hash[kmer];

}


size_t KmerHash::get_size() {

	return kmer_hash.size();

}


bool KmerHash::empty() {

	return kmer_hash.empty();

}


bool KmerHash::reuse(const kmer_int_type kmer) {

	kmer_hash_type_iterator it = find_kmer(kmer);
	return (it != kmer_hash.end());
	
}


bool KmerHash::exists(kmer_int_type kmer) {

	return (kmer_abundance(kmer) > 0);
	
}


bool KmerHash::exists(const string& kmer) {

	kmer_int_type it = kmer_to_int(kmer);
	return (exists(it));
	
}


size_t KmerHash::kmer_abundance(kmer_int_type kmer) {

	kmer_hash_type_iterator it = find_kmer(kmer);

	if (it != kmer_hash.end())
		return it -> second;
	else
		return 0;

}


size_t KmerHash::kmer_abundance(const string& kmer) {

	kmer_int_type it = kmer_to_int(kmer);
	return (kmer_abundance(it));

}


kmer_int_type KmerHash::get_seed_kmer() {

	kmer_int_type max_kmer = kmer_hash.begin()->first;
	int max_cov = 0;
	
	kmer_hash_type_iterator it;

	kmer_int_type temp_kmer;
	vector<kmer_pair> candidates;
	for (it = kmer_hash.begin(); it != kmer_hash.end(); ++it) {
		if (static_cast<int>(it->second) > max_cov) {
			temp_kmer = it -> first;
			get_reverse_candidates(temp_kmer, candidates);
			if (candidates.size() > 1)
				continue;
			get_forward_candidates(temp_kmer, candidates);
			if (candidates.size() > 1)
				continue;
			max_cov = static_cast<int>(it->second);
			max_kmer = it -> first;
		}
	}

	return max_kmer;
}





void KmerHash::get_reverse_candidates(kmer_int_type seed_kmer, vector<kmer_pair>& candidates) {

	candidates.clear();
	kmer_int_type reverse_suffix = seed_kmer >> 2;

	for (kmer_int_type i = 0; i < 4; ++i) {
		kmer_pair candidate;
		candidate.first = (i << (kmer_length*2-2)) | reverse_suffix;
		candidate.second = kmer_abundance(candidate.first);
		if (candidate.second) {
			candidates.push_back(candidate);
		} 
	}

	kmer_sorter_by_count_desc_t sorter(*this);
	sort(candidates.begin(), candidates.end(), sorter);

}

void KmerHash::get_forward_candidates(kmer_int_type seed_kmer, vector<kmer_pair>& candidates) {
	
	candidates.clear();
	kmer_int_type forward_prefix = (seed_kmer << (33-kmer_length)*2) >> (32-kmer_length)*2;
	for (kmer_int_type i = 0; i < 4; ++i) {
		kmer_pair candidate;
		candidate.first = forward_prefix | i;
		candidate.second = kmer_abundance(candidate.first);
		if (candidate.second) {
			candidates.push_back(candidate);
		} 
	}

	kmer_sorter_by_count_desc_t sorter(*this);
	sort(candidates.begin(), candidates.end(), sorter);

}



void KmerHash::get_forward_candidates(vector<kmer_int_type>& seeds, vector<kmer_pair>& candidates) {

	candidates.clear();
	for (int j = 0; j < seeds.size(); ++j) {
		kmer_int_type seed_kmer = seeds[j];
		kmer_int_type forward_prefix = (seed_kmer << (33-kmer_length)*2) >> (32-kmer_length)*2;
		for (kmer_int_type i = 0; i < 4; ++i) {
			kmer_pair candidate;
			candidate.first = forward_prefix | i;
			candidate.second = kmer_abundance(candidate.first);
			if (candidate.second) {
				candidates.push_back(candidate);
			} 
		}
	}

	//kmer_sorter_by_count_desc_t sorter(*this);
	//sort(candidates.begin(), candidates.end(), sorter);

}



void KmerHash::get_reverse_candidates(vector<kmer_int_type>& seeds, vector<kmer_pair>& candidates) {

	candidates.clear();

	for (int j = 0; j < seeds.size(); ++j) {

		kmer_int_type seed_kmer = seeds[j];
		kmer_int_type reverse_suffix = seed_kmer >> 2;
		for (kmer_int_type i = 0; i < 4; ++i) {
			kmer_pair candidate;
			candidate.first = (i << (kmer_length*2-2)) | reverse_suffix;
			candidate.second = kmer_abundance(candidate.first);
			if (candidate.second) {
				candidates.push_back(candidate);
			} 
		}

	}
	
	//kmer_sorter_by_count_desc_t sorter(*this);
	//sort(candidates.begin(), candidates.end(), sorter);

}


void KmerHash::get_hash(map<string, int>& data) {


	size_t data_size = data.size();


	if (data_size == 0) {
		cout << "Error, no reads exsit!" << endl;
		return;
	}

	cout << "constructing kmer hash ..." << endl;

	time_t start_time = time(NULL);
		
	map<string, int>::iterator it;

	for (it = data.begin(); it != data.end(); ++it) {
		const string& read = it -> first;
		if (read.length() < kmer_length)
			continue;
		for (int j = 0; j <= read.length()-kmer_length;++j) {
			const string& kmer = read.substr(j, kmer_length);
			int cng = contains_non_gatc(kmer);
			if (cng != -1) {
                j = j + cng;
                continue;
            }
			kmer_int_type kmer_int = kmer_to_int(kmer);

			if (exists(kmer_int)) {
				kmer_hash[kmer_int] = kmer_hash[kmer_int] + it -> second;
				continue;
			}
			kmer_hash[kmer_int] = it -> second;
		}
	}
	time_t end_time = time(NULL);
	cout << "Kmer Hash has been constructed, total " << kmer_hash.size() << " kmers! (elapsed time: " << (end_time - start_time) << " s)" << endl;


}


void KmerHash::get_homologous_kmers(map<string, int>& data, vector<int>& kmer_count) {

	int max_cov = 0;
	kmer_hash_type_iterator it;
	for (it = kmer_hash.begin(); it != kmer_hash.end(); ++it) {
		if (static_cast<int>(it->second) > max_cov) {
			max_cov = it -> second;
		}
	}

	//const int x = max_cov + 1;
	//vector<int> kmer_count(x, 0);
	for (int i = 0; i < max_cov+1; ++i) {
		kmer_count.push_back(0);
	}

	for (it = kmer_hash.begin(); it != kmer_hash.end(); ++it) {
		int id = it -> second;
		kmer_count[id] = kmer_count[id] + 1;
	}


	int large_num = 0;
	int large_sum = 0;

	for (int i = 2; i <= max_cov; ++i) {
		if (kmer_count[i] > 0) {
			large_num = large_num + kmer_count[i];
			large_sum = large_sum + kmer_count[i]*i;
		}
		
	}
	
	g_ave_kmer_depth = large_sum / large_num;

	large_num = 0;
	large_sum = 0;

	for (int i = g_ave_kmer_depth; i <= max_cov; ++i) {
		if (kmer_count[i] > 0) {
			large_num = large_num + 1;
		}	
	}


	int window_len = 3 * (max_cov - g_ave_kmer_depth + 1) / large_num; 
	//if (window_len < 5)
	//	window_len = 5;
	//int step_len = window_len*0.5;
	//if (step_len < 1)
	//	step_len = 1;

	vector<int> window_vec;
	int idx = g_ave_kmer_depth;
	int temp = 0;
	while(idx < kmer_count.size()) {
		temp = 0;
		int i = 0;
		int j = 0;
		while(i < window_len && idx + i < kmer_count.size()) {
			temp = temp + kmer_count[idx + i];
			if (kmer_count[idx+i] > 0)
				j++;
			i++;
		}
		if (j > 0)
			temp = temp / j;
		window_vec.push_back(temp);
		idx = idx + window_len;
	}




	large_sum = 0;
	for (int i = 0; i < window_vec.size(); ++i)
		large_sum = large_sum + window_vec[i];

	int ave_value = large_sum / window_vec.size();


	idx = window_vec.size()-1;

	while (idx * 2 > window_vec.size()) {

		while(window_vec[idx] < ave_value && idx * 2 > window_vec.size()) {
			idx--;
		}
		if (window_vec[idx] < ave_value) {
			for (int i = idx; i < window_vec.size(); ++i) {
				if (window_vec[i] > window_vec[idx])
					idx = i;
			}
		}
		int peak_end = idx;
		while(peak_end + 1 < window_vec.size() && window_vec[peak_end] >= window_vec[peak_end+1]) {
			peak_end++;
		}
		
		if (peak_end + 2 < window_vec.size() && window_vec[peak_end+1] >= window_vec[peak_end+2]) {
			peak_end++;
			while(peak_end + 1 < window_vec.size() && window_vec[peak_end] >= window_vec[peak_end+1]) {
				peak_end++;
			}
		}
		
		int peak_start = idx;
		while(peak_start > 0 && window_vec[peak_start-1] >= window_vec[peak_start]) {
			peak_start--;
		}
		
		if (peak_start > 1  && window_vec[peak_start-2] >= window_vec[peak_start-1] ) {
			peak_start--;
			while(peak_start > 0 && window_vec[peak_start-1] >= window_vec[peak_start]) {
				peak_start--;
			}
		}
		
		while(peak_start > 0 && window_vec[peak_start-1] <= window_vec[peak_start]) {
			peak_start--;
		}
		
		if (peak_start > 1 && window_vec[peak_start-2] <= window_vec[peak_start-1]) {
			peak_start--;
			while(peak_start > 0 && window_vec[peak_start-1] <= window_vec[peak_start]) {
				peak_start--;
			}
		}
		
		g_homo_start = peak_start*window_len + g_ave_kmer_depth;
		g_homo_end = peak_end*window_len + window_len + g_ave_kmer_depth;

		if (check_homologous(data))
			break;
		else {
			g_homo_end = 0;
		}
		idx = peak_start;
	}

	if (g_homo_end == 0) {
		g_homo_end = max_cov;
		g_homo_start = max_cov*0.7;
	}
	
}


bool KmerHash::check_homologous(map<string, int>& data) {

	int homo_reads_sum = 0;
	map<string, int>::iterator dr;
	for (dr = data.begin(); dr != data.end(); ++dr) {
		const string& read = dr -> first;
		if (read.length() < kmer_length)
			continue;
		int is_homo = 0;
		for (int j = 0; j <= read.length()-kmer_length;++j) {
			const string& kmer = read.substr(j, kmer_length);
			int cng = contains_non_gatc(kmer);
			if (cng != -1) {
                j = j + cng;
                continue;
            }
			kmer_int_type kmer_int = kmer_to_int(kmer);
			if (kmer_hash[kmer_int] >= g_homo_start && kmer_hash[kmer_int] <= g_homo_end) {
				is_homo++;
				if (is_homo > 1) {
					homo_reads_sum++;
					break;
				}
			} 
		}
	}

	if (homo_reads_sum  < 2)
		return false;

	
	map<kmer_int_type, bool> homo_kmers;
	kmer_hash_type_iterator it;
	for (it = kmer_hash.begin(); it != kmer_hash.end(); ++it) {
		if (it->second >= g_homo_start && it -> second <= g_homo_end) {
			homo_kmers[it -> first] = false;
		}
	}

	map<kmer_int_type, bool>::iterator it2;
	vector<kmer_int_type> candidates;
	int max_graph_connection = 0;
	for (it2 = homo_kmers.begin(); it2 != homo_kmers.end(); ++it2) {
		if (it2 -> second)
			continue;
		int graph_connection = 1;
		vector<kmer_int_type> connected_kmers;
		vector<kmer_int_type> candidates;
		kmer_int_type seed = it2 -> first;
		get_forward_homo_candidates(seed, candidates);
		for (int i = 0; i < candidates.size(); ++i) {
			if (homo_kmers[candidates[i]] == false) {
				homo_kmers[candidates[i]] = true;
				graph_connection++;
				connected_kmers.push_back(candidates[i]);
			}
		}
		get_reverse_homo_candidates(seed, candidates);
		for (int i = 0; i < candidates.size(); ++i) {
			if (homo_kmers[candidates[i]] == false) {
				homo_kmers[candidates[i]] = true;
				graph_connection++;
				connected_kmers.push_back(candidates[i]);
			}
		}
		while(connected_kmers.size() > 0) {
			seed = connected_kmers.back();
			connected_kmers.pop_back();
			get_forward_homo_candidates(seed, candidates);
			for (int i = 0; i < candidates.size(); ++i) {
				if (homo_kmers[candidates[i]] == false) {
					homo_kmers[candidates[i]] = true;
					graph_connection++;
					connected_kmers.push_back(candidates[i]);
				}
			}
			get_reverse_homo_candidates(seed, candidates);
			for (int i = 0; i < candidates.size(); ++i) {
				if (homo_kmers[candidates[i]] == false) {
					homo_kmers[candidates[i]] = true;
					graph_connection++;
					connected_kmers.push_back(candidates[i]);
				}
			}

			if (graph_connection > max_graph_connection)
				max_graph_connection = graph_connection;

			if (graph_connection > homo_kmers.size()*0.8 && graph_connection < 2000) {
				return false;
			}
			
		}

	}

	return true;

} 



void KmerHash::get_reverse_homo_candidates(kmer_int_type seed_kmer, vector<kmer_int_type>& candidates) {

	candidates.clear();
	kmer_int_type reverse_suffix = seed_kmer >> 2;

	for (kmer_int_type i = 0; i < 4; ++i) {
		kmer_pair candidate;
		candidate.first = (i << (kmer_length*2-2)) | reverse_suffix;
		candidate.second = kmer_abundance(candidate.first);
		if (candidate.second >= g_homo_start && candidate.second <= g_homo_end) {
			candidates.push_back(candidate.first);
		} 
	}

}

void KmerHash::get_forward_homo_candidates(kmer_int_type seed_kmer, vector<kmer_int_type>& candidates) {
	
	candidates.clear();
	kmer_int_type forward_prefix = (seed_kmer << (33-kmer_length)*2) >> (32-kmer_length)*2;
	for (kmer_int_type i = 0; i < 4; ++i) {
		kmer_pair candidate;
		candidate.first = forward_prefix | i;
		candidate.second = kmer_abundance(candidate.first);
		if (candidate.second >= g_homo_start && candidate.second <= g_homo_end) {
			candidates.push_back(candidate.first);
		} 
	}

}

