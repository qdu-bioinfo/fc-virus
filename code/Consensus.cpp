
#include "Consensus.h"
#include<iostream>

using namespace std;


string forward_extend_by_homo_reads(KmerHash& kmer_hash, HomoKmer& homo_hash, vector<string>& homo_reads, vector<int>& homo_reads_cov, kmer_int_type seed) {

	vector<kmer_int_type> homo_kmers;
	for (int i = 0; i <= homo_reads[seed].length()-g_kmer_length; ++i) {
		const string& kmer = homo_reads[seed].substr(i, g_kmer_length);
		kmer_int_type kmer_int = kmer_to_int(kmer);
		if(homo_hash[kmer_int].size() > 0) {
			homo_kmers.push_back(kmer_int);
		}
	}


	string contig = "";

	kmer_int_type homo1, homo2;
	string kmer1, kmer2;
	vector<size_t> homo_reads1;
	int check_size = 0;
	int homo_seed = -1;
	map<string, int>::iterator it;
	string max_str = "";
	int max_value = 0;
	int max_seed = 0;

	while(homo_kmers.size() > 1) {

		for (int i = 0; i < homo_kmers.size()-1; ++i) {

			homo1 = homo_kmers[i];
			homo2 = homo_kmers[i+1];
			kmer1 = kmer_to_string(homo1, g_kmer_length);
			kmer2 = kmer_to_string(homo2, g_kmer_length);		
			
			check_size = homo_reads1.size();
			if (check_size > 2000)
				check_size = 2000;

			homo_reads1 = homo_hash[homo1];
			map<string, int> candi_map;
			int temp_sum = 0;

			for (int j = 0; j < homo_reads1.size(); ++j) {
				
				string::size_type pos2 = homo_reads[homo_reads1[j]].find(kmer2);
				if (pos2 != string::npos) {
					string::size_type pos1 = homo_reads[homo_reads1[j]].find(kmer1);
					if (pos1 < pos2) {
						temp_sum++;
						const string& candi_str = homo_reads[homo_reads1[j]].substr(pos1, pos2-pos1);
						
						if (candi_map.find(candi_str) == candi_map.end()) {
							candi_map[candi_str] = homo_reads_cov[homo_reads1[j]];
						} else {
							candi_map[candi_str] = candi_map[candi_str] + homo_reads_cov[homo_reads1[j]];
						}
						if (i + 2 == homo_kmers.size() && max_seed < homo_reads_cov[homo_reads1[j]] && pos2 + g_kmer_length < homo_reads[homo_reads1[j]].length()) {
							for (int k = pos2+1; k <= homo_reads[homo_reads1[j]].length() - g_kmer_length; ++k) {
								const string& kmer = homo_reads[homo_reads1[j]].substr(k, g_kmer_length);
								kmer_int_type kmer_int = kmer_to_int(kmer);
								if (homo_hash[kmer_int].size() > 0 && homo_reads_cov[homo_reads1[j]] == kmer_hash[kmer_int]) {
									homo_seed = homo_reads1[j];
									max_seed = homo_reads_cov[homo_reads1[j]];
									break;
								}
							}
						}
						if (temp_sum > check_size)
							break;
					}
				}
			}

			max_value = 0;
			max_str = "";
			for (it = candi_map.begin(); it != candi_map.end(); ++it) {
				if (it -> second > max_value) {
					max_value = it -> second;
					max_str = it -> first;
				}
			}
			
			contig = contig + max_str;

		}

		////new candidates
		
		if (homo_seed == -1) {
			homo_reads1 = homo_hash[homo2];
			for (int j = 0; j < homo_reads1.size(); ++j) {
				string::size_type pos2 = homo_reads[homo_reads1[j]].find(kmer2);
				if (pos2 + g_kmer_length < homo_reads[homo_reads1[j]].length() && max_seed < homo_reads_cov[homo_reads1[j]]) {
					for (int k = pos2+1; k <= homo_reads[homo_reads1[j]].length() - g_kmer_length; ++k) {
						const string& kmer = homo_reads[homo_reads1[j]].substr(k, g_kmer_length);
						kmer_int_type kmer_int = kmer_to_int(kmer);
						if (homo_hash[kmer_int].size() > 0 && homo_reads_cov[homo_reads1[j]] == kmer_hash[kmer_int]) {
							homo_seed = homo_reads1[j];
							max_seed = homo_reads_cov[homo_reads1[j]];
							break;
						}
					}
				}
			}
			
		}
		
		
		if (homo_seed == -1)
			break;

		for (int i = 0; i < homo_kmers.size(); ++i) {
			homo_hash[homo_kmers[i]].clear();
		}

		homo_kmers.clear();
		string::size_type ss = homo_reads[homo_seed].find(kmer2);
		if (ss == string::npos)
			break;
		homo_kmers.push_back(homo2);
		for (int i = ss+1; i <= homo_reads[homo_seed].length()-g_kmer_length; ++i) {
			const string& kmer = homo_reads[homo_seed].substr(i, g_kmer_length);
			kmer_int_type kmer_int = kmer_to_int(kmer);
			if(homo_hash[kmer_int].size() > 0) {
				homo_kmers.push_back(kmer_int);
			}
		}

		homo_seed = -1;
		max_seed = 0;

	}


	contig = contig + kmer2;

	return contig;

}




string reverse_extend_by_homo_reads(KmerHash& kmer_hash, HomoKmer& homo_hash, vector<string>& homo_reads, vector<int>& homo_reads_cov, kmer_int_type seed, vector<kmer_int_type>& seed_kmers, vector<size_t>& reads1) {

	if (seed_kmers.size() < 2) 
		return "";

	kmer_int_type homo1 = seed_kmers[0];
	kmer_int_type homo2 = seed_kmers[1];
	string kmer1 = kmer_to_string(homo1, g_kmer_length);
	string kmer2 = kmer_to_string(homo2, g_kmer_length);
	vector<size_t> homo_reads1 = reads1;

	int max_value = 0;
	int max_seed = 0;
	int homo_seed = -1;
	string::size_type end_pos;
	

	for (int j = 0; j < homo_reads1.size(); ++j) {
				
		string::size_type pos2 = homo_reads[homo_reads1[j]].find(kmer2);
		if (pos2 != string::npos) {
			string::size_type pos1 = homo_reads[homo_reads1[j]].find(kmer1);
			if (pos1 < pos2 && pos1 != string::npos && max_seed < homo_reads_cov[homo_reads1[j]]) {
				if (pos1 > 0) {
					for (int k = 0; k < pos1; ++k) {
						const string& kmer = homo_reads[homo_reads1[j]].substr(k, g_kmer_length);
						kmer_int_type kmer_int = kmer_to_int(kmer);
						if (homo_hash[kmer_int].size() > 0 && kmer_hash[kmer_int] > max_seed) {
							homo_seed = homo_reads1[j];
							max_seed = homo_reads_cov[homo_reads1[j]];
							end_pos = pos1;
							break;
						}
					}
				}
			}
		}
	}
	
	if (homo_seed == -1) {
		for (int j = 0; j < homo_reads1.size(); ++j) {
			string::size_type pos1 = homo_reads[homo_reads1[j]].find(kmer1);
			if (pos1 > 0 && pos1 != string::npos && max_seed < homo_reads_cov[homo_reads1[j]]) {
				for (int k = 0; k < pos1; ++k) {
					const string& kmer = homo_reads[homo_reads1[j]].substr(k, g_kmer_length);
					kmer_int_type kmer_int = kmer_to_int(kmer);
					if (homo_hash[kmer_int].size() > 0 && homo_reads_cov[homo_reads1[j]] == kmer_hash[kmer_int]) {
						homo_seed = homo_reads1[j];
						max_seed = homo_reads_cov[homo_reads1[j]];
						end_pos = pos1;
						break;
					}
				}
			}	
		}
	}
	
	if (homo_seed == -1)
		return "";


	string contig = "";
	int check_size = 0;

	map<string, int>::iterator it;
	string max_str = "";


	vector<kmer_int_type> homo_kmers;
	homo_kmers.push_back(homo1);
	for (int i = end_pos - 1; i >= 0; --i) {
		const string& kmer = homo_reads[homo_seed].substr(i, g_kmer_length);
		kmer_int_type kmer_int = kmer_to_int(kmer);
		if(homo_hash[kmer_int].size() > 0) {
			homo_kmers.push_back(kmer_int);
		}
	}


	while(homo_kmers.size() > 1) {

		for (int i = 0; i < homo_kmers.size()-1; ++i) {

			homo1 = homo_kmers[i+1];
			homo2 = homo_kmers[i];
			kmer1 = kmer_to_string(homo1, g_kmer_length);
			kmer2 = kmer_to_string(homo2, g_kmer_length);		
			
			check_size = homo_reads1.size();
			if (check_size > 2000)
				check_size = 2000;

			homo_reads1 = homo_hash[homo1];
			map<string, int> candi_map;
			int temp_sum = 0;

			for (int j = 0; j < homo_reads1.size(); ++j) {
				
				string::size_type pos2 = homo_reads[homo_reads1[j]].find(kmer2);
				if (pos2 != string::npos) {
					string::size_type pos1 = homo_reads[homo_reads1[j]].find(kmer1);
					if (pos1 < pos2) {
						temp_sum++;
						const string& candi_str = homo_reads[homo_reads1[j]].substr(pos1, pos2-pos1);
						if (candi_map.find(candi_str) == candi_map.end()) {
							candi_map[candi_str] = homo_reads_cov[homo_reads1[j]];
						} else {
							candi_map[candi_str] = candi_map[candi_str] + homo_reads_cov[homo_reads1[j]];
						}
						/// new seed update
						if (i + 2 == homo_kmers.size() && pos1 > 0 && max_seed < homo_reads_cov[homo_reads1[j]]) {
							for (int k = 0; k < pos1; ++k) {
								const string& kmer = homo_reads[homo_reads1[j]].substr(k, g_kmer_length);
								kmer_int_type kmer_int = kmer_to_int(kmer);
								if (homo_hash[kmer_int].size() > 0 && homo_reads_cov[homo_reads1[j]] == kmer_hash[kmer_int]) {
									homo_seed = homo_reads1[j];
									max_seed = homo_reads_cov[homo_reads1[j]];
									break;
								}
							}
						}
						if (temp_sum > check_size)
							break;
					}
				}
			}

			max_value = 0;
			max_str = "";
			for (it = candi_map.begin(); it != candi_map.end(); ++it) {
				if (it -> second > max_value) {
					max_value = it -> second;
					max_str = it -> first;
				}
			}
			
			contig = max_str + contig;

		}

		////new candidates
		
		if (homo_seed == -1) {
			homo_reads1 = homo_hash[homo2];
			for (int j = 0; j < homo_reads1.size(); ++j) {
				string::size_type pos1 = homo_reads[homo_reads1[j]].find(kmer1);
				if (pos1 > 0 && pos1 != string::npos && max_seed < homo_reads_cov[homo_reads1[j]]) {
					for (int k = 0; k < pos1; ++k) {
						const string& kmer = homo_reads[homo_reads1[j]].substr(k, g_kmer_length);
						kmer_int_type kmer_int = kmer_to_int(kmer);
						if (homo_hash[kmer_int].size() > 0 && homo_reads_cov[homo_reads1[j]] == kmer_hash[kmer_int]) {
							homo_seed = homo_reads1[j];
							max_seed = homo_reads_cov[homo_reads1[j]];
							break;
						}
					}
				}
			}
			
		}
		
		//homo_seed = -1;
		//if (homo_seed == -1)
			break;

		for (int i = 0; i < homo_kmers.size(); ++i) {
			homo_hash[homo_kmers[i]].clear();
		}

		homo_kmers.clear();
		string::size_type ss = homo_reads[homo_seed].find(kmer1);
		if (ss == string::npos)
			break;
		homo_kmers.push_back(homo1);
		for (int i = ss-1; i >= 0; --i) {
			const string& kmer = homo_reads[homo_seed].substr(i, g_kmer_length);
			kmer_int_type kmer_int = kmer_to_int(kmer);
			if(homo_hash[kmer_int].size() > 0) {
				homo_kmers.push_back(kmer_int);
			}
		}

		homo_seed = -1;
		max_seed = 0;

	}

	return contig;

}





float coverage_estimation(KmerHash& kmer_hash, kmer_int_type seed) {

	map<kmer_int_type, bool> used_kmers;
	bool is_candi = true;
	vector<int> cov_vec;
	vector<kmer_pair> candidates;
	kmer_int_type kmer_int = seed;

	while (is_candi && cov_vec.size() < 2000) {
		kmer_hash.get_reverse_candidates(kmer_int, candidates);
		if (candidates.empty())
			break;
		is_candi = false;
		for (int i = 0; i < candidates.size(); ++i) {
			if (!used_kmers[candidates[i].first]) {
				is_candi = true;
				kmer_int = candidates[i].first;
				int total_cov = 0;
				for (int j = i; j < candidates.size(); ++j) {
					used_kmers[candidates[j].first] = true;
					total_cov = total_cov + candidates[j].second;
				}
				if (candidates.size() == 1) 
					cov_vec.push_back(total_cov);
				cov_vec.push_back(total_cov);
				break;
			}
		}
	}
	used_kmers.clear();

	is_candi = true;
	while(is_candi && cov_vec.size() < 2000) {

		kmer_hash.get_forward_candidates(kmer_int, candidates);
		if (candidates.empty())
			break;
		is_candi = false;
		for (int i = 0; i < candidates.size(); ++i) {
			if (!used_kmers[candidates[i].first]) {
				is_candi = true;
				kmer_int = candidates[i].first;
				int total_cov = 0;
				for (int j = i; j < candidates.size(); ++j) {
					used_kmers[candidates[j].first] = true;
					total_cov = total_cov + candidates[j].second;
				}
				if (candidates.size() == 1)
					cov_vec.push_back(total_cov);
				break;
			}
		}
	}

	sort(cov_vec.begin(), cov_vec.end());

	int idx = cov_vec.size()*0.1;

	float ave_cov = 0.0;
	int check_num = 0;
	for (int i = idx; i < cov_vec.size()-idx; ++i) {
		check_num++;
		ave_cov = ave_cov + cov_vec[i];
	}

	if (check_num > 1)
		ave_cov = ave_cov * 1.0 / check_num;
	
	return ave_cov;
}

string reverse_extend_by_kmers(KmerHash& kmer_hash, kmer_int_type seed, map<kmer_int_type, bool>& used_kmers, float& ave_cov) {


	kmer_int_type kmer_int = seed;
	string contig = kmer_to_string(kmer_int, g_kmer_length);
	vector<kmer_pair> candidates;
	int base_num, total_cov;
	char base;
	bool is_candi = true;

	while(is_candi) {
		kmer_hash.get_reverse_candidates(kmer_int, candidates);
		if (candidates.empty())
			break;
		is_candi = false;
		for (int i = 0; i < candidates.size(); ++i) {
			if (!used_kmers[candidates[i].first]) {
				is_candi = true;
				kmer_int = candidates[i].first;
				total_cov = candidates[i].second;
				for (int j = candidates.size()-1; j > i; --j) {
					if (total_cov + candidates[j].second > ave_cov)
						break;
					total_cov = total_cov + candidates[j].second;
					kmer_hash[candidates[j].first] = 0;
				}
				/*
				if (i + 1 < candidates.size()) {
					kmer_hash[candidates[candidates.size()-1].first] = 0;
					total_cov = total_cov + candidates[candidates.size()-1].second;
				}
				*/
				kmer_hash[kmer_int] = 0;
				used_kmers[kmer_int] = true;
				base_num = (kmer_int >> (g_kmer_length*2-2)) & 3ll;
				base = int_to_base(base_num);
				contig = base + contig;
				break;
			}
		}
	}

	return contig;
}



string forward_extend_by_kmers(KmerHash& kmer_hash, kmer_int_type seed, map<kmer_int_type, bool>& used_kmers, float& ave_cov) {

	kmer_int_type kmer_int = seed;
	string contig = kmer_to_string(kmer_int, g_kmer_length);
	vector<kmer_pair> candidates;
	int base_num, total_cov;
	char base;
	bool is_candi = true;
	while(is_candi) {

		kmer_hash.get_forward_candidates(kmer_int, candidates);
		if (candidates.empty())
			break;
		is_candi = false;
		for (int i = 0; i < candidates.size(); ++i) {
			if (!used_kmers[candidates[i].first]) {
				is_candi = true;
				kmer_int = candidates[i].first;
				total_cov = candidates[i].second;
				for (int j = candidates.size()-1; j > i; --j) {
					if (total_cov + candidates[j].second > ave_cov)
						break;
					total_cov = total_cov + candidates[j].second;
					kmer_hash[candidates[j].first] = 0;
				}
				/*
				if (i + 1 < candidates.size()) {
					kmer_hash[candidates[candidates.size()-1].first] = 0;
					total_cov = total_cov + candidates[candidates.size()-1].second;
				}
				*/
				used_kmers[kmer_int] = true;
				kmer_hash[kmer_int] = 0;
				base_num = kmer_int & 3ll;
				base = int_to_base(base_num);
				contig = contig + base;
				break;
			}
		}
	}

	return contig;
}




float compute_consensus_depth(KmerHash& kmer_hash, string& contig) {

	int total_cov;
	vector<int> cov_vec;
	map<kmer_int_type, bool> used_kmers;
	for (int i = 0; i <= contig.length()-g_kmer_length; ++i) {
		total_cov = 0;
		const string& kmer = contig.substr(i, g_kmer_length);
		kmer_int_type kmer_int = kmer_to_int(kmer);
		used_kmers[kmer_int] = true;
		cov_vec.push_back(kmer_hash[kmer_int]);
	}



	bool is_candi = true;
	vector<kmer_pair> candidates;
	string kmer = contig.substr(0, g_kmer_length);
	kmer_int_type kmer_int = kmer_to_int(kmer);

	while (is_candi && cov_vec.size() < 2000) {
		kmer_hash.get_reverse_candidates(kmer_int, candidates);
		if (candidates.empty())
			break;
		is_candi = false;
		for (int i = 0; i < candidates.size(); ++i) {
			if (!used_kmers[candidates[i].first]) {
				is_candi = true;
				kmer_int = candidates[i].first;
				int total_cov = 0;
				for (int j = i; j < candidates.size(); ++j) {
					used_kmers[candidates[j].first] = true;
					total_cov = total_cov + candidates[j].second;
				}
				if (candidates.size() == 1) 
					cov_vec.push_back(total_cov);
				cov_vec.push_back(total_cov);
				break;
			}
		}
	}
	//used_kmers.clear();

	is_candi = true;
	kmer = contig.substr(contig.length()-g_kmer_length);
	kmer_int = kmer_to_int(kmer);

	while(is_candi && cov_vec.size() < 2000) {

		kmer_hash.get_forward_candidates(kmer_int, candidates);
		if (candidates.empty())
			break;
		is_candi = false;
		for (int i = 0; i < candidates.size(); ++i) {
			if (!used_kmers[candidates[i].first]) {
				is_candi = true;
				kmer_int = candidates[i].first;
				int total_cov = 0;
				for (int j = i; j < candidates.size(); ++j) {
					used_kmers[candidates[j].first] = true;
					total_cov = total_cov + candidates[j].second;
				}
				if (candidates.size() == 1)
					cov_vec.push_back(total_cov);
				break;
			}
		}
	}

	sort(cov_vec.begin(), cov_vec.end());
	int idx = cov_vec.size()*0.1;

	float ave_cov = 0.0;
	int check_num = 0;
	for (int i = idx; i < cov_vec.size()-idx; ++i) {
		check_num++;
		ave_cov = ave_cov + cov_vec[i];
	}

	if (check_num > 1)
		ave_cov = ave_cov * 1.0 / check_num;
	

	return ave_cov;


}


void update_kmers_status(KmerHash& kmer_hash, string& contig, float ave_cov) {

	vector<kmer_pair> candidates;
	string kmer = contig.substr(0, g_kmer_length);
	kmer_int_type kmer_int = kmer_to_int(kmer);
	float temp_cov = 0;

	if (kmer_hash[kmer_int] > 2* ave_cov) {
			kmer_hash[kmer_int] = kmer_hash[kmer_int] - ave_cov;
		} else {
			kmer_hash[kmer_int] = 0;
		}

	for (int i = 0; i < contig.length()-g_kmer_length; ++i) {
		kmer = contig.substr(i, g_kmer_length);
		kmer_int = kmer_to_int(kmer);
		kmer_hash.get_forward_candidates(kmer_int, candidates);
		kmer = contig.substr(i+1, g_kmer_length);
		kmer_int = kmer_to_int(kmer);
		temp_cov = kmer_hash[kmer_int];
		float used_cov = 0;
		if (candidates.size() > 0) {
			
			if (temp_cov <= candidates[0].second) {
		
				kmer = kmer_to_string(candidates[0].first, g_kmer_length);
				contig = contig.substr(0, i) + kmer + contig.substr(i+g_kmer_length);
				if (candidates[0].second > 2*ave_cov) {
					kmer_hash[candidates[0].first] = kmer_hash[candidates[0].first] - ave_cov;
					used_cov = ave_cov;
				} else {
					kmer_hash[candidates[0].first] = 0;
					used_cov = candidates[0].second;
				}
				
				
			} else {
				if (kmer_hash[kmer_int] > 2* ave_cov) {
					kmer_hash[kmer_int] = kmer_hash[kmer_int] - ave_cov;
					used_cov = ave_cov;
				} else {
					kmer_hash[kmer_int] = 0;
					used_cov = kmer_hash[kmer_int];
				}
					
				
			}
			
			int idx = candidates.size()-1;
			while(idx >= 0 && used_cov < ave_cov) {
				if (candidates[idx].second + used_cov > 2*ave_cov) {
					kmer_hash[candidates[idx].first] = kmer_hash[candidates[idx].first] + used_cov - ave_cov;
				} else {
					kmer_hash[candidates[idx].first] = 0;
				}
				used_cov = used_cov + candidates[idx].second;
				idx--;
			}	
		}
		
	}



}



void construct_consensus(KmerHash& kmer_hash, HomoKmer& homo_hash, vector<string>& homo_reads, vector<int>& homo_reads_cov, int seed) {


	vector<kmer_int_type> homo_kmers;
	for (int i = 0; i <= homo_reads[seed].length()-g_kmer_length; ++i) {
		const string& kmer = homo_reads[seed].substr(i, g_kmer_length);
		kmer_int_type kmer_int = kmer_to_int(kmer);
		if(homo_hash[kmer_int].size() > 0) {
			homo_kmers.push_back(kmer_int);
			if (homo_kmers.size() == 2)
				break;
		}
	}

	vector<size_t> reads1 = homo_hash[homo_kmers[0]];
	string right = forward_extend_by_homo_reads(kmer_hash, homo_hash, homo_reads, homo_reads_cov, seed);
	string left = reverse_extend_by_homo_reads(kmer_hash, homo_hash, homo_reads, homo_reads_cov, seed, homo_kmers, reads1);
	string trunk = left + right;
	//string trunk = right;

	float ave_cov = compute_consensus_depth(kmer_hash, trunk);
	//update_kmers_status(kmer_hash, trunk, ave_cov);


	map<kmer_int_type, bool> used_kmers;
	kmer_int_type kmer_int;


	string kmer = trunk.substr(0, g_kmer_length);
	kmer_int = kmer_to_int(kmer);
	left = reverse_extend_by_kmers(kmer_hash, kmer_int, used_kmers, ave_cov);
	kmer = trunk.substr(trunk.length()-g_kmer_length);
	kmer_int = kmer_to_int(kmer);
	right = forward_extend_by_kmers(kmer_hash, kmer_int, used_kmers, ave_cov);
	//trunk = left + right.substr(g_kmer_length);
	trunk = left + trunk.substr(g_kmer_length) + right.substr(g_kmer_length);

	cout << "The length of reconstructed consensus is " << trunk.length() << " bp." << endl;


	string file_name = out_dir + "FC-Virus.fa";
	fstream trunk_file;
	trunk_file.open(file_name.c_str(), fstream::out);
	trunk_file << ">FC-Virus_Consensus_" << trunk.length() << endl;
	trunk_file << trunk << endl;
	trunk_file.close();

}

