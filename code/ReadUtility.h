#ifndef READUTILITY_H
#define READUTILITY_H


#include "GeneralSet.h"
#include"KmerUtility.h"
#include <vector>
#include <string>


using namespace std;


typedef vector<unsigned long long> read_int_type;

void load_reads(string file, map<string, int>& data, bool rev);
void load_reads_fa(string file, map<string, int>& data, bool rev);
void load_reads_fq(string file, map<string, int>& data, bool rev);

#endif
