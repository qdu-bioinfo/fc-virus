#include "ReadUtility.h"
#include "GeneralSet.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>


using namespace std;

const int MAX_STR = 1024;


void load_reads(string file, map<string, int>& data, bool rev) {

	if (g_file_type == "fa")
		load_reads_fa(file, data, rev);
	else
		load_reads_fq(file, data, rev);


}



void load_reads_fa(string file, map<string, int>& data, bool rev) {

	time_t s_time = time(NULL);
	
	fstream in;
	in.open(file.c_str(), fstream::in);

	if (!in.is_open()) {

		cout << "Error! Can't open file " << file << endl;
		exit(1);

	}

	cout << "Loading reads from file " << file << " ..." << endl;

	char temp[MAX_STR];
	string read;
	in.getline(temp, MAX_STR);
	in.getline(temp, MAX_STR);
	read = temp;
	while (!in.eof()) {

		in.getline(temp, MAX_STR);	

		if(temp[0] == '>') {
		
            if (rev)
               	read = revcomp(read);
			if (data.find(read) != data.end())
				data[read] = data[read] + 1;
			else
				data[read] = 1;
			//data.push_back(read);			
			in.getline(temp, MAX_STR);
			read = temp;			
		} else {
			read = read + temp;		
		}

	}

	if (rev){
		read = revcomp(read);
	}
	//input_data.push_back(read);
	if (data.find(read) != data.end())
		data[read] = data[read] + 1;
	else
		data[read] = 1;

	in.close();

	time_t e_time = time(NULL);
	cout << "Success! (total cost time: " << (e_time - s_time) << " s)" << endl;

}


void load_reads_fq(string file, map<string, int>& data, bool rev) {

	time_t s_time = time(NULL);
	
	fstream in;
	in.open(file.c_str(), fstream::in);

	if (!in.is_open()) {

		cout << "Error! Can't open file " << file << endl;
		exit(1);

	}

	cout << "Loading reads from file " << file << " ..." << endl;

	char temp[MAX_STR];
	string read;
	in.getline(temp, MAX_STR);
	while (!in.eof()) {
		
		if(temp[0] == '@') {
			in.getline(temp, MAX_STR);
			read = temp;
			in.getline(temp, MAX_STR);
			while (temp[0] != '+') {
				read = read + temp;	
				in.getline(temp, MAX_STR);
			}
            if (rev)
               	read = revcomp(read);
			if (data.find(read) != data.end())
				data[read] = data[read] + 1;
			else
				data[read] = 1;
			//input_data.push_back(read);					
		} 
		in.getline(temp, MAX_STR);

	}

	in.close();

	time_t e_time = time(NULL);
	cout << "Success! (total cost time: " << (e_time - s_time) << " s)" << endl;

}


