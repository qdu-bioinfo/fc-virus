# FC-Virus

FC-Virus is a powerful tool designed for full-length virus genome assembly. The package includes precompiled binaries for direct execution.

## Requirements

- **Operating System:** Linux
- **C++ Standard:** C++11 or higher

## Installation Instructions

You can install FC-Virus in the following ways:

### 1. Install via Bioconda
```
We have already submitted the application to Bioconda. Please stay tuned. Once it’s approved, you’ll be able to install it using the following command:
```
```
bioconda install FC-Virus
```
### 2. Use Precompiled Binaries:
Download FC-Virus from GitHub or clone the repository:
```
https://github.com/qdu-bioinfo/FC-Virus.git
```
 Ensure unzip is installed and unzip the downloaded file
 ```
unzip FC-Virus-master.zip
```
 Navigate to the FC-Virus directory
 ```
cd /path/to/FC-Virus
```
 Install by running make
 ```
make
```
 Set up environment variables
 ```
vim ~/.bashrc
```
 Add the following line to ~/.bashrc
 ```
export PATH="$PATH:/path/to/FC-Virus/bin"
```
 Apply the changes
 ```
source ~/.bashrc
```
### 3. Compile from source
To compile from source, navigate to the code directory:
```
cd /path/to/FC-Virus/code
```
Compile using g++:
```
g++ -o FC-Virus -std=c++11 main.cpp GeneralSet.cpp kmer.cpp ex_r.cpp ex_l.cpp
```
Then run FC-Virus:
```
./FC-Virus -o /path/to/output -k 25 --left /path/to/left.fq --right /path/to/right.fq
```

Example command
```
./FC-Virus -t fq --left ../path/to/forward.fastq --right ../path/to/reverse.fastq -o ../path/to/outfile/
```
