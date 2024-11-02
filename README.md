# FC-Virus

FC-Virus is a powerful tool designed for full-length virus genome assembly. The package includes precompiled binaries for direct execution.

## Requirements

- **Operating System:** Linux
- **C++ Standard:** C++11 or higher

## Installation Instructions

You can install fc-virus in the following ways:

### 1. Install via Bioconda
```
We have already submitted the application to Bioconda. Please stay tuned. Once it’s approved, you’ll be able to install it using the following command:
```
```
bioconda install fc-virus
```
### 2. Use Precompiled Binaries:
Download fc-virus from GitHub or clone the repository:
```
git clone https://github.com/qdu-bioinfo/fc-virus.git
```
 Navigate to the fc-virus directory
 ```
cd /path/to/fc-virus
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
export PATH="$PATH:/path/to/fc-virus/bin"
```
 Apply the changes
 ```
source ~/.bashrc
```
### 3. Compile from source
To compile from source, navigate to the code directory:
```
cd /path/to/fc-virus/code
```
Compile using g++:
```
g++ -o fc-virus GeneralSet.cpp ReadUtility.cpp KmerUtility.cpp KmerHash.cpp HomoKmer.cpp Consensus.cpp FC-Virus.cpp
```
Then run fc-virus:
```
./fc-virus -t fq --left forward.fastq --right  reverse.fastq -o ./outfile/
```

Example command
```
./fc-virus -t fq --left ./path/to/forward.fastq --right ./path/to/reverse.fastq -o ./path/to/outfile/
```
