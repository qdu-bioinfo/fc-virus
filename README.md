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
=======
git clone https://github.com/qdu-bioinfo/FC-Virus.git
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
g++ -o FC-Virus GeneralSet.cpp ReadUtility.cpp KmerUtility.cpp KmerHash.cpp HomoKmer.cpp Consensus.cpp FC-Virus.cpp
```
Then run FC-Virus:
```
./FC-Virus -t fq --left forward.fastq --right  reverse.fastq -o ./outfile/
```

Example command
```
./FC-Virus -t fq --left ./path/to/forward.fastq --right ./path/to/reverse.fastq -o ./path/to/outfile/
```
