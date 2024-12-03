# FC-Virus

FC-Virus is a genome assembly algorithm designed to accurately reconstruct full-length consensus sequences for viral quasispecies. A key advantage of FC-Virus is its ability to produce a single consensus sequence that delivers the same assembly outcomes as multiple contigs from other assemblers.

## Requirements

- **Operating System:** Linux
- **C++ Standard:** C++11 or higher
- **Boost**
  ```
   wget http://downloads.sourceforge.net/project/boost/boost/1.80.0/boost_1_80_0.tar.gz
   tar xfz boost_1_80_0.tar.gz
   rm boost_1_80_0.tar.gz
   cd boost_1_80_0
   ./bootstrap.sh --prefix=/usr/local --with-libraries=program_options,regex,filesystem,system
   export
   ./b2 install
   cd /home
   rm -rf boost_1_80_0
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
  ```

## Installation Instructions

You can install fc-virus in the following ways:

### 1. Install via Bioconda
```
conda install fc
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
