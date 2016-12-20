# SIFT4G

SIFT (Sorting Intolerant From Tolerant) For Genomes (http://www.nature.com/nprot/journal/v11/n1/full/nprot.2015.123.html)

## Requirements
- g++ (4.8+)
- GNU Make
- nvcc (2.\*+) (optional)

\*note: It was only tested on Linux (Ubuntu 14.04).

## Installation

To build SIFT4G run the following commands from your terminal:

    git clone --recursive https://github.com/rvaser/sift4g.git sift4g
    cd sift4g/
    make

Running the 'make' command will create the bin folder which contains the sift4g executable.

**Troubleshooting** If you accidentally left out '--recursive' from git clone, run the following commands before running 'make':

    git submodule init
    git submodule update

**Optional** If you have a CUDA enabled graphics card (and the nvcc compiler) run 'make gpu' instead of 'make'.

## Usage

To run the default version of SIFT4G run the following command:

    ./bin/sift4g -q <query file> -d <database file>

To see all available parameters run the command bellow:

    ./bin/sift4g -h
