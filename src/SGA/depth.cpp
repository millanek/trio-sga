//
//  depth.cpp
//  sga_git
//
//  Created by Milan Malinsky on 03/09/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//


#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include "depth.h"
#include "Util.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "TrioCorrectProcess.h"
#include "CorrectionThresholds.h"
#include "KmerDistribution.h"
#include "BWTIntervalCache.h"
#include "LRAlignment.h"

// Functions

//
// Getopt
//
#define SUBPROGRAM "depth"
static const char *DEPTH_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Milan Malinsky.\n"
"\n"
"Copyright 2014 Wellcome Trust Sanger Institute\n";

static const char *DEPTH_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... CONTIGSFILE\n"
"Print statistics about contigs (histogram of contig lenghts, N30, N50, N70, N90)\n"
"\n"
"      --help                           display this help and exit\n"
"      -s, --genome-size=NUM            (required) An estimate of genome size in bp\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static std::string prefix;
    static std::string readsFile;
    static int genomeSize;
}

static const char* shortopts = "p:s:";

enum { OPT_HELP = 1, OPT_GENOME_SIZE };

static const struct option longopts[] = {
    { "prefix",             required_argument, NULL, 'p' },
    { "genome-size",        required_argument, NULL, 's' },
    { "help",               no_argument,       NULL, OPT_HELP },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int depthMain(int argc, char** argv)
{
    parseDepthOptions(argc, argv);
    std::ifstream* readsFile = new std::ifstream(opt::readsFile.c_str());
    //std::string contigSizesFileName = opt::prefix + "_sizes.txt"; std::ofstream* contigSizesFile = new std::ofstream(contigSizesFileName.c_str());
    std::string line;
    int totalLength = 0;
 
    
    getline(*readsFile, line);
    // Try to figure out if this is fasta or a fastq file
    bool fastq = false; bool fasta = false;
    if (line[0] == '>') { fasta = true; }
    else if (line[0] == '@') { fastq = true; }
    else {
        std::cerr << "Input file type not recognised. Exitting..." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Now calculate read coverage
    if (fasta) {
        while (getline(*readsFile, line)) {
            if (line[0] == '>') {
                continue;
            } else {
                std::string::size_type cLength = line.length();
                totalLength += cLength;
            }
        }
    }
    int row = 4;
    if (fastq) {
        while (getline(*readsFile, line)) {
            if (row % 4 == 0) {
                std::string::size_type cLength = line.length();
                totalLength += cLength;
            }
            row++;
        }
    }
    
    std::cout << "Total amount of sequence:\t" << totalLength << std::endl;
    std::cout << "Approximate coverage:\t" << (double)totalLength/opt::genomeSize << std::endl;
    
    return 0;
}

//
// Handle command line arguments
//
void parseDepthOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'p': arg >> opt::prefix; break;
            case 's': arg >> opt::genomeSize; break;
            case '?': die = true; break;
            case OPT_HELP:
                std::cout << DEPTH_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 1)
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }
    else if (argc - optind > 1)
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }
    
    
    if (die)
    {
        std::cout << "\n" << DEPTH_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::readsFile = argv[optind++];
    
    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }
}
