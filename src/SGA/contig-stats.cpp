//
//  contig-stats.cpp
//  sga_git
//
//  Created by Milan Malinsky on 02/09/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//

#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include "Util.h"
#include "correct_trio.h"
#include "correct.h"
#include "SuffixArray.h"
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
#include "contig-stats.h"

// Functions

//
// Getopt
//
#define SUBPROGRAM "contig-stats"
static const char *STATS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Milan Malinsky.\n"
"\n"
"Copyright 2014 Wellcome Trust Sanger Institute\n";

static const char *STATS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... CONTIGSFILE\n"
"Print statistics about contigs (histogram of contig lenghts, N30, N50, N70, N90)\n"
"\n"
"      --help                           display this help and exit\n"
"      -b, --bin-size=NUM               bin size for the histogram of contig lengths (dafault 100bp)\n"
"      -m, --min-size=NUM               minimum contig size for it to be counted (dafault 100bp)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static std::string prefix;
    static double binSize = 100.0;
    static int minSize = 100;
    static std::string readsFile;
}

static const char* shortopts = "p:b:m:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "prefix",             required_argument, NULL, 'p' },
    { "min-size",           required_argument, NULL, 'm' },
    { "bin-size",           required_argument, NULL, 'b' },
    { "help",               no_argument,       NULL, OPT_HELP },
    { NULL, 0, NULL, 0 }
};


bool greaterThan (int i,int j) { return (i>j); }


//
// Main
//
int contigStatsMain(int argc, char** argv)
{
    parseContigStatsOptions(argc, argv);
    std::ifstream* readsFile = new std::ifstream(opt::readsFile.c_str());
    std::string contigSizesFileName = opt::prefix + "_sizes.txt"; std::ofstream* contigSizesFile = new std::ofstream(contigSizesFileName.c_str());
    std::string NXXFileName = opt::prefix + "_NXX.txt"; std::ofstream* NXXFile = new std::ofstream(NXXFileName.c_str());
    std::string line;
    std::vector<std::string::size_type> contigLengths;
    std::vector<std::string::size_type> flooredContigLengths;
    std::map<int,int> lengthCounts;
    
    while (getline(*readsFile, line)) {
        if (line[0] == '>') {
            continue;
        } else {
            std::string::size_type cLength = line.length();
            if (cLength >= opt::minSize) { // Only use contigs that have at least minSize length
                int flooredCLength = opt::binSize * floor(cLength/(double)opt::binSize);
                contigLengths.push_back(line.length());
                flooredContigLengths.push_back(flooredCLength);
            }
        }
    }
    int maxLcontig = (int)*std::max_element(flooredContigLengths.begin(), flooredContigLengths.end());
    int start; if (opt::binSize > 100) { start = 0; } else { start = 0; }
    // initialise lengthCounts
    for (int i = start; i <= maxLcontig; i = i + opt::binSize) {
        lengthCounts[i] = 0;
    }
    // and now fill them in
    for (int i = 0; i != flooredContigLengths.size(); i++) {
        lengthCounts[(int)flooredContigLengths[i]]++;
    }
    // and print
    for (int i = start; i <= maxLcontig; i = i + opt::binSize) {
        *contigSizesFile << i << "\t" << i+opt::binSize << "\t" << lengthCounts[i] << std::endl;
    }
    
    // Calculate N30, N50, N70, N90
    std::sort(contigLengths.begin(), contigLengths.end(), greaterThan);
    int totalLength = std::accumulate(contigLengths.begin(), contigLengths.end(), 0);
    int N30 = totalLength * 0.3; int N50 = totalLength * 0.5; int N70 = totalLength * 0.7; int N90 = totalLength * 0.9;
    bool N30done = false; bool N50done = false; bool N70done = false; bool N90done = false;
    
    int runningTotal = 0;
    *NXXFile << "Statistics for: " << opt::prefix << std::endl;
    *NXXFile << "Total assembly length:\t" << totalLength << std::endl;
    for (std::vector<std::string::size_type>::iterator it=contigLengths.begin(); it != contigLengths.end(); ++it) {
        runningTotal += *it;
        if (runningTotal >= N30 && !N30done) {
            *NXXFile << "N30:\t" << *it << std::endl;
            N30done = true;
        }
        if (runningTotal >= N50 && !N50done) {
            *NXXFile << "N50:\t" << *it << std::endl;
            N50done = true;
        }
        if (runningTotal >= N70 && !N70done) {
            *NXXFile << "N70:\t" << *it << std::endl;
            N70done = true;
        }
        if (runningTotal >= N90 && !N90done) {
            *NXXFile << "N90:\t" << *it << std::endl;
            N90done = true;
        }
    }
    
    //print_vector(contigLengths, *contigSizesFile, '\n');
    
    return 0;
}

//
// Handle command line arguments
//
void parseContigStatsOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'p': arg >> opt::prefix; break;
            case 'b': arg >> opt::binSize; break;
            case 'm': arg >> opt::minSize; break;
            case '?': die = true; break;
            case OPT_HELP:
                std::cout << STATS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << STATS_VERSION_MESSAGE;
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
    
    if (opt::binSize <= 50)
    {
        std::cerr << SUBPROGRAM ": bin size cannot be <= 50bp\n";
        die = true;
    }
    

    
    
    
    if (die)
    {
        std::cout << "\n" << STATS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::readsFile = argv[optind++];
    
    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }
}
