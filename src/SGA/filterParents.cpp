//
//  filterParents.cpp
//  sga_git
//
//  Created by Milan Malinsky on 17/10/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
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
#include "FilterParentProcess.h"
#include "CorrectionThresholds.h"
#include "KmerDistribution.h"
#include "BWTIntervalCache.h"
#include "LRAlignment.h"
#include "filterParents.h"

// Functions
int learnKmerParameters(const BWT* pBWT);
//
// Getopt
//
#define SUBPROGRAM "correctTrio"
static const char *FILTER_PARENTS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Milan Malinsky.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *FILTER_PARENTS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... PARENT_READSFILE OFFSPRING_READSFILE\n"
"Correct sequencing errors in the reads from PARENT_READSFILE, and only keep reads consistent with OFFSPRING_READSFILE\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -n, --do-not-correct             do not error-correct the reads, just check for consistency with the offspring\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the parent index files (default: prefix of the reads file)\n"
"      --offspring-prefix=PREFIX        use PREFIX for the names of the offspring index files (default: prefix of the reads file)\n"
"      -o, --outfile=FILE               write the corrected reads to FILE (default: READSFILE.ec.fa)\n"
"      -t, --threads=NUM                use NUM threads for the computation (default: 1)\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"          --metrics=FILE               collect error correction metrics (error rate by position in read, etc) and write them to FILE\n"
"\nKmer correction parameters:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 31)\n"
"      -i, --kmer-rounds=N              Perform N rounds of k-mer correction, correcting up to N bases (default: 10)\n"
"      -x, --parent-kmer-threshold=N    Manually set the kmer threshold for the mother\n"
"                                       (default: learn the k-mer correction threshold from the data)\n"
"          --offspring-kmer-threshold=N    Manually set the kmer threshold for the father\n"
"                                       (default: learn the k-mer correction threshold from the data)\n"
//"          --learn                      Attempt to learn the k-mer correction threshold (experimental). Overrides -x parameter.\n"
"          --do-not-discard             Low quality reads that couldn't be corrected will be kept in the output file(s)\n"
"                                       Useful if read filtering is done later via 'sga filter' command using index of corrected reads"
"                                       (DOES THIS MAKE A DIFFERENCE?)\n"
"                                       Also, when dealing with paired-end, ensures that the output file(s) still contain both reads of the pair\n"
"          --paired                     The reads are paired-end with records interleaved within a single file\n"
"                                       The read pairing information is going to be used in phasing\n"
"                                       i.e. separating reads to maternal_READSFILE_ec.fa and paternal_READSFILE_ec.fa\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string parentPrefix;
    static std::string offspringPrefix;
    static std::string parentReadFile;
    static std::string offspringReadFile;
    static std::string outFile;
    static std::string discardFile;
    static std::string metricsFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
    
    static bool bPaired = false;
    static bool bCorrect = true;
    static std::string consistentFile;
    static std::string inconsistentFile;
    
    static int numKmerRounds = 10;
    static int kmerLength = 31;
    static int parentKmerThreshold = 0;
    //    static bool bLearnKmerParams = false;
    static int intervalCacheLength = 10;
    
}

static const char* shortopts = "p:d:t:o:k:i:x:vn";

enum { OPT_HELP = 1, OPT_VERSION, OPT_METRICS, OPT_DISCARD, OPT_PAIR, OPT_OFFSPRING_THRESH, OPT_OFFSPRING_PREFIX };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "do-not-correct",       no_argument,       NULL, 'n' },
    { "threads",       required_argument, NULL, 't' },
    { "outfile",       required_argument, NULL, 'o' },
    { "prefix",        required_argument, NULL, 'p' },
    { "offspring-prefix",        required_argument, NULL, OPT_OFFSPRING_PREFIX },
    { "sample-rate",   required_argument, NULL, 'd' },
    { "kmer-size",     required_argument, NULL, 'k' },
    { "kmer-rounds",   required_argument, NULL, 'i' },
    { "parent-kmer-threshold",   required_argument, NULL, 'x' },
    { "offspring-kmer-threshold", required_argument, NULL, 'c' },
    //   { "learn",         no_argument,       NULL, OPT_LEARN },
    { "paired",         no_argument,       NULL, OPT_PAIR },
    { "do-not-discard",       no_argument,       NULL, OPT_DISCARD },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { "metrics",       required_argument, NULL, OPT_METRICS },
    { NULL, 0, NULL, 0 }
};



BWTIndexSet loadIndicesFilter(const std::string& readFile) {
    BWT* pBWT = new BWT(stripFilename(readFile) + BWT_EXT, opt::sampleRate);
    BWT* pRBWT = NULL;
    SampledSuffixArray* pSSA = NULL;
    
    BWTIntervalCache* pIntervalCache = new BWTIntervalCache(opt::intervalCacheLength, pBWT);
    
    BWTIndexSet indexSet;
    indexSet.pBWT = pBWT;
    indexSet.pRBWT = pRBWT;
    indexSet.pSSA = pSSA;
    indexSet.pCache = pIntervalCache;
    
    return indexSet;
}

void deleteIndicesFilter(BWTIndexSet& indexSet) {
    delete indexSet.pBWT;
    if(indexSet.pCache != NULL)
        delete indexSet.pCache;
    if(indexSet.pRBWT != NULL)
        delete indexSet.pRBWT;
    if(indexSet.pSSA != NULL)
        delete indexSet.pSSA;
}

//
// Main
//
int filterParentsMain(int argc, char** argv)
{
    parseFilterParentsOptions(argc, argv);
    
    std::cout << "Correcting sequencing errors for " << opt::parentReadFile << "\n";
    std::cout << "and keeping only reads consistent with: " << opt::offspringReadFile << "\n";
    
    // Load indices (  not allowing custom prefixes // BWT* pBWT = new BWT(opt::offspringPrefix + BWT_EXT, opt::sampleRate);)
    BWTIndexSet offspringIndexSet = loadIndicesFilter(opt::offspringPrefix + ".fastq");
    BWTIndexSet parentIndexSet;
    if (opt::bCorrect) {
        parentIndexSet = loadIndicesFilter(opt::parentPrefix + ".fastq");
    }

    
    // Learn the parameters of the kmer corrector
    int parentThreshold = 0;
    if (opt::bCorrect) {
        parentThreshold = (opt::parentKmerThreshold > 0) ? opt::parentKmerThreshold : learnKmerParameters(parentIndexSet.pBWT);
    }
    
    // Open outfiles and start a timer
    std::ostream* pWriter = createWriter(opt::outFile);
    std::ostream* pInconsistentWriter = createWriter(opt::discardFile);
    
    Timer* pTimer = new Timer(PROGRAM_IDENT);
    offspringIndexSet.pBWT->printInfo();
    
    // Set the error correction parameters
    FilterParentParameters fpParams;
    fpParams.offspringIndices = offspringIndexSet;
    fpParams.parentIndices = parentIndexSet;
    fpParams.parentThreshold = parentThreshold;
    fpParams.bCorrect = opt::bCorrect;
    fpParams.offspringThreshold = 2;
    fpParams.numKmerRounds = opt::numKmerRounds;
    fpParams.kmerLength = opt::kmerLength;
    
    // Setup post-processor
    bool bCollectMetrics = !opt::metricsFile.empty();
    
    FilterParentPostProcess* postProcessor;

    postProcessor = new FilterParentPostProcess(pWriter, pInconsistentWriter, bCollectMetrics, opt::bPaired);
    
    if(opt::numThreads <= 1)
    {
        // Serial mode
        FilterParentProcess processor(fpParams);
        if (opt::bPaired) {
            SequenceProcessFramework::processSequencesSerial<SequenceWorkItemPair,
            FilterParentPairResult,
            FilterParentProcess,
            FilterParentPostProcess>(opt::parentReadFile, &processor, postProcessor);
        } else {
            SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
            FilterParentResult,
            FilterParentProcess,
            FilterParentPostProcess>(opt::parentReadFile, &processor, postProcessor);
        }
    }
    else
    {
        // Parallel mode
        std::vector<FilterParentProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            FilterParentProcess* pProcessor = new FilterParentProcess(fpParams);
            processorVector.push_back(pProcessor);
        }
        
        if (opt::bPaired) {
            SequenceProcessFramework::processSequencesParallelOpenMP<SequenceWorkItemPair,
            FilterParentPairResult,
            FilterParentProcess,
            FilterParentPostProcess>(opt::parentReadFile, processorVector, postProcessor);
        } else {
            SequenceProcessFramework::processSequencesParallelOpenMP<SequenceWorkItem,
            FilterParentResult,
            FilterParentProcess,
            FilterParentPostProcess>(opt::parentReadFile, processorVector, postProcessor);
        }
        
        for(int i = 0; i < opt::numThreads; ++i)
        {
            delete processorVector[i];
        }
    }
    
    if(bCollectMetrics)
    {
        std::ostream* pMetricsWriter = createWriter(opt::metricsFile);
        postProcessor->writeMetrics(pMetricsWriter);
        delete pMetricsWriter;
    }
    
    deleteIndicesFilter(offspringIndexSet);
    deleteIndicesFilter(parentIndexSet);
    
    delete pTimer;
    
    delete pWriter;
    if(pInconsistentWriter != NULL)
        delete pInconsistentWriter;
    
    return 0;
}



//
// Handle command line arguments
//
void parseFilterParentsOptions(int argc, char** argv)
{
    std::string algo_str;
    bool bDiscardReads = true;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'p': arg >> opt::parentPrefix; break;
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'd': arg >> opt::sampleRate; break;
            case 'k': arg >> opt::kmerLength; break;
            case 'x': arg >> opt::parentKmerThreshold; break;
            case 'n': opt::bCorrect = false; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'i': arg >> opt::numKmerRounds; break;
                //  case OPT_LEARN: opt::bLearnKmerParams = true; break;
            case OPT_OFFSPRING_PREFIX: arg >> opt::offspringPrefix; break;
            case OPT_DISCARD: bDiscardReads = false; break;
            case OPT_PAIR: opt::bPaired = true; break;
            case OPT_METRICS: arg >> opt::metricsFile; break;
            case OPT_HELP:
                std::cout << FILTER_PARENTS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << FILTER_PARENTS_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 2)
    {
        std::cerr << SUBPROGRAM ": missing arguments - need reads from the parent and the offspring (in this order)\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }
    
    // Validate parameters
    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }
    
    if(opt::numKmerRounds <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of kmer rounds: " << opt::numKmerRounds << ", must be at least 1\n";
        die = true;
    }
    
    if(opt::kmerLength <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer length: " << opt::kmerLength << ", must be greater than zero\n";
        die = true;
    }
    
    if(opt::parentKmerThreshold < 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer threshold: " << opt::parentKmerThreshold << ", cannot be negative\n";
        die = true;
    }
    
    if (die)
    {
        std::cout << "\n" << FILTER_PARENTS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    
    // Parse the input filenames
    opt::parentReadFile = argv[optind++];
    opt::offspringReadFile = argv[optind++];
    
    if(opt::parentPrefix.empty())
    {
        opt::parentPrefix = stripFilename(opt::parentReadFile);
    }
    
    if(opt::offspringPrefix.empty())
    {
        opt::offspringPrefix = stripFilename(opt::offspringReadFile);
    }
    
    // Set the correction threshold
    //if(opt::kmerThreshold <= 0)
    //{
    //    std::cerr << "Invalid kmer support threshold: " << opt::kmerThreshold << "\n";
    //    exit(EXIT_FAILURE);
    //}
    // CorrectionThresholds::Instance().setBaseMinSupport(opt::kmerThreshold);
    
    std::string out_prefix = stripFilename(opt::offspringReadFile);
    if(!opt::outFile.empty()) {
        out_prefix = stripFilename(opt::outFile);
    } else {
        std::cout << "Need to provide a name for the outfile\n" << FILTER_PARENTS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    opt::discardFile = out_prefix + ".discarded.fa";
    
    
    if(bDiscardReads)
    {
        opt::discardFile = out_prefix + ".discard.fa";
    }
    else
    {
        opt::discardFile.clear();
    }
}