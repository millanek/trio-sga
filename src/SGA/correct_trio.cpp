//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// correct - Correct sequencing errors in reads using the FM-index
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
#include "TrioCorrectProcess.h"
#include "CorrectionThresholds.h"
#include "KmerDistribution.h"
#include "BWTIntervalCache.h"
#include "LRAlignment.h"

// Functions
int learnKmerParameters(const BWT* pBWT);

//
// Getopt
//
#define SUBPROGRAM "correctTrio"
static const char *CORRECT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Milan Malinsky.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *CORRECT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... OFFSPRING_READSFILE MOTHER_READSFILE FATHER_READSFILE\n"
"Correct sequencing errors in all the reads in OFFSPRING_READSFILE, using information from trio sequencing\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -p, --prefix=PREFIX              use PREFIX for the names of the index files (default: prefix of the input file)\n"
"      -m, --mother-prefix=PREFIX              use PREFIX for the names of the mother index files (default: prefix of the mother input file)\n"
"      -f, --father-prefix=PREFIX              use PREFIX for the names of the father index files (default: prefix of the father input file)\n"
"      -o, --outfile=FILE               write the corrected reads to FILE (default: READSFILE.ec.fa)\n"
"      -t, --threads=NUM                use NUM threads for the computation (default: 1)\n"
"      -d, --sample-rate=N              use occurrence array sample rate of N in the FM-index. Higher values use significantly\n"
"                                       less memory at the cost of higher runtime. This value must be a power of 2 (default: 128)\n"
"          --metrics=FILE               collect error correction metrics (error rate by position in read, etc) and write them to FILE\n"
"\nKmer correction parameters:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 31)\n"
"      -i, --kmer-rounds=N              Perform N rounds of k-mer correction, correcting up to N bases (default: 10)\n"
"      -x, --off-kmer-threshold=N       Manually set the kmer correction threshold for offspring\n"
"                                       (default: learn the k-mer correction threshold from the data)\n"
"      -c, --off-cond-threshold=N       Manually set the correction threshold for offspring kmers not found in parents\n"
"                                       (default: off-kmer-threshold+2)\n"
"          --mother-kmer-threshold=N    Manually set the kmer threshold for the mother\n"
"                                       (default: learn the k-mer correction threshold from the data)\n"
"          --father-kmer-threshold=N    Manually set the kmer threshold for the father\n"
"                                       (default: learn the k-mer correction threshold from the data)\n"
//"          --learn                      Attempt to learn the k-mer correction threshold (experimental). Overrides -x parameter.\n"
"          --do-not-discard             Low quality reads that couldn't be corrected will be kept in the output file(s)\n"
"                                       Useful if read filtering is done later via 'sga filter' command using index of corrected reads" 
"                                       (DOES THIS MAKE A DIFFERENCE?)\n"        
"                                       Also, when dealing with paired-end, ensures that the output file(s) still contain both reads of the pair\n"
"          --paired                     The reads are paired-end with records interleaved within a single file\n"
"                                       The read pairing information is going to be used in phasing\n"
"                                       i.e. separating reads to maternal_READSFILE_ec.fa and paternal_READSFILE_ec.fa\n"
"          --phase                      Separate offspring reads into two files for assembling maternal and paternal chromosomes separately\n"
"                                       Filenames will be: maternal_READSFILE_ec.fa paternal_READSFILE_ec.fa\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string offspringPrefix;
    static std::string motherPrefix;
    static std::string fatherPrefix;
    static std::string offspringReadFile;
    static std::string motherReadFile;
    static std::string fatherReadFile;
    static std::string outFile;
    static std::string discardFile;
    static std::string metricsFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
    
    static bool bPhase = false;
    static bool bPaired = false;
    static std::string maternalFile;
    static std::string paternalFile;
    static std::string neitherParentFile;
    static std::string inconsistentPhaseFile;
    
    static int kmerLength = 31;
    static int kmerThreshold = 0;
    static int conditionalKmerThreshold = 0;
    static int motherThreshold = 0;
    static int fatherThreshold = 0;
    static int numKmerRounds = 10;
//    static bool bLearnKmerParams = false;
    static int intervalCacheLength = 10;
    
}

static const char* shortopts = "p:d:t:o:k:i:x:c:m:f:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_METRICS, OPT_DISCARD, OPT_PHASE, OPT_PAIR, OPT_MOTHER_THRESH, OPT_FATHER_THRESH  };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't' },
    { "outfile",       required_argument, NULL, 'o' },
    { "prefix",        required_argument, NULL, 'p' },
    { "mother-prefix",        required_argument, NULL, 'm' },
    { "father-prefix",        required_argument, NULL, 'f' },
    { "sample-rate",   required_argument, NULL, 'd' },
    { "kmer-size",     required_argument, NULL, 'k' },
    { "kmer-rounds",   required_argument, NULL, 'i' },
    { "off-kmer-threshold",   required_argument, NULL, 'x' },
    { "off-cond-threshold", required_argument, NULL, 'c' },
    { "mother-kmer-threshold",   required_argument, NULL, OPT_MOTHER_THRESH },
    { "father-kmer-threshold",   required_argument, NULL, OPT_FATHER_THRESH },
 //   { "learn",         no_argument,       NULL, OPT_LEARN },
    { "paired",         no_argument,       NULL, OPT_PAIR },
    { "phase",         no_argument,       NULL, OPT_PHASE },
    { "do-not-discard",       no_argument,       NULL, OPT_DISCARD },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { "metrics",       required_argument, NULL, OPT_METRICS },
    { NULL, 0, NULL, 0 }
};



BWTIndexSet loadIndicesCT(const std::string& readFile) {
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

void deleteIndicesCT(BWTIndexSet& indexSet) {
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
int correctTrioMain(int argc, char** argv)
{
    parseCorrectTrioOptions(argc, argv);
    
    std::cout << "Correcting sequencing errors for " << opt::offspringReadFile << "\n";
    std::cout << "Taking into account reads from mother: " << opt::motherReadFile << " and father: "<< opt::fatherReadFile << "\n";
    
    // Load indices (  not allowing custom prefixes // BWT* pBWT = new BWT(opt::offspringPrefix + BWT_EXT, opt::sampleRate);)
    BWTIndexSet offspringIndexSet = loadIndicesCT(opt::offspringPrefix + ".fastq");
    BWTIndexSet motherIndexSet = loadIndicesCT(opt::motherPrefix + ".fastq");
    BWTIndexSet fatherIndexSet = loadIndicesCT(opt::fatherPrefix + ".fastq");
     
    
    
    // Learn the parameters of the kmer corrector
    int offspringThreshold = (opt::kmerThreshold > 0) ? opt::kmerThreshold : learnKmerParameters(offspringIndexSet.pBWT);
    int conditionalThreshold = (opt::conditionalKmerThreshold > 0) ? opt::conditionalKmerThreshold : offspringThreshold + 2; 
    int motherThreshold = (opt::motherThreshold > 0) ? opt::kmerThreshold : learnKmerParameters(motherIndexSet.pBWT);
    int fatherThreshold = (opt::fatherThreshold > 0) ? opt::kmerThreshold : learnKmerParameters(fatherIndexSet.pBWT);
    
    // Open outfiles and start a timer
    std::ostream* pWriter = (!opt::bPhase ? createWriter(opt::outFile) : NULL);
    std::ostream* pMatWriter = (opt::bPhase ? createWriter(opt::maternalFile) : NULL);
    std::ostream* pPatWriter = (opt::bPhase ? createWriter(opt::paternalFile) : NULL);
    std::ostream* pDiscardWriter = (!opt::discardFile.empty() ? createWriter(opt::discardFile) : NULL);
    std::ostream* pNeitherParentWriter = (opt::bPhase ? createWriter(opt::neitherParentFile) : NULL);
    std::ostream* pInconsistentPhaseWriter = (opt::bPhase ? createWriter(opt::inconsistentPhaseFile) : NULL);
    
    Timer* pTimer = new Timer(PROGRAM_IDENT);
    offspringIndexSet.pBWT->printInfo();
    
    // Set the error correction parameters
    TrioCorrectParameters ecParams;
    ecParams.offspringIndices = offspringIndexSet;
    ecParams.motherIndices = motherIndexSet;
    ecParams.fatherIndices = fatherIndexSet;
    ecParams.offspringThreshold = offspringThreshold;
    ecParams.motherThreshold = motherThreshold;
    ecParams.fatherThreshold = fatherThreshold;
    ecParams.numKmerRounds = opt::numKmerRounds;
    ecParams.kmerLength = opt::kmerLength;
    ecParams.offspringConditionalThreshold = conditionalThreshold;
    
    // Setup post-processor
    bool bCollectMetrics = !opt::metricsFile.empty();
    
    TrioCorrectPostProcess* postProcessor;
    if (opt::bPhase) {
        postProcessor = new TrioCorrectPostProcess(pMatWriter, pPatWriter, pDiscardWriter, pNeitherParentWriter,pInconsistentPhaseWriter,bCollectMetrics, opt::bPaired); 
    } else {
        postProcessor = new TrioCorrectPostProcess(pWriter, pDiscardWriter, bCollectMetrics, opt::bPaired);
    }
    if(opt::numThreads <= 1)
    {
        // Serial mode
        TrioCorrectProcess processor(ecParams); 
        if (opt::bPaired) {
            SequenceProcessFramework::processSequencesSerial<SequenceWorkItemPair,
            TrioCorrectPairResult, 
            TrioCorrectProcess, 
            TrioCorrectPostProcess>(opt::offspringReadFile, &processor, postProcessor);
        } else {
            SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
            TrioCorrectResult, 
            TrioCorrectProcess, 
            TrioCorrectPostProcess>(opt::offspringReadFile, &processor, postProcessor);
        }
    }
    else
    {
        // Parallel mode
        std::vector<TrioCorrectProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            TrioCorrectProcess* pProcessor = new TrioCorrectProcess(ecParams);
            processorVector.push_back(pProcessor);
        }
        
        if (opt::bPaired) {
            SequenceProcessFramework::processSequencesParallelOpenMP<SequenceWorkItemPair,
            TrioCorrectPairResult, 
            TrioCorrectProcess, 
            TrioCorrectPostProcess>(opt::offspringReadFile, processorVector, postProcessor);
        } else {    
            SequenceProcessFramework::processSequencesParallelOpenMP<SequenceWorkItem,
            TrioCorrectResult, 
            TrioCorrectProcess, 
            TrioCorrectPostProcess>(opt::offspringReadFile, processorVector, postProcessor);
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
    
    deleteIndicesCT(offspringIndexSet);
    deleteIndicesCT(motherIndexSet);
    deleteIndicesCT(fatherIndexSet);
    
    delete pTimer;
    
    delete pWriter;
    if(pDiscardWriter != NULL)
        delete pDiscardWriter;
    
    return 0;
}


// 
// Handle command line arguments
//
void parseCorrectTrioOptions(int argc, char** argv)
{
    std::string algo_str;
    bool bDiscardReads = true;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'p': arg >> opt::offspringPrefix; break;
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'd': arg >> opt::sampleRate; break;
            case 'k': arg >> opt::kmerLength; break;
            case 'x': arg >> opt::kmerThreshold; break;
            case 'c': arg >> opt::conditionalKmerThreshold; break;    
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'i': arg >> opt::numKmerRounds; break;
            case 'm': arg >> opt::motherPrefix; break;
            case 'f': arg >> opt::fatherPrefix; break;
          //  case OPT_LEARN: opt::bLearnKmerParams = true; break;
            case OPT_DISCARD: bDiscardReads = false; break;
            case OPT_PHASE: opt::bPhase = true; break;   
            case OPT_PAIR: opt::bPaired = true; break;
            case OPT_METRICS: arg >> opt::metricsFile; break;
            case OPT_MOTHER_THRESH: arg >> opt::motherThreshold; break;
            case OPT_FATHER_THRESH: arg >> opt::fatherThreshold; break;
            case OPT_HELP:
                std::cout << CORRECT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CORRECT_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 3) 
    {
        std::cerr << SUBPROGRAM ": missing arguments - need reads from offspring mother father (in this order)\n";
        die = true;
    } 
    else if (argc - optind > 3) 
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
    
    if(opt::kmerThreshold < 0)
    {
        std::cerr << SUBPROGRAM ": invalid kmer threshold: " << opt::kmerThreshold << ", cannot be negative\n";
        die = true;
    }
    
    if (die) 
    {
        std::cout << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    
    // Parse the input filenames
    opt::offspringReadFile = argv[optind++];
    opt::motherReadFile = argv[optind++];
    opt::fatherReadFile = argv[optind++];
    
    if(opt::offspringPrefix.empty())
    {
        opt::offspringPrefix = stripFilename(opt::offspringReadFile);
    }
    if(opt::motherPrefix.empty())
    {
        opt::motherPrefix = stripFilename(opt::motherReadFile);
    }
    
    if(opt::fatherPrefix.empty())
    {
        opt::fatherPrefix = stripFilename(opt::fatherReadFile);
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
    }
    if(opt::bPhase) {
        opt::maternalFile = "maternal_" + out_prefix + ".ec.fa";
        opt::paternalFile = "paternal_" + out_prefix + ".ec.fa";
        opt::neitherParentFile = "neither_" + out_prefix + ".fa";
        opt::inconsistentPhaseFile = "inconsistent_" + out_prefix + ".fa";
    } else {
        opt::outFile = out_prefix + ".ec.fa";
    }    
    if(bDiscardReads)
    {
        opt::discardFile = out_prefix + ".discard.fa";
    }
    else
    {
        opt::discardFile.clear();
    }
}