//
//  PCR-pair-removal.cpp
//  sga_git
//
//  Created by Milan Malinsky on 04/02/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "Util.h"
#include "SampledSuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "Timer.h"
#include "PCR-pair-removal.h"
#include "RemovePCRDuplicates.h"
#include "SequenceProcessFramework.h"
#include "Timer.h"
#include "SeqReader.h"
#include "Alphabet.h"
#include "BitVector.h"
#include "BWTDiskConstruction.h"


//
// Getopt
//
#define SUBPROGRAM "PCR-pair-removal"
static const char *SR_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Milan Malinsky.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *SR_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] READS_FILE\n"
"Sample reads from READS_FILE with a given probability\n"
"Useful for downsapling to a lower coverage and for randomly splitting reads into two files\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"\nInput/Output options:\n"
"      -o, --out=FILE                   write sampled reads to FILE\n"
"      -t, --threads=NUM                use NUM threads for the computation (default: 1)\n"
"      --pe-mode=INT                    0 - do not treat reads as paired\n"
"                                       1 - reads are paired and the records are interleaved within a single file (default)\n"
"      -d, --discard=DISCARD_FILE       Removed read pairs are output separately into DISCARD_FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static std::string outFile;
    static std::string discardFile;
    static unsigned int peMode = 1;
    static std::string readFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
}

static const char* shortopts = "o:d:t:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PEMODE };

static const struct option longopts[] = {
    { "verbose",                no_argument,       NULL, 'v' },
    { "out",                    required_argument, NULL, 'o' },
    { "threads",                required_argument, NULL, 't' },
    { "pe-mode",                required_argument, NULL, OPT_PEMODE },
    { "discard",                required_argument, NULL, 'd' },
    { "help",                   no_argument,       NULL, OPT_HELP },
    { "version",                no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int PCRpairRemovalMain(int argc, char** argv)
{
    parsePCRpairRemovalOptions(argc, argv);
    
    std::cerr << "Parameters:\n";
    
    std::cerr << "PE Mode: " << opt::peMode << "\n";
    std::cerr << "Unique pairs=>Outfile: " << opt::outFile << "\n";
    if(!opt::discardFile.empty()) {
        std::cerr << "PCR duplicates=>Outfile: " << opt::discardFile + "\n";
    }
       
    // Load indices (  not allowing custom prefixes // BWT* pBWT = new BWT(opt::offspringPrefix + BWT_EXT, opt::sampleRate);)
    BWT* pBWT = new BWT(stripFilename(opt::readFile) + BWT_EXT, opt::sampleRate);
    SampledSuffixArray* pSAI = new SampledSuffixArray(stripFilename(opt::readFile) + SAI_EXT, SSA_FT_SAI);
    
    
    // Open outfiles and start a timer
    std::ostream* pWriter = createWriter(opt::outFile);
    std::ostream* pDiscardWriter = (!opt::discardFile.empty() ? createWriter(opt::discardFile) : NULL);
    Timer* pTimer = new Timer(PROGRAM_IDENT);
    pBWT->printInfo();
    
    
    BitVector* pSharedBV = new BitVector(pBWT->getNumStrings());
    // Set the PCR removal parameters
    PCRRemoveParameters remParams;
    remParams.pBWT = pBWT;
    remParams.pSAI = pSAI;
    remParams.pSharedBV = pSharedBV;
    
    
    PCRRemovePostProcess* postProcessor = new PCRRemovePostProcess(pWriter, pDiscardWriter, opt::peMode);
    
    
    if(opt::numThreads <= 1)
    {
        // Serial mode
        PCRRemoveProcess processor(remParams); 
        if (opt::peMode) {
            SequenceProcessFramework::processSequencesSerial<SequenceWorkItemPair,
            PCRRemoveResult, 
            PCRRemoveProcess, 
            PCRRemovePostProcess>(opt::readFile, &processor, postProcessor);
        } 
    }
    else
    {
        // Parallel mode
        std::vector<PCRRemoveProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            PCRRemoveProcess* pProcessor = new PCRRemoveProcess(remParams);
            processorVector.push_back(pProcessor);
        }
        
        if (opt::peMode) {
            SequenceProcessFramework::processSequencesParallelOpenMP<SequenceWorkItemPair,
            PCRRemoveResult, 
            PCRRemoveProcess, 
            PCRRemovePostProcess>(opt::readFile, processorVector, postProcessor);
        }    
        
        for(int i = 0; i < opt::numThreads; ++i)
        {
            delete processorVector[i];
        }
    }
   
    /*
    if(bCollectMetrics)
    {
        std::ostream* pMetricsWriter = createWriter(opt::metricsFile);
        postProcessor->writeMetrics(pMetricsWriter);
        delete pMetricsWriter;
    }
    */
    
    // Cleanup
    delete postProcessor;
    if(pSharedBV != NULL)
        delete pSharedBV;
    delete pBWT;
    delete pSAI;
    delete pWriter;
    if(pDiscardWriter != NULL)
        delete pDiscardWriter;
    
    std::cout << "Updating the BWT and Suffix Array index files\n";
    // Rebuild the FM-index without the discarded reads
    std::string out_prefix = stripFilename(opt::outFile);
    removeReadsFromIndices(stripFilename(opt::readFile), opt::discardFile, out_prefix, BWT_EXT, SAI_EXT, false, opt::numThreads);
    delete pTimer;

    return 0;
}


//
// Handle command line arguments
//
void parsePCRpairRemovalOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'o': arg >> opt::outFile; break;
            case 'd': arg >> opt::discardFile; break;
            case 't': arg >> opt::numThreads; break;                
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_PEMODE: arg >> opt::peMode; break;
            case OPT_HELP:
                std::cout << SR_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << SR_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 1)
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }
    
    if (die)
    {
        std::cout << "\n" << SR_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    opt::readFile = argv[optind++];
    
    // Validate parameters
    if(opt::numThreads <= 0) {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }
    if (opt::peMode > 1) {
        std::cerr << SUBPROGRAM ": invalid pe-mode: " << opt::numThreads << "\n";
        die = true;
    }
    
    if(opt::outFile.empty()) {
        std::string out_prefix = stripExtension(opt::readFile);
        opt::outFile = out_prefix + ".PCRrem.fa";
    } else {
        std::string out_prefix = stripExtension(opt::outFile);
        opt::outFile = out_prefix + ".PCRrem.fa";
    }
    
    if(opt::discardFile.empty()) {
        std::string prefix = stripExtension(opt::readFile);
        opt::discardFile = prefix + ".PCRdiscard.fa";
    } else {
        std::string prefix = stripExtension(opt::discardFile);
        opt::discardFile = prefix + ".PCRdiscard.fa";
    }
    
    if (die)
    {
        std::cout << "\n" << SR_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    if(opt::peMode > 1)
    {
        std::cerr << SUBPROGRAM ": error pe-mode must be 0 or 1 (found: " << opt::peMode << ")\n";
        exit(EXIT_FAILURE);
    }
    
}
