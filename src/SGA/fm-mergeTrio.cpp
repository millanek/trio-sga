//
//  fm-mergeTrio.cpp
//  sga_git
//
//  Created by Milan Malinsky on 19/04/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//

#include "fm-mergeTrio.h"
// fm-merge - Merge reads/sequences using the FM-index
// without explicitly constructing the full graph
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "fm-merge.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "OverlapProcess.h"
#include "ReadInfoTable.h"
#include "FMMergeTrioProcess.h"

//
// Getopt
//
#define SUBPROGRAM "fm-merge-trio"
static const char *FMMERGETRIO_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Milan Malinsky.\n"
"\n"
"Copyright 2014 Wellcome Trust Sanger Institute\n";

static const char *FMMERGETRIO_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE MOTHER_READSFILE FATHER_READSFILE\n"
"Merge unambiguously sequences from the READSFILE using the FM-index.\n"
"This program requires rmdup to be run before it.\n"
"\n"
"      --help                           display this help and exit\n"
"      -c, --contrib-parent={1,2}       !!(required)!! 1 if assembling the maternally contributed chromosome, 2 for paternal\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM worker threads (default: no threading)\n"
"      -m, --min-overlap=LEN            minimum overlap required between two reads to merge (default: 45)\n"
"      -o, --outfile=FILE               write the merged sequences to FILE (default: basename.merged.fa)\n"
"      -p, --try-phasing                (experimental) Try to use the parent's reads to do phasing during assembly\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int numThreads = 1;
    static int contribParent = 0;
    static std::string readsFile;
    static std::string motherReadFile;
    static std::string fatherReadFile;
    static std::string outFile;
    static std::string prefix;
    static bool bPhase = false;
    static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
}

static const char* shortopts = "p:m:d:e:t:l:s:o:vixc:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "contrib-parent",     required_argument,       NULL, 'c' },
    { "try-phasing",     no_argument,       NULL, 'p' },
    { "verbose",     no_argument,       NULL, 'v' },
    { "threads",     required_argument, NULL, 't' },
    { "min-overlap", required_argument, NULL, 'm' },
    { "outfile",     required_argument, NULL, 'o' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};



BWTIndexSet loadForwardReverseBWT(const std::string& readFile) {
    BWT* pBWT = new BWT(stripFilename(readFile) + BWT_EXT);
    BWT* pRBWT = new BWT(stripFilename(readFile) + RBWT_EXT);
    SampledSuffixArray* pSSA = NULL;
    
    BWTIndexSet indexSet;
    indexSet.pBWT = pBWT;
    indexSet.pRBWT = pRBWT;
    indexSet.pSSA = pSSA;
    
    return indexSet;
}

//
// Main
//
int FMMergeTrioMain(int argc, char** argv)
{
    parseFMMergeTrioOptions(argc, argv);
    
    BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
    BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT);
    OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT,0.0f, 0,0,true);
    pOverlapper->setExactModeOverlap(true);
    pOverlapper->setExactModeIrreducible(true);
    Timer* pTimer = new Timer(PROGRAM_IDENT);
    pBWT->printInfo();
    
    BWTIndexSet motherIndices = loadForwardReverseBWT(opt::motherReadFile);
    BWTIndexSet fatherIndices = loadForwardReverseBWT(opt::fatherReadFile);
    
    // Construct a bitvector indicating what reads have been used
    // All the processes read from this vector and only the post processor
    // writes to it.
    BitVector markedReads(pBWT->getNumStrings());
    
    std::ostream* pWriter = createWriter(opt::outFile);
    FMMergePostProcess postProcessor(pWriter, &markedReads);
    ContributingParent cp;
    if (opt::contribParent == 1) { cp = MATERNAL; }
    if (opt::contribParent == 2) { cp = PATERNAL; }
    
    if(opt::numThreads <= 1)
    {
        printf("[%s] starting serial-mode read merging\n", PROGRAM_IDENT);
        FMMergeTrioProcess processor(pOverlapper, opt::minOverlap, &markedReads,motherIndices,fatherIndices, cp, opt::verbose, opt::bPhase);
        SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
        FMMergeResult,
        FMMergeTrioProcess,
        FMMergePostProcess>(opt::readsFile, &processor, &postProcessor);
    }
    else
    {
        printf("[%s] starting parallel-mode read merging computation with %d threads\n", PROGRAM_IDENT, opt::numThreads);
        
        std::vector<FMMergeTrioProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            FMMergeTrioProcess* pProcessor = new FMMergeTrioProcess(pOverlapper, opt::minOverlap, &markedReads,motherIndices,fatherIndices, cp, opt::verbose, opt::bPhase);
            processorVector.push_back(pProcessor);
        }
        
        SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
        FMMergeResult,
        FMMergeTrioProcess,
        FMMergePostProcess>(opt::readsFile, processorVector, &postProcessor);
        
        for(size_t i = 0; i < processorVector.size(); ++i)
        {
            delete processorVector[i];
            processorVector[i] = NULL;
        }
    }
    
    // Check that every bit was set in the bit vector
    size_t numSet = 0;
    size_t numTotal = pBWT->getNumStrings();
    for(size_t i = 0; i < numTotal; ++i)
    {
        if(markedReads.test(i))
            ++numSet;
    }
    
    // Get the number of strings in the BWT, this is used to pre-allocated the read table
    delete pOverlapper;
    delete pBWT;
    delete pRBWT;
    delete pWriter;
    
    // Cleanup
    delete pTimer;
    if(opt::numThreads > 1)
        pthread_exit(NULL);
    
    return 0;
}

//
// Handle command line arguments
//
void parseFMMergeTrioOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'm': arg >> opt::minOverlap; break;
            case 'p': opt::bPhase = true; break;
            case 'o': arg >> opt::outFile; break;
            case 't': arg >> opt::numThreads; break;
            case 'c': arg >> opt::contribParent; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << FMMERGETRIO_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << FMMERGETRIO_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::contribParent <= 0) {
        std::cerr << SUBPROGRAM ": you need to specify the contributing parent for assembly (1-maternal, 2-paternal): " << opt::numThreads << "\n";
        die = true;
    }
    
    if (opt::contribParent >= 3) {
        std::cerr << SUBPROGRAM ": you need to specify the contributing parent for assembly (1-maternal, 2-paternal): " << opt::numThreads << "\n";
        die = true;
    }
    
    if (argc - optind < 3)
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }
    else if (argc - optind > 3)
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }
    
    if(opt::numThreads <= 0)
    {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
        die = true;
    }
    
    if (die)
    {
        std::cout << "\n" << FMMERGETRIO_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::readsFile = argv[optind++];
    opt::prefix = stripFilename(opt::readsFile);
    opt::motherReadFile = argv[optind++];
    opt::fatherReadFile = argv[optind++];
    
    if(opt::outFile.empty())
        opt::outFile = opt::prefix + ".merged.fa";
}
