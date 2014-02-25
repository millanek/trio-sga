//
//  sample-reads.cpp
//  sga_git
//
//  Created by Milan Malinsky on 19/01/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "Util.h"
#include "sample-reads.h"
#include "Timer.h"
#include "SeqReader.h"
#include "Alphabet.h"


//
// Getopt
//
#define SUBPROGRAM "sample-reads"
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
"      -o, --out=FILE                   write sampled reads to FILE (default: stdout)\n"
"      --pe-mode=INT                    0 - do not treat reads as paired\n"
"                                       1 - reads are paired and the records are interleaved within a single file (default)\n"
"      -p, --probability=FLOAT          Randomly sample reads or pairs with acceptance probability FLOAT.\n"
"      -d, --discard=DISCARD_FILE       Reads not sampled are output separately into DISCARD_FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string outFile;
    static std::string discardFile;
    static unsigned int peMode = 1;
    static double sampleFreq = 1.0f;
    static std::string readFile;
}

static const char* shortopts = "o:p:d:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PEMODE };

static const struct option longopts[] = {
    { "verbose",                no_argument,       NULL, 'v' },
    { "out",                    required_argument, NULL, 'o' },
    { "probability",            required_argument, NULL, 'p' },
    { "pe-mode",                required_argument, NULL, OPT_PEMODE },
    { "discard",            required_argument, NULL, 'd' },
    { "help",                   no_argument,       NULL, OPT_HELP },
    { "version",                no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

static int64_t s_numReadsRead = 0;
static int64_t s_numReadsKept = 0;
static int64_t s_numBasesRead = 0;
static int64_t s_numBasesKept = 0;
static int64_t s_numInvalidPE = 0;

//
// Main
//
int sampleReadsMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga sample-reads");
    parseSampleReadsOptions(argc, argv);
    
    std::cerr << "Parameters:\n";
    

    std::cerr << "Sample freq: " << opt::sampleFreq << "\n";
    std::cerr << "PE Mode: " << opt::peMode << "\n";
    std::cerr << "Sampled reads=>Outfile: " << (opt::outFile.empty() ? "stdout" : opt::outFile) << "\n";
    if(!opt::discardFile.empty()) {
        std::cerr << "Unsampled reads=>Outfile: " << opt::discardFile + "\n";
    }

    
    // Seed the RNG
    srand(time(NULL));
    
    std::ostream* pWriter;
    if(opt::outFile.empty())
    {
        pWriter = &std::cout;
    }
    else
    {
        std::ofstream* pFile = new std::ofstream(opt::outFile.c_str());
        assertFileOpen(*pFile, opt::outFile);
        pWriter = pFile;
    }
    
    // Create a filehandle to write unsampled reads to, if necessary
    std::ostream* pDiscardWriter = NULL;
    if(!opt::discardFile.empty())
        pDiscardWriter = createWriter(opt::discardFile);
    
    
    std::string filename = opt::readFile;
    std::cerr << "Processing " << filename << "\n\n";
    if(opt::peMode == 0)
    {
        SeqReader reader(filename, SRF_NO_VALIDATION);
        SeqRecord record;
        
        while(reader.get(record))
        {
            std::string seqStr = record.seq.toString();
            ++s_numReadsRead;
            s_numBasesRead += seqStr.size();

            if(samplePassTF())
            {
                record.write(*pWriter);
                ++s_numReadsKept;
                s_numBasesKept += record.seq.length();
            } else {
                if (!opt::discardFile.empty()) {
                    record.write(*pDiscardWriter);
                }
            }   
        }
    }
    else
    {
        assert(opt::peMode == 1);
        
        SeqReader* pReader1;
        SeqReader* pReader2;
        
        // Read from a single file
        pReader1 = new SeqReader(filename, SRF_NO_VALIDATION);
        pReader2 = pReader1;
        std::cerr << "Processing interleaved pe file " << filename << "\n";
        
        SeqRecord record1;
        SeqRecord record2;
        while(pReader1->get(record1) && pReader2->get(record2))
        {
            // If the names of the records are the same, append a /1 and /2 to them
            if(record1.id == record2.id)
            {
                record1.id.append("/1");
                record2.id.append("/2");
            }
            
            // Ensure the read names are sensible
            std::string expectedID2 = getPairID(record1.id);
            std::string expectedID1 = getPairID(record2.id);
            
            if(expectedID1 != record1.id || expectedID2 != record2.id)
            {
                std::cerr << "Warning: Pair IDs do not match (expected format /1,/2 or /A,/B)\n";
                std::cerr << "Read1 ID: " << record1.id << "\n";
                std::cerr << "Read2 ID: " << record2.id << "\n";
                s_numInvalidPE += 2;
            }
            
            std::string seqStr = record1.seq.toString();
            ++s_numReadsRead;
            s_numBasesRead += seqStr.size();            
            seqStr = record2.seq.toString();
            ++s_numReadsRead;
            s_numBasesRead += seqStr.size();
            
            if(samplePassTF()) {
                record1.write(*pWriter);
                record2.write(*pWriter);
                s_numReadsKept += 2;
                s_numBasesKept += record1.seq.length();
                s_numBasesKept += record2.seq.length();
            } else {
                if (!opt::discardFile.empty()) {
                    record1.write(*pDiscardWriter);
                    record2.write(*pDiscardWriter);
                }
            }    
        }
        
        if(pReader2 != pReader1)
        {
            // only delete reader2 if it is a distinct pointer
            delete pReader2;
            pReader2 = NULL;
        }
        delete pReader1;
        pReader1 = NULL;
            
        
        
    }
    
    if(pWriter != &std::cout)
        delete pWriter;
    
    if(!opt::discardFile.empty())
        delete pDiscardWriter;
    
    std::cerr << "\nPreprocess stats:\n";
    std::cerr << "Reads parsed:\t" << s_numReadsRead << "\n";
    std::cerr << "Reads kept:\t" << s_numReadsKept << " (" << (double)s_numReadsKept / (double)s_numReadsRead << ")\n";
    std::cerr << "Bases parsed:\t" << s_numBasesRead << "\n";
    std::cerr << "Bases kept:\t" << s_numBasesKept << " (" << (double)s_numBasesKept / (double)s_numBasesRead << ")\n";
    std::cerr << "Number of incorrectly paired reads that were discarded: " << s_numInvalidPE << "\n";
    delete pTimer;
    return 0;
}

// return true if the random value is lower than the acceptance value
bool samplePassTF()
{
    if(opt::sampleFreq >= 1.0f)
        return true; // no sampling
    
    double r = rand() / (RAND_MAX + 1.0f);
    return r < opt::sampleFreq;
}

//
// Handle command line arguments
//
void parseSampleReadsOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'o': arg >> opt::outFile; break;
            case 'd': arg >> opt::discardFile; break;
            case 'p': arg >> opt::sampleFreq; break;
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
    
    if(opt::peMode > 1)
    {
        std::cerr << SUBPROGRAM ": error pe-mode must be 0 or 1 (found: " << opt::peMode << ")\n";
        exit(EXIT_FAILURE);
    }
    
    

    
    opt::readFile = argv[optind++];
    
    
}
