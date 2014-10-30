//
//  removeUnpairedReads.cpp
//  sga_git
//
//  Created by Milan Malinsky on 29/10/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//

#include "removeUnpairedReads.h"
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include "Util.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "ASQG.h"
#include "gzstream.h"

// Functions

//
// Getopt
//
#define SUBPROGRAM "remove-unpaired-reads"
static const char *RUPR_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Milan Malinsky.\n"
"\n"
"Copyright 2014 Wellcome Trust Sanger Institute\n";

static const char *RUPR_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Remove reads that are not properly paired from READSFILE\n"
"\n"
"      --help                           display this help and exit\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static std::string prefix;
    static std::string readsFile;
}

static const char* shortopts = "p:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "prefix",             required_argument, NULL, 'p' },
    { "help",               no_argument,       NULL, OPT_HELP },
    { NULL, 0, NULL, 0 }
};


enum LastRead {
    READ_ONE,
    READ_TWO
};

//
// Main
//
int removeUPreadsMain(int argc, char** argv)
{
    parseRemoveUPreadsOptions(argc, argv);
    std::ifstream* readsFile = new std::ifstream(opt::readsFile.c_str());
    std::string line;
    std::vector<std::string> currentPair;
    int lineNo = 0;
    LastRead lr; bool bProperlyPaired = false;
    
    
    while (getline(*readsFile, line)) {
        if (lineNo % 4 == 0) {
            if (line[line.length()-1] == '1') {
                if (lr == READ_TWO && bProperlyPaired) {
                    assert(currentPair.size() >= 8);
                    for (int i = (int)currentPair.size()-8; i != currentPair.size(); i++) {
                        std::cout << currentPair[i] << std::endl;
                    }
                }
                currentPair.clear();
                lr = READ_ONE;
                bProperlyPaired = false;
            } else if (line[line.length()-1] == '2') {
                
                if (lr == READ_ONE) {
                    //std::cerr << line.substr(0,line.size()-1) << std::endl;
                    //std::cerr << currentPair[currentPair.size()-4].substr(0,currentPair[currentPair.size()-4].size()-1) << std::endl;
                    if (line.substr(0,line.size()-1) == currentPair[currentPair.size()-4].substr(0,currentPair[currentPair.size()-4].size()-1)) { // Try to match the read IDs
                        bProperlyPaired = true;
                    } else {
                        bProperlyPaired = false;
                    }
                } else {
                    if (bProperlyPaired) {
                        assert(currentPair.size() >= 8);
                        for (int i = (int)currentPair.size()-8; i != currentPair.size(); i++) {
                            std::cout << currentPair[i] << std::endl;
                        }
                        currentPair.clear();
                    }
                    bProperlyPaired = false;
                }
                lr = READ_TWO;
            } else {
                std::cerr << line << std::endl;
                assert(false);
            }
        }
        currentPair.push_back(line);
        lineNo++;
    }
    
    // Print the last pair
    if (bProperlyPaired) {
        assert(currentPair.size() >= 8);
        for (int i = (int)currentPair.size()-8; i != currentPair.size(); i++) {
            std::cout << currentPair[i] << std::endl;
        }
    }
    
    
    return 0;
}

//
// Handle command line arguments
//
void parseRemoveUPreadsOptions(int argc, char** argv)
{
    std::string algo_str;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'p': arg >> opt::prefix; break;
            case '?': die = true; break;
            case OPT_HELP:
                std::cout << RUPR_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << RUPR_VERSION_MESSAGE;
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
        std::cout << "\n" << RUPR_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::readsFile = argv[optind++];
    
    if(opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::readsFile);
    }
}
