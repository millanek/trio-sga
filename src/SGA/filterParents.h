//
//  filterParents.h
//  sga_git
//
//  Created by Milan Malinsky on 17/10/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//

#ifndef CORRECT_TRIO
#define CORRECT_TRIO
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int filterParentsMain(int argc, char** argv);

// options
void parseFilterParentsOptions(int argc, char** argv);

BWTIndexSet loadIndices(const std::string& readFile);
void deleteIndices(BWTIndexSet& indexSet);

#endif
