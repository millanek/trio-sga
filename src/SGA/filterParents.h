//
//  filterParents.h
//  sga_git
//
//  Created by Milan Malinsky on 17/10/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//

#ifndef FILTER_PARENTS
#define FILTER_PARENTS
#include <getopt.h>
#include "config.h"
#include "BWTAlgorithms.h"

// functions

//
int filterParentsMain(int argc, char** argv);

// options
void parseFilterParentsOptions(int argc, char** argv);


#endif
