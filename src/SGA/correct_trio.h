//
//  correct_trio.h
//  sga_git
//
//  Created by Milan Malinsky on 14/01/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//
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
int correctTrioMain(int argc, char** argv);

// options
void parseCorrectTrioOptions(int argc, char** argv);

#endif
