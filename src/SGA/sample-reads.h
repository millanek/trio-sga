//
//  sample-reads.h
//  sga_git
//
//  Created by Milan Malinsky on 19/01/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef sga_git_sample_reads_h
#define sga_git_sample_reads_h
#include <getopt.h>
#include "config.h"

void parseSampleReadsOptions(int argc, char** argv);
int sampleReadsMain(int argc, char** argv);
bool samplePassTF();

#endif
