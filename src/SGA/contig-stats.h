//
//  contig-stats.h
//  sga_git
//
//  Created by Milan Malinsky on 02/09/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//

#ifndef __sga_git__contig_stats__
#define __sga_git__contig_stats__

// Print an arbitrary vector to a file
template <class T> void print_vector(T vector, std::ofstream& outFile, char delim = '\t') {
    for (int i = 0; i < vector.size(); i++) {
        if (i == (vector.size()-1))
            outFile << vector[i] << std::endl;
        else
            outFile << vector[i] << delim;
    }
}


void parseContigStatsOptions(int argc, char** argv);
int contigStatsMain(int argc, char** argv);

#endif /* defined(__sga_git__contig_stats__) */
