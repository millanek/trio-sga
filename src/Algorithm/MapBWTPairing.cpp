//
//  MapBWTPairing.cpp
//  sga_git
//
//  Created by Milan Malinsky on 04/02/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "MapBWTPairing.h"
#include "BWTAlgorithms.h"

//
//
//
MapBWTPairingProcess::MapBWTPairingProcess(MapBWTPairingParameters params) : m_params(params)
{
    
}

//
MapBWTPairingProcess::~MapBWTPairingProcess()
{
    
}


//
MapBWTPairingResult MapBWTPairingProcess::process(const SequenceWorkItemPair& workItemPair)
{
    MapBWTPairingResult result;
    
    
    std::vector<int> dupCheckR1 = readDuplicateCheck(workItemPair.first); 
    std::vector<int> dupCheckR2; 
    
    if (dupCheckR1.empty()) { 
        std::cout << "This is not a PCR duplicate pair...Shall I check for duplicates of the other read in the pair?\n"; 
    } else {
        dupCheckR2 = readDuplicateCheck(workItemPair.second);
    } 
    
    
    
    
    /*if(result.kmerPassed && result.dupPassed && m_params.checkHPRuns)
     result.hpPassed = performHomopolymerCheck(workItem);
     else
     result.hpPassed = true;
     
     if(m_params.checkDegenerate && result.dupPassed && result.kmerPassed && result.hpPassed)
     result.degenPassed = performDegenerateCheck(workItem);
     else
     result.degenPassed = true; */
    
    return result;
}


std::vector<int> MapBWTPairingProcess::readDuplicateCheck(const SequenceWorkItem& workItem) {
    std::vector<int> result;
    std::string seq = workItem.read.seq.toString();
    
    // Look up the BWT interval of the sequence (the same strand/orientation only)
    BWTInterval fwdInterval = BWTAlgorithms::findInterval(m_params.pBWT, seq);
    BWTAlgorithms::updateInterval(fwdInterval, '$', m_params.pBWT);
    
    assert(fwdInterval.isValid());
    
    if (fwdInterval.size() == 1) {
        // This read hasn't got any duplicates
        return result;
    } else {
        // We've got the lexicographic index(es) for duplicate reads. Turn it into a proper ID(s)
        for(int64_t i = fwdInterval.lower; i <= fwdInterval.upper; ++i) { 
            int64_t idx = m_params.pSAI->lookupLexoRank(i);
            result.push_back(idx);
        }
    }
    
    return result;
    
}




MapBWTPairingPostProcess::MapBWTPairingPostProcess(std::ostream* pCorrectedWriter, 
                                           std::ostream* pDiscardWriter, bool bPaired) : 
m_pCorrectedWriter(pCorrectedWriter),
m_pDiscardWriter(pDiscardWriter),
m_bPaired(bPaired),
m_PCRdupPassed(0), m_PCRdupFail(0)
{
    
}


void MapBWTPairingPostProcess::process(const SequenceWorkItemPair& workItemPair, const MapBWTPairingResult& result) {
    
    if (result.bPCRDup) {
        m_PCRdupFail++;
    } else {
        m_PCRdupPassed++;
    }
    SeqRecord f_seq = workItemPair.first.read;
    SeqRecord s_seq = workItemPair.second.read;
    
}
