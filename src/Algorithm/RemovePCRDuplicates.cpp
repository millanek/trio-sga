//
//  RemovePCRDuplicates.cpp
//  sga_git
//
//  Created by Milan Malinsky on 03/02/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "RemovePCRDuplicates.h"
#include "BWTAlgorithms.h"

//
//
//
PCRRemoveProcess::PCRRemoveProcess(PCRRemoveParameters params) : m_params(params)
{
    
}

//
PCRRemoveProcess::~PCRRemoveProcess()
{
    
}


//
PCRRemoveResult PCRRemoveProcess::process(const SequenceWorkItemPair& workItemPair)
{
    PCRRemoveResult result;
    result.bPCRDup = readDuplicateCheck(workItemPair); 
    return result;
}


bool PCRRemoveProcess::readDuplicateCheck(const SequenceWorkItemPair& workItemPair) {
    std::string seq1 = workItemPair.first.read.seq.toString();
    std::string seq2 = workItemPair.second.read.seq.toString();
    
    seq1 = seq1.substr(0,50);
    seq2 = seq2.substr(0,50);
    
    // Look up the BWT interval of the sequence (the same strand/orientation only)
    BWTInterval fwdInterval1 = BWTAlgorithms::findInterval(m_params.pBWT, seq1);
    BWTInterval fwdInterval2 = BWTAlgorithms::findInterval(m_params.pBWT, seq2);
    BWTAlgorithms::updateInterval(fwdInterval1, '$', m_params.pBWT);
    BWTAlgorithms::updateInterval(fwdInterval2, '$', m_params.pBWT);
    
    assert(fwdInterval1.isValid() && fwdInterval2.isValid());
    
    // We've got the lexicographic indexes for both reads. Turn it into a proper read IDs (in the original file).
    // Use the information that paired reads are next to each other (interleaved) in the fastq file
    std::vector<size_t> readIndices;
    for(int64_t i = fwdInterval1.lower; i <= fwdInterval1.upper; ++i) { 
        int64_t idx = m_params.pSAI->lookupLexoRank(i);
        // Only consider matching reads that are also first of a pair (the /1 sequence)
        if (idx % 2 == 0)
            readIndices.push_back(idx);
        if (idx % 2 == 1)
            readIndices.push_back(idx-1);
    }
    for(int64_t i = fwdInterval2.lower; i <= fwdInterval2.upper; ++i) { 
        int64_t idx = m_params.pSAI->lookupLexoRank(i);
        // Only consider matching reads that are also second in a pair (the /2 sequence)
        if (idx % 2 == 1) 
            readIndices.push_back(idx-1);
    }
    // Get the lowest read ID denoting this read pair 
    std::sort(readIndices.begin(),readIndices.end());
    std::vector<size_t>::iterator it = std::adjacent_find(readIndices.begin(), readIndices.end());
    assert(it != readIndices.end()); // There must be at least one read pair (the one currently being processed)
    size_t lowestReadIDforPair = (*it);
    
    // Check if the bit representing the lowest read ID for this read pair is set in the shared bit vector
    if(!m_params.pSharedBV->test(lowestReadIDforPair) )
    {
        // This read pair is not a duplicate
        
        // Attempt to atomically set the bit from false to true
        if(m_params.pSharedBV->updateCAS(lowestReadIDforPair, false, true))
        {
            // Call succeed, return that this read pair not a duplicate
            return false;
        }
        else
        {
            // Call failed, some other thread set the bit before
            // this thread. Return that the read pair is a duplicate
            return true;
        }
    }
    else
    {
        // this read pair is duplicate
        return true;
    }
    
}




PCRRemovePostProcess::PCRRemovePostProcess(std::ostream* pCorrectedWriter, 
                     std::ostream* pDiscardWriter, bool bPaired) : 
m_pCorrectedWriter(pCorrectedWriter),
m_pDiscardWriter(pDiscardWriter),
m_bPaired(bPaired),
m_PCRdupPassed(0), m_PCRdupFail(0)
{
        
}

//
PCRRemovePostProcess::~PCRRemovePostProcess()
{
    std::cout << "Reads passed PCR duplicate removal: " << m_PCRdupPassed << "\n";
    std::cout << "Reads failed PCR duplicate removal: " << m_PCRdupFail << "\n";
}


void PCRRemovePostProcess::process(const SequenceWorkItemPair& workItemPair, const PCRRemoveResult& result) {
    
    SeqRecord f_seq = workItemPair.first.read;
    SeqRecord s_seq = workItemPair.second.read;
    
    if (result.bPCRDup) {
        std::stringstream newID1; std::stringstream newID2;
        newID1 << f_seq.id << ",seqrank=" << workItemPair.first.idx;
        f_seq.id = newID1.str();
        newID2 << s_seq.id << ",seqrank=" << workItemPair.second.idx;
        s_seq.id = newID2.str();
        f_seq.write(*m_pDiscardWriter);
        s_seq.write(*m_pDiscardWriter);
        m_PCRdupFail = m_PCRdupFail + 2;
    } else {
        f_seq.write(*m_pCorrectedWriter);
        s_seq.write(*m_pCorrectedWriter);
        m_PCRdupPassed = m_PCRdupPassed + 2;
    }

    
}


