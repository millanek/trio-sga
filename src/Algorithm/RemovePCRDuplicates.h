//
//  RemovePCRDuplicates.h
//  sga_git
//
//  Created by Milan Malinsky on 03/02/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef sga_git_RemovePCRDuplicates_h
#define sga_git_RemovePCRDuplicates_h

#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "SampledSuffixArray.h"
#include "Util.h"
#include "BitVector.h"


// Parameter object 
struct PCRRemoveParameters
{
    //
    const BWT* pBWT;
    const SampledSuffixArray* pSAI;
    BitVector* pSharedBV;
    
};



class PCRRemoveResult
{
public:
    PCRRemoveResult() : bPCRDup(false) {}

    bool bPCRDup;
};


//
class PCRRemoveProcess
{
public:
    ~PCRRemoveProcess();
    PCRRemoveProcess(const PCRRemoveParameters params);
    
    PCRRemoveResult process(const SequenceWorkItemPair& workItemPair);
    
    bool readDuplicateCheck(const SequenceWorkItemPair& workItemPair); 
    
private:
    
    PCRRemoveParameters m_params;
};

// Write the results from the overlap step to an ASQG file
class PCRRemovePostProcess
{
public:
    PCRRemovePostProcess(std::ostream* pCorrectedWriter, 
                           std::ostream* pDiscardWriter, bool bPaired);
    
    ~PCRRemovePostProcess();
    
    void process(const SequenceWorkItemPair& workItemPair, const PCRRemoveResult& result);
    
private:
    
   // void collectMetrics(const std::string& originalSeq, 
    //                    const std::string& correctedSeq, const std::string& qualityStr);
    
    
    std::ostream* m_pCorrectedWriter;
    std::ostream* m_pDiscardWriter;
    bool m_bPaired; // Indicate if we are processing this read set as paired-end
    
    size_t m_PCRdupPassed;
    size_t m_PCRdupFail;
    
};


#endif
