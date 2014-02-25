//
//  MapBWTPairing.h
//  sga_git
//
//  Created by Milan Malinsky on 04/02/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef sga_git_MapBWTPairing_h
#define sga_git_MapBWTPairing_h
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "SampledSuffixArray.h"
#include "Util.h"


// Parameter object 
struct MapBWTPairingParameters
{
    //
    const BWT* pBWT;
    const SampledSuffixArray* pSAI;
    
};



class MapBWTPairingResult
{
public:
    MapBWTPairingResult() : bPCRDup(false) {}
    
    bool bPCRDup;
};


//
class MapBWTPairingProcess
{
public:
    ~MapBWTPairingProcess();
    MapBWTPairingProcess(const MapBWTPairingParameters params);
    
    MapBWTPairingResult process(const SequenceWorkItemPair& workItemPair);
    
    std::vector<int> readDuplicateCheck(const SequenceWorkItem& item); 
    
private:
    
    MapBWTPairingParameters m_params;
};

// Write the results from the overlap step to an ASQG file
class MapBWTPairingPostProcess
{
public:
    MapBWTPairingPostProcess(std::ostream* pCorrectedWriter, 
                         std::ostream* pDiscardWriter, bool bPaired);
    
    ~MapBWTPairingPostProcess();
    
    void process(const SequenceWorkItemPair& workItemPair, const MapBWTPairingResult& result);
    
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
