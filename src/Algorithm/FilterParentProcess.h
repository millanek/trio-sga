//
//  FilterParentProcess.h
//  sga_git
//
//  Created by Milan Malinsky on 17/10/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//


// for a sequence work item
//
#ifndef FILTERPARENTPROCESS_H
#define FILTERPARENTPROCESS_H

#include "Util.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "Metrics.h"
#include "BWTIndexSet.h"
#include "SampledSuffixArray.h"

typedef std::map<std::string, int> KmerCountMap;


enum ECFlag
{
    ECF_NOTCORRECTED,
    ECF_CORRECTED,
    ECF_AMBIGIOUS,
    ECF_DUPLICATE
};


struct KmerCachesFilter
{
    KmerCountMap kmerCacheOffspring;
    KmerCountMap kmerCacheParent;
};

enum OffspringStatus
{
    OFFSPRING_ABOVE_THRESHOLD,
    OFFSPRING_ONCE,
    OFFSPRING_ABSENT
};

class FilterParentResult
{
public:
    FilterParentResult() : kmerQC(false),os(OFFSPRING_ABSENT) {}
    
    DNAString correctSequence;
    ECFlag flag;
    
    // Metrics
    bool kmerQC;
    OffspringStatus os;
};

struct FilterParentPairResult
{
    FilterParentResult firstResult;
    FilterParentResult secondResult;
    OffspringStatus os;
};


// Parameter object for the error corrector
struct FilterParentParameters
{
    //
    BWTIndexSet offspringIndices;
    BWTIndexSet parentIndices;
    
    bool bCorrect;
    // k-mer based corrector params
    int numKmerRounds;
    int kmerLength;
    size_t parentThreshold;
    size_t offspringThreshold;
};

//
class FilterParentProcess
{
public:
    ~FilterParentProcess();
    FilterParentProcess(const FilterParentParameters params);
    
    FilterParentResult process(const SequenceWorkItem& item);
    FilterParentPairResult process(const SequenceWorkItemPair& workItemPair);
    
    FilterParentResult kmerCorrection(const SequenceWorkItem& item, KmerCachesFilter& kmerCaches);
    
private:
    
    bool attemptKmerCorrection(size_t i, size_t k_idx, size_t minCount, std::string& readSequence);
    int getKmerCount(const std::string& kmer, const BWTIndexSet& indexSet, KmerCountMap& kmerCache);
    int maximum(int x, int y, int z);
    OffspringStatus checkInOffspring(const DNAString& correctSequence, KmerCachesFilter& kmerCaches);
    FilterParentParameters m_params;
};

// Write the results from the overlap step to an ASQG file
class FilterParentPostProcess
{
public:
    FilterParentPostProcess(std::ostream* pCorrectedWriter, std::ostream* pDiscardWriter, bool bCollectMetrics, bool bPaired);
    
    ~FilterParentPostProcess();
    
    void process(const SequenceWorkItem& item, const FilterParentResult& result);
    void process(const SequenceWorkItemPair& workItemPair, const FilterParentPairResult& results);
    void writeMetrics(std::ostream* pWriter);
    
private:
    
    void collectMetrics(const std::string& originalSeq,
                        const std::string& correctedSeq, const std::string& qualityStr);
    
    // Helper functions for writing phased paired-end data
    void writeSingleRead(const FilterParentResult& result, const SeqRecord& seq);
    void writeBothReads(const FilterParentResult& result, const SeqRecord& f_seq, const SeqRecord& s_seq);
    
    std::ostream* m_pCorrectedWriter;
    std::ostream* m_pDiscardWriter;
    bool m_bCollectMetrics;
    bool m_bPaired; // Indicate if we are processing this read set as paired-end
    
    ErrorCountMap<char> m_qualityMetrics;
    ErrorCountMap<int64_t> m_positionMetrics;
    ErrorCountMap<char> m_originalBaseMetrics;
    ErrorCountMap<std::string> m_precedingSeqMetrics;
    
    size_t m_totalBases;
    size_t m_totalErrors;
    size_t m_readsKept;
    size_t m_readsDiscarded;
    
    size_t m_kmerQCPassed;
    size_t m_qcFail;

    
};

#endif


