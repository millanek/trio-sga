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
    FilterParentResult() : kmerQC(false),inFather(PARENT_ABSENT),inMother(PARENT_ABSENT) {}
    
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
    FilterParentPostProcess(std::ostream* pCorrectedWriter,
                           std::ostream* pDiscardWriter, bool bCollectMetrics, bool bPaired);
    FilterParentPostProcess(std::ostream* pMaternalWriter, std::ostream* pPaternalWriter,
                           std::ostream* pDiscardWriter, std::ostream* pNeitherParentWriter, std::ostream* pInconsistentPhaseWriter,
                           bool bCollectMetrics, bool bPaired);
    
    ~TrioCorrectPostProcess();
    
    void process(const SequenceWorkItem& item, const FilterParentResult& result);
    void process(const SequenceWorkItemPair& workItemPair, const FilterParentPairResult& results);
    void writeMetrics(std::ostream* pWriter);
    
private:
    
    void collectMetrics(const std::string& originalSeq,
                        const std::string& correctedSeq, const std::string& qualityStr);
    
    // Helper functions for writing phased paired-end data
    void writeSingleRead(const FilterParentResult& result, const SeqRecord& seq);
    void writeBothReadsBasedOnOne(const TrioCorrectResult& result, const SeqRecord& f_seq, const SeqRecord& s_seq);
    
    std::ostream* m_pCorrectedWriter;
    std::ostream* m_pDiscardWriter;
    std::ostream* m_pInconsistentPhaseWriter; // Output sequences with inconsistent phasing (read 1 in parent A, read 2 in parent B)
    bool m_bCollectMetrics;
    bool m_bPaired; // Indicate if we are processing this read set as paired-end
    bool m_bPhase;
    
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
    size_t m_reads_paternal;
    size_t m_reads_maternal;
    size_t m_passed_QC_but_neither_m_or_f; // Reads that passed offspring QC but are neither in the mother or the father
    size_t m_passed_QC_and_once_in_m_or_f;
    size_t m_passed_QC_but_inconsistent_phase; // Reads that passed offspring QC, are below threshold in the mother and the father, but al least once in the mother or the father
    
};

#endif


