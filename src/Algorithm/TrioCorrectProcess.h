//
//  TrioCorrectProcess.h
//  sga_git
//
//  Created by Milan Malinsky on 15/01/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

// ErrorCorrectProcess - Wrapper to perform error correction
// for a sequence work item
//
#ifndef TRIOCORRECTPROCESS_H
#define TRIOCORRECTPROCESS_H

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

enum ParentCheckResult
{
    PARENT_ABOVE_THRESHOLD,
    PARENT_ONCE,
    PARENT_ABSENT
};

struct KmerCaches
{
    KmerCountMap kmerCacheOffspring;
    KmerCountMap kmerCacheFather;
    KmerCountMap kmerCacheMother;
};

class TrioCorrectResult
{
public:
    TrioCorrectResult() : kmerQC(false),inFather(PARENT_ABSENT),inMother(PARENT_ABSENT) {}
    
    DNAString correctSequence;
    ECFlag flag;
    
    // Metrics
    bool kmerQC;
    ParentCheckResult inFather;
    ParentCheckResult inMother;
};

struct TrioCorrectPairResult
{
    TrioCorrectResult firstResult;
    TrioCorrectResult secondResult;
};


// Parameter object for the error corrector
struct TrioCorrectParameters
{
    //
    BWTIndexSet offspringIndices;
    BWTIndexSet motherIndices;
    BWTIndexSet fatherIndices;
    
    // k-mer based corrector params
    int numKmerRounds;
    int kmerLength;
    size_t offspringThreshold;
    size_t offspringConditionalThreshold;
    size_t motherThreshold;
    size_t fatherThreshold;
    
};


struct ParentStatus
{
    ParentCheckResult inFather;
    ParentCheckResult inMother;
};


struct ParentCorrectionStatus
{
    bool correctedInFather;
    bool correctedInMother;
};

//
class TrioCorrectProcess
{
public:
    ~TrioCorrectProcess();
    TrioCorrectProcess(const TrioCorrectParameters params);
    
    TrioCorrectResult process(const SequenceWorkItem& item);
    TrioCorrectPairResult process(const SequenceWorkItemPair& workItemPair);
    
    TrioCorrectResult kmerCorrection(const SequenceWorkItem& item, KmerCaches& kmerCaches, ParentCorrectionStatus& parentCorrectionStatus);
    void phase(TrioCorrectResult& result, KmerCaches& kmerCaches); 
    void phase(TrioCorrectPairResult& result, KmerCaches& kmerCaches);  
    
private:
    
    bool attemptKmerCorrection(size_t i, size_t k_idx, size_t minCount, std::string& readSequence, ParentCorrectionStatus& pCS);
    int getKmerCount(const std::string& kmer, const BWTIndexSet& indexSet, KmerCountMap& kmerCache);
    int maximum(int x, int y, int z);
    ParentStatus checkInParents(const DNAString& correctSequence, KmerCaches& kmerCaches);
    TrioCorrectParameters m_params;
};

// Write the results from the overlap step to an ASQG file
class TrioCorrectPostProcess
{
public:
    TrioCorrectPostProcess(std::ostream* pCorrectedWriter, 
                            std::ostream* pDiscardWriter, bool bCollectMetrics, bool bPaired);
    TrioCorrectPostProcess(std::ostream* pMaternalWriter, std::ostream* pPaternalWriter,
                           std::ostream* pDiscardWriter, std::ostream* pNeitherParentWriter, std::ostream* pInconsistentPhaseWriter, 
                           bool bCollectMetrics, bool bPaired);
    
    ~TrioCorrectPostProcess();
    
    void process(const SequenceWorkItem& item, const TrioCorrectResult& result);
    void process(const SequenceWorkItemPair& workItemPair, const TrioCorrectPairResult& results);
    void writeMetrics(std::ostream* pWriter);
    
private:
    
    void collectMetrics(const std::string& originalSeq, 
                        const std::string& correctedSeq, const std::string& qualityStr);
    
    // Helper functions for writing phased paired-end data
    void writeMaternalPaternal(const SeqRecord& f_seq, const SeqRecord& s_seq);
    void writeSingleReadBasedOnPhase(const TrioCorrectResult& result, const SeqRecord& seq);
    void writeBothReadsBasedOnPhaseOfOne(const TrioCorrectResult& result, const SeqRecord& f_seq, const SeqRecord& s_seq);
    void findIfOnceInEitherParent(const TrioCorrectResult& result);
    bool checkInconsistentPhase(const TrioCorrectPairResult& results);
    
    std::ostream* m_pCorrectedWriter;
    std::ostream* m_pMaternalWriter;
    std::ostream* m_pPaternalWriter;
    std::ostream* m_pDiscardWriter;
    std::ostream* m_pNeitherParentWriter;   // Output sequences that passed QC but are present in neither parent
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

