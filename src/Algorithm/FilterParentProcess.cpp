//
//  FilterParentProcess.cpp
//  sga_git
//
//  Created by Milan Malinsky on 17/10/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//

#include "FilterParentProcess.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include "multiple_alignment.h"
#include "KmerOverlaps.h"
#include "StringThreader.h"
#include "BWTAlgorithms.h"

//#define KMER_TESTING 1

//
//
//
FilterParentProcess::FilterParentProcess(const FilterParentParameters params) : m_params(params)
{
    
}

//
FilterParentProcess::~FilterParentProcess()
{
    
}

//
FilterParentResult FilterParentProcess::process(const SequenceWorkItem& workItem)
{
    KmerCachesFilter kmerCaches;
    FilterParentResult result = kmerCorrection(workItem,kmerCaches);
    result.os = checkInOffspring(result.correctSequence, kmerCaches);
    kmerCaches.kmerCacheOffspring.clear();
    kmerCaches.kmerCacheParent.clear();
    //if(!result.kmerQC)
    //    std::cout << workItem.read.id << " failed error correction QC\n";
    return result;
}

FilterParentPairResult FilterParentProcess::process(const SequenceWorkItemPair& workItemPair)
{
    KmerCachesFilter kmerCaches;
    FilterParentPairResult result;
    result.firstResult = kmerCorrection(workItemPair.first,kmerCaches);
    result.secondResult = kmerCorrection(workItemPair.second,kmerCaches);
    result.firstResult.os = checkInOffspring(result.firstResult.correctSequence, kmerCaches);
    result.secondResult.os = checkInOffspring(result.secondResult.correctSequence, kmerCaches);
    if (result.firstResult.os == OFFSPRING_ABSENT || result.secondResult.os == OFFSPRING_ABSENT) {
        result.os = OFFSPRING_ABSENT;
    } else if (result.firstResult.os == OFFSPRING_ONCE || result.secondResult.os == OFFSPRING_ONCE) {
        result.os = OFFSPRING_ONCE;
    } else {
        result.os = OFFSPRING_ABOVE_THRESHOLD;
    }
    
    kmerCaches.kmerCacheOffspring.clear();
    kmerCaches.kmerCacheParent.clear();
    
    return result;
}

OffspringStatus FilterParentProcess::checkInOffspring(const DNAString& correctSequence, KmerCachesFilter& kmerCaches) {
    OffspringStatus offspringStatus;
    std::string seq = correctSequence.toString();
    
    int n = (int)seq.size();
    int nk = n - m_params.kmerLength + 1;
    
    std::vector<int> thresholdSolidVector(n, 0);
    std::vector<int> onceSolidVector(n, 0);
    
    for(int i = 0; i < nk; ++i)
    {
        std::string kmer = seq.substr(i, m_params.kmerLength);
        
        size_t countOffspring = getKmerCount(kmer, m_params.offspringIndices, kmerCaches.kmerCacheOffspring);
        
        if(countOffspring >= m_params.offspringThreshold) {
            for(int j = i; j < i + m_params.kmerLength; ++j)
                thresholdSolidVector[j] = 1;
        } else if (countOffspring >= 1) {
            for(int j = i; j < i + m_params.kmerLength; ++j)
                onceSolidVector[j] = 1;
        }
    }
    offspringStatus = OFFSPRING_ABOVE_THRESHOLD;
    for(int i = 0; i < n; ++i) {
        if(thresholdSolidVector[i] != 1)
            offspringStatus = OFFSPRING_ONCE;
    }
    
    if (offspringStatus == OFFSPRING_ONCE) {
        for(int i = 0; i < n; ++i) {
            if(onceSolidVector[i] != 1)
                offspringStatus = OFFSPRING_ABSENT;
        }
    }
    
    return offspringStatus;
}


// Correct a read with a k-mer based corrector
FilterParentResult FilterParentProcess::kmerCorrection(const SequenceWorkItem& workItem, KmerCachesFilter& kmerCaches)
{
    assert(m_params.parentIndices.pBWT != NULL);
    assert(m_params.parentIndices.pCache != NULL);
    
    FilterParentResult result;
    
    // KmerCountMap kmerCacheOffspring;
    // KmerCountMap kmerCacheFather;
    // KmerCountMap kmerCacheMother;
    
    SeqRecord currRead = workItem.read;
    std::string readSequence = workItem.read.seq.toString();
    
#ifdef KMER_TESTING
    std::cout << "Kmer correcting read " << workItem.read.id << "\n";
#endif
    
    if((int)readSequence.size() < m_params.kmerLength)
    {
        // The read is shorter than the kmer length, nothing can be done
        result.correctSequence = readSequence;
        result.kmerQC = false;
        return result;
    }
    
    int n = (int)readSequence.size();
    int nk = n - m_params.kmerLength + 1;
    
    // Are all kmers in the read well-represented?
    bool allSolid = false;
    bool done = false;
    int rounds = 0;
    int maxAttempts = m_params.numKmerRounds;
    
    // For each kmer, calculate the minimum phred score seen in the bases
    // of the kmer
    std::vector<int> minPhredVector(nk, 0);
    for(int i = 0; i < nk; ++i)
    {
        int end = i + m_params.kmerLength - 1;
        int minPhred = std::numeric_limits<int>::max();
        for(int j = i; j <= end; ++j)
        {
            int ps = workItem.read.getPhredScore(j);
            if(ps < minPhred)
                minPhred = ps;
        }
        minPhredVector[i] = minPhred;
    }
    
    while(!done && nk > 0)
    {
        // Compute the kmer counts across the read
        // and determine the positions in the read that are not covered by any solid kmers
        // These are the candidate incorrect bases
        std::vector<int> countVectorParent(nk, 0);
        std::vector<int> solidVector(n, 0);
        
        for(int i = 0; i < nk; ++i)
        {
            std::string kmer = readSequence.substr(i, m_params.kmerLength);
            
            //size_t countOffspring = getKmerCount(kmer, m_params.offspringIndices, kmerCaches.kmerCacheOffspring);
            size_t countParent = getKmerCount(kmer, m_params.parentIndices, kmerCaches.kmerCacheParent);
            
            // Get the phred score for the last base of the kmer
            int phred = minPhredVector[i];
            countVectorParent[i] = (int)countParent;
            //            std::cout << i << "\t" << phred << "\t" << count << "\n";
            
            // Adjust offspringThreshold based on phred scores
            unsigned int parentPhredBasedThreshold;
            if (phred >= 20) {
                parentPhredBasedThreshold = (int)m_params.parentThreshold;
            } else {
                parentPhredBasedThreshold = (int)m_params.parentThreshold + 1;
            }
            if(countParent >= parentPhredBasedThreshold)
            {
                for(int j = i; j < i + m_params.kmerLength; ++j)
                    solidVector[j] = 1;
            }
#ifdef KMER_TESTING
            std::cout << "Count parent: " << countParent << " Count offspring: " << countOffspring << "\n";
#endif
        }
        
        allSolid = true;
        for(int i = 0; i < n; ++i)
        {
#ifdef KMER_TESTING
            std::cout << "Position[" << i << "] = " << solidVector[i] << "\n";
#endif
            if(solidVector[i] != 1)
                allSolid = false;
        }
        
#ifdef KMER_TESTING
        std::cout << "Read " << workItem.read.id << (allSolid ? " is solid\n" : " has potential errors\n");
#endif
        
        // Stop if all kmers are well represented or we have exceeded the number of correction rounds
        if(allSolid || rounds++ > maxAttempts)
            break;
        
        // Attempt to correct the leftmost potentially incorrect base
        bool corrected = false;
        for(int i = 0; i < n; ++i)
        {
            if(solidVector[i] != 1)
            {
                // Attempt to correct the base using the leftmost covering kmer
                int phred = workItem.read.getPhredScore(i);
                int parentPhredBasedThreshold;
                if (phred >= 20) {
                    parentPhredBasedThreshold = (int)m_params.parentThreshold;
                } else {
                    parentPhredBasedThreshold = (int)m_params.parentThreshold + 1;
                }
                
                int left_k_idx = (i + 1 >= m_params.kmerLength ? i + 1 - m_params.kmerLength : 0);
                corrected = attemptKmerCorrection(i, left_k_idx, std::max(countVectorParent[left_k_idx], parentPhredBasedThreshold), readSequence);
                if(corrected)
                    break;
                
                // base was not corrected, try using the rightmost covering kmer
                size_t right_k_idx = std::min(i, n - m_params.kmerLength);
                corrected = attemptKmerCorrection(i, right_k_idx, std::max(countVectorParent[right_k_idx], parentPhredBasedThreshold), readSequence);
                if(corrected)
                    break;
            }
        }
        
        // If no base in the read was corrected, stop the correction process
        if(!corrected)
        {
            assert(!allSolid);
            done = true;
        }
    }
    
    if(allSolid)
    {
        result.correctSequence = readSequence;
        result.kmerQC = true;
    }
    else
    {
        result.correctSequence = workItem.read.seq.toString();
        result.kmerQC = false;
    }
    
    return result;
}

int FilterParentProcess::getKmerCount(const std::string& kmer, const BWTIndexSet& indexSet, KmerCountMap& kmerCache) {
    // First check if this kmer is in the cache
    // If its not, find its count from the fm-index and cache it
    int count = 0;
    KmerCountMap::iterator iter = kmerCache.find(kmer);
    if(iter != kmerCache.end())
    {
        count = iter->second;
    }
    else
    {
        count = (int)BWTAlgorithms::countSequenceOccurrences(kmer, indexSet);
        kmerCache.insert(std::make_pair(kmer, count));
    }
    return count;
}


// Attempt to correct the base at position idx in readSequence. Returns true if a correction was made
// The correction is made only if the count of the corrected kmer is at least minCount
bool FilterParentProcess::attemptKmerCorrection(size_t i, size_t k_idx, size_t minCount, std::string& readSequence)
{
    assert(i >= k_idx && i < k_idx + m_params.kmerLength);
    size_t base_idx = i - k_idx;
    char originalBase = readSequence[i];
    std::string kmer = readSequence.substr(k_idx, m_params.kmerLength);
    size_t bestCount = 0;
    char bestBase = '$';
    std::string local_phase = "";
    
#if KMER_TESTING
    std::cout << "i: " << i << " k-idx: " << k_idx << " " << kmer << " " << reverseComplement(kmer) << "\n";
#endif
    
    
    for(int j = 0; j < DNA_ALPHABET::size; ++j)
    {
        char currBase = ALPHABET[j];
        if(currBase == originalBase)
            continue;
        kmer[base_idx] = currBase;
        size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.parentIndices);
        
#if KMER_TESTING
        printf("base: %c offspring: %zu mother: %zu father: %zu\n", currBase, count, countMother, countFather);
#endif
        if(count > bestCount && count >= minCount)
        {
            // Multiple corrections exist, do not correct
            if(bestBase != '$')
                return false;
            
            bestCount = count;
            bestBase = currBase;
        }
    }
    
    if(bestCount >= minCount)
    {
        assert(bestBase != '$');
        readSequence[i] = bestBase;
        return true;
    }
    return false;
}




int FilterParentProcess::maximum(int x, int y, int z) {
    int max = x;
    
    if (y > max) {
        max = y;
    }
    
    if (z > max) {
        max = z;
    }
    
    return max;
}

//
//
//
FilterParentPostProcess::FilterParentPostProcess(std::ostream* pCorrectedWriter,
                                               std::ostream* pDiscardWriter,
                                               bool bCollectMetrics, bool bPaired) :
m_pCorrectedWriter(pCorrectedWriter),
m_pDiscardWriter(pDiscardWriter),
m_bCollectMetrics(bCollectMetrics),
m_bPaired(bPaired),
m_bPhase(false),
m_totalBases(0), m_totalErrors(0),
m_readsKept(0), m_readsDiscarded(0),
m_kmerQCPassed(0),
m_qcFail(0),
m_reads_paternal(0), m_reads_maternal(0),
m_passed_QC_but_neither_m_or_f(0),
m_passed_QC_and_once_in_m_or_f(0),
m_passed_QC_but_inconsistent_phase(0)
{
    
}

FilterParentPostProcess::FilterParentPostProcess(std::ostream* pMaternalWriter, std::ostream* pPaternalWriter,
                                               std::ostream* pDiscardWriter,std::ostream* pNeitherParentWriter, std::ostream* pInconsistentPhaseWriter, bool bCollectMetrics,
                                               bool bPaired) :
m_pCorrectedWriter(NULL),
m_pDiscardWriter(pDiscardWriter),
m_pInconsistentPhaseWriter(pInconsistentPhaseWriter),
m_bCollectMetrics(bCollectMetrics),
m_bPaired(bPaired),
m_bPhase(true),
m_totalBases(0), m_totalErrors(0),
m_readsKept(0), m_readsDiscarded(0),
m_kmerQCPassed(0),
m_qcFail(0),
m_reads_paternal(0), m_reads_maternal(0),
m_passed_QC_but_neither_m_or_f(0),
m_passed_QC_and_once_in_m_or_f(0),
m_passed_QC_but_inconsistent_phase(0)
{
    
}




//
void FilterParentPostProcess::writeMetrics(std::ostream* pWriter)
{
    *pWriter << "ErrorCorrect -- Corrected " << m_totalErrors << " out of " << m_totalBases <<
    " bases (" << (double)m_totalErrors / m_totalBases << ")\n";
    *pWriter << "Reads that passed kmer QC: " << m_kmerQCPassed << "\n";
    *pWriter << "Reads that failed kmer QC: " << m_qcFail << "\n";
    *pWriter << "Kept " << m_readsKept << " reads. Discarded " << m_readsDiscarded <<
    " reads (" << (double)m_readsDiscarded / (m_readsKept + m_readsDiscarded)<< ")\n";
    if (m_bPhase) {
        *pWriter << "Reads that passed kmer QC but are neither in mother or father data: " << m_passed_QC_but_neither_m_or_f << "\n";
        *pWriter << "Reads that passed kmer QC but the phasing is inconsistent: " << m_passed_QC_but_inconsistent_phase << "\n";
        *pWriter<< "Reads that passed kmer QC, are below threshold in both parents, but at least once in the mother or the father: " << m_passed_QC_and_once_in_m_or_f << "\n";
        *pWriter << "Reads for paternal chromosome assembly: " << m_reads_paternal << "(out of " << m_readsKept << ")\n";
        *pWriter << "Reads for maternal chromosome assembly: " << m_reads_maternal << "(out of " << m_readsKept << ")\n";
    }
    m_positionMetrics.write(pWriter, "Bases corrected by position\n", "pos");
    m_originalBaseMetrics.write(pWriter, "\nOriginal base that was corrected\n", "base");
    m_precedingSeqMetrics.write(pWriter, "\nkmer preceding the corrected base\n", "kmer");
    m_qualityMetrics.write(pWriter, "\nBases corrected by quality value\n\n", "quality");
    
    std::cout << "ErrorCorrect -- Corrected " << m_totalErrors << " out of " << m_totalBases <<
    " bases (" << (double)m_totalErrors / m_totalBases << ")\n";
    std::cout << "Reads that passed kmer QC: " << m_kmerQCPassed << "\n";
    std::cout << "Reads that failed kmer QC: " << m_qcFail << "\n";
    std::cout << "Kept " << m_readsKept << " reads. Discarded " << m_readsDiscarded <<
    " reads (" << (double)m_readsDiscarded / (m_readsKept + m_readsDiscarded)<< ")\n";
    if (m_bPhase) {
        std::cout << "Reads that passed kmer QC but are neither in mother or father data: " << m_passed_QC_but_neither_m_or_f << "\n";
        std::cout << "Reads that passed kmer QC but the phasing is inconsistent: " << m_passed_QC_but_inconsistent_phase << "\n";
        std::cout << "Reads that passed kmer QC, are below threshold in both parents, but at least once in the mother or the father: " << m_passed_QC_and_once_in_m_or_f << "\n";
        std::cout << "Reads for paternal chromosome assembly: " << m_reads_paternal << "(out of " << m_readsKept << ")\n";
        std::cout << "Reads for maternal chromosome assembly: " << m_reads_maternal << "(out of " << m_readsKept << ")\n";
    }
}

//
void FilterParentPostProcess::process(const SequenceWorkItem& item, const FilterParentResult& result)
{
    
    // Determine if the read should be discarded
    bool readQCPass = true;
    if(result.kmerQC)
    {
        m_kmerQCPassed += 1;
    }
    else
    {
        readQCPass = false;
        m_qcFail += 1;
    }
    
    // Collect metrics for the reads that were actually corrected
    if(m_bCollectMetrics && readQCPass)
    {
        collectMetrics(item.read.seq.toString(),
                       result.correctSequence.toString(),
                       item.read.qual);
    }
    
    SeqRecord record = item.read;
    record.seq = result.correctSequence;
    
    if (m_bPhase) {
        if (readQCPass) {
            ++m_readsKept;
            bool neither_m_or_f = true; // Passed QC but maybe neither in mother or the father --- I should count these
            if (result.inMother) {
                record.write(*m_pMaternalWriter);
                ++m_reads_maternal;
                neither_m_or_f = false;
            }
            if (result.inFather) {
                record.write(*m_pPaternalWriter);
                ++m_reads_paternal;
                neither_m_or_f = false;
            }
            if (neither_m_or_f) {
                ++m_passed_QC_but_neither_m_or_f;
                // We do not discard anything that passed kmerQC (these could be de-novo mutations in the offspring)
                record.write(*m_pMaternalWriter);
                record.write(*m_pPaternalWriter);
                ++m_reads_maternal;
                ++m_reads_paternal;
            }
        } else {
            if (m_pDiscardWriter == NULL) {
                record.write(*m_pMaternalWriter);
                record.write(*m_pPaternalWriter);
                ++m_reads_maternal;
                ++m_reads_paternal;
                ++m_readsKept;
            } else {
                record.write(*m_pDiscardWriter);
                ++m_readsDiscarded;
            }
        }
    } else {
        if(readQCPass || m_pDiscardWriter == NULL)
        {
            record.write(*m_pCorrectedWriter);
            ++m_readsKept;
        }
        else
        {
            record.write(*m_pDiscardWriter);
            ++m_readsDiscarded;
        }
    }
}

void FilterParentPostProcess::process(const SequenceWorkItemPair& workItemPair, const FilterParentPairResult& results)
{
    if (results.firstResult.kmerQC) { m_kmerQCPassed += 1; }
    else { m_qcFail += 1; }
    if (results.secondResult.kmerQC) { m_kmerQCPassed += 1; }
    else { m_qcFail += 1; }
    
    // Collect metrics for the reads that were actually corrected
    if (m_bCollectMetrics && results.firstResult.kmerQC) {
        collectMetrics(workItemPair.first.read.seq.toString(), results.firstResult.correctSequence.toString(),
                       workItemPair.first.read.qual);
    }
    if (m_bCollectMetrics && results.secondResult.kmerQC) {
        collectMetrics(workItemPair.second.read.seq.toString(), results.secondResult.correctSequence.toString(),
                       workItemPair.second.read.qual);
    }
    
    SeqRecord f_seq = workItemPair.first.read;
    SeqRecord s_seq = workItemPair.second.read;
    f_seq.seq = results.firstResult.correctSequence;
    s_seq.seq = results.secondResult.correctSequence;
    
    if (m_bPhase) {
        if((results.firstResult.kmerQC && results.secondResult.kmerQC) || m_pDiscardWriter == NULL) {
            m_readsKept = m_readsKept + 2;
            if (results.firstResult.kmerQC && results.secondResult.kmerQC) {
                bool pair_neither_m_or_f = true; // Both reads passed QC but maybe the fragment appears to be neither from the mother nor the father
                if ((results.firstResult.inFather == PARENT_ABOVE_THRESHOLD) && (results.secondResult.inFather == PARENT_ABOVE_THRESHOLD)) {
                    f_seq.write(*m_pPaternalWriter); s_seq.write(*m_pPaternalWriter);
                    m_reads_paternal = m_reads_paternal + 2;
                    pair_neither_m_or_f = false;
                }
                if ((results.firstResult.inMother == PARENT_ABOVE_THRESHOLD) && (results.secondResult.inMother == PARENT_ABOVE_THRESHOLD)) {
                    f_seq.write(*m_pMaternalWriter); s_seq.write(*m_pMaternalWriter);
                    m_reads_maternal = m_reads_maternal + 2;
                    pair_neither_m_or_f = false;
                }
                if (pair_neither_m_or_f) {
                    if ((results.firstResult.inFather == PARENT_ABSENT) && (results.firstResult.inMother == PARENT_ABSENT)) {
                        ++m_passed_QC_but_neither_m_or_f;
                        f_seq.write(*m_pNeitherParentWriter);
                    }
                    if ((results.secondResult.inFather == PARENT_ABSENT) && (results.secondResult.inMother == PARENT_ABSENT)) {
                        ++m_passed_QC_but_neither_m_or_f;
                        s_seq.write(*m_pNeitherParentWriter);
                    }
                    findIfOnceInEitherParent(results.firstResult);
                    findIfOnceInEitherParent(results.secondResult);
                    // We discard pairs with inconsistent phase
                    // Anything else that passed kmerQC could be de-novo mutations in the offspring, so we keep it
                    if (checkInconsistentPhase(results)) {
                        m_passed_QC_but_inconsistent_phase = m_passed_QC_but_inconsistent_phase + 2;
                        f_seq.write(*m_pInconsistentPhaseWriter); s_seq.write(*m_pInconsistentPhaseWriter);
                        f_seq.write(*m_pDiscardWriter); s_seq.write(*m_pDiscardWriter);
                    }
                    else {
                        writeMaternalPaternal(f_seq, s_seq);
                    }
                }
            } else {
                // Not passing kmer QC implies that at least one of the read pairs failed to pass correction threshold in all
                // three datasets (offspring, mother, father) but we do not discard them
                if (results.firstResult.kmerQC) {  // First one passed
                    writeBothReadsBasedOnPhaseOfOne(results.firstResult, f_seq, s_seq);
                } else if (results.secondResult.kmerQC) {  // Second one passed - so phase the pair based on it
                    writeBothReadsBasedOnPhaseOfOne(results.secondResult, f_seq, s_seq);
                } else {  // Both reads failed QC, so just write them to both outputs
                    writeMaternalPaternal(f_seq, s_seq);
                }
            }
        } else if (results.firstResult.kmerQC) {
            writeSingleReadBasedOnPhase(results.firstResult, f_seq);
            s_seq.write(*m_pDiscardWriter);
            ++m_readsDiscarded;
        } else if (results.secondResult.kmerQC) {
            writeSingleReadBasedOnPhase(results.secondResult, s_seq);
            f_seq.write(*m_pDiscardWriter);
            ++m_readsDiscarded;
        } else {
            f_seq.write(*m_pDiscardWriter);
            s_seq.write(*m_pDiscardWriter);
            m_readsDiscarded = m_readsDiscarded + 2;
        }
    } else {
        if((results.firstResult.kmerQC && results.secondResult.kmerQC) || m_pDiscardWriter == NULL)
        {
            f_seq.write(*m_pCorrectedWriter);
            s_seq.write(*m_pCorrectedWriter);
            m_readsKept = m_readsKept + 2;
        } else if (results.firstResult.kmerQC) {
            f_seq.write(*m_pCorrectedWriter);
            s_seq.write(*m_pDiscardWriter);
            ++m_readsKept;
            ++m_readsDiscarded;
        } else if (results.secondResult.kmerQC) {
            f_seq.write(*m_pDiscardWriter);
            s_seq.write(*m_pCorrectedWriter);
            ++m_readsKept;
            ++m_readsDiscarded;
        } else {
            f_seq.write(*m_pDiscardWriter);
            s_seq.write(*m_pDiscardWriter);
            m_readsDiscarded = m_readsDiscarded + 2;
        }
    }
}


void FilterParentPostProcess::writeSingleRead(const FilterParentResult& result, const SeqRecord& seq) {
    ++m_readsKept;
    bool neither_m_or_f = true;
    if (result.inMother == PARENT_ABOVE_THRESHOLD) {
        seq.write(*m_pMaternalWriter);
        ++m_reads_maternal;
        neither_m_or_f = false;
    }
    if (result.inFather == PARENT_ABOVE_THRESHOLD) {
        seq.write(*m_pPaternalWriter);
        ++m_reads_paternal;
        neither_m_or_f = false;
    }
    if (neither_m_or_f) {
        if ((result.inFather == PARENT_ABSENT) && (result.inMother == PARENT_ABSENT)) {
            ++m_passed_QC_but_neither_m_or_f;
            seq.write(*m_pMaternalWriter); seq.write(*m_pPaternalWriter); // Write to both streams -- could be de-novo mutation
            seq.write(*m_pNeitherParentWriter);
        } else {
            ++m_passed_QC_and_once_in_m_or_f;
        }
    }
}

void FilterParentPostProcess::writeBothReadsBasedOnOne(const FilterParentResult& result, const SeqRecord& f_seq, const SeqRecord& s_seq) {
    bool neither_m_or_f = true;
    if (result.inMother == PARENT_ABOVE_THRESHOLD) {
        f_seq.write(*m_pMaternalWriter); s_seq.write(*m_pMaternalWriter);
        m_reads_maternal = m_reads_maternal + 2;
        neither_m_or_f = false;
    }
    if (result.inFather == PARENT_ABOVE_THRESHOLD) {
        f_seq.write(*m_pPaternalWriter); s_seq.write(*m_pPaternalWriter);
        m_reads_paternal = m_reads_paternal + 2;
        neither_m_or_f = false;
    }
    if (neither_m_or_f) {
        if ((result.inFather == PARENT_ABSENT) && (result.inMother == PARENT_ABSENT)) {
            ++m_passed_QC_but_neither_m_or_f;
            f_seq.write(*m_pMaternalWriter); s_seq.write(*m_pMaternalWriter);
            f_seq.write(*m_pPaternalWriter); s_seq.write(*m_pPaternalWriter);
        } else {
            ++m_passed_QC_and_once_in_m_or_f;
        }
    }
}

void FilterParentPostProcess::collectMetrics(const std::string& originalSeq,
                                            const std::string& correctedSeq,
                                            const std::string& qualityStr)
{
    size_t precedingLen = 2;
    for(size_t i = 0; i < originalSeq.length(); ++i)
    {
        char qc = !qualityStr.empty() ? qualityStr[i] : '\0';
        char ob = originalSeq[i];
        
        ++m_totalBases;
        
        m_positionMetrics.incrementSample(i);
        
        if(!qualityStr.empty())
            m_qualityMetrics.incrementSample(qc);
        
        m_originalBaseMetrics.incrementSample(ob);
        
        std::string precedingMer;
        if(i > precedingLen)
        {
            precedingMer = originalSeq.substr(i - precedingLen, precedingLen);
            m_precedingSeqMetrics.incrementSample(precedingMer);
        }
        
        if(originalSeq[i] != correctedSeq[i])
        {
            m_positionMetrics.incrementError(i);
            if(!qualityStr.empty())
                m_qualityMetrics.incrementError(qc);
            m_originalBaseMetrics.incrementError(ob);
            
            if(!precedingMer.empty())
            {
                m_precedingSeqMetrics.incrementError(precedingMer);
            }
            ++m_totalErrors;
        }
    }
}
