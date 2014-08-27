//
//  FMMergeTrioProcess.h
//  sga_git
//
//  Created by Milan Malinsky on 18/04/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//

#ifndef __sga_git__FMMergeTrioProcess__
#define __sga_git__FMMergeTrioProcess__


#include "FMMergeProcess.h"

enum ContributingParent
{
    MATERNAL = 1,
    PATERNAL = 2
};

enum ParentPresence
{
    ABSENT_IN_PARENT,
    PRESENT_IN_PARENT
};

enum BothParentPresenceSummary
{
    PRESENT_IN_BOTH,
    ABSENT_IN_MOTHER,
    ABSENT_IN_FATHER,
    ABSENT_IN_BOTH
};

enum UniqueIfPhased
{
    KEEP_FIRST,
    KEEP_SECOND,
    INCONCLUSIVE
};

class VariantInfoVectors {
public:
    VariantInfoVectors() {}
    
    std::vector<OverlapBlock> overlapVec;
    std::vector<std::string> fullStringVec;
    std::vector<std::string> vertexStringVec;
    std::vector<std::string> vertexIDVec;
    std::vector<bool> isRCVec;
};

class ParentCheckSequenceResult
{
public:
    ParentCheckSequenceResult() : inFather(ABSENT_IN_PARENT),inMother(ABSENT_IN_PARENT),bothParentsSummary(ABSENT_IN_BOTH) {}
    
   // DNAString sequence;

    void setPresenceInFather(ParentPresence f) {
        inFather = f;
        updateBothParentsSummary();
    }
    void setPresenceInMother(ParentPresence m) {
        inMother = m;
        updateBothParentsSummary();
    }
    ParentPresence getPresenceInFather() {
        return (inFather);
    }
    ParentPresence getPresenceInMother() {
        return (inMother);
    }
    
    BothParentPresenceSummary bothParentsSummary;
    
private:
    ParentPresence inFather;
    ParentPresence inMother;
    void updateBothParentsSummary() {
        if (inFather == PRESENT_IN_PARENT && inMother == PRESENT_IN_PARENT) {
            bothParentsSummary = PRESENT_IN_BOTH;
        }
        if (inFather == ABSENT_IN_PARENT && inMother == PRESENT_IN_PARENT) {
            bothParentsSummary = ABSENT_IN_FATHER;
        }
        if (inFather == PRESENT_IN_PARENT && inMother == ABSENT_IN_PARENT) {
            bothParentsSummary = ABSENT_IN_MOTHER;
        }
        if (inFather == ABSENT_IN_PARENT && inMother == ABSENT_IN_PARENT) {
            bothParentsSummary = ABSENT_IN_BOTH;
        }
    }
};

// Compute the overlap blocks for reads
class FMMergeTrioProcess
{
public:
    FMMergeTrioProcess(const OverlapAlgorithm* pOverlapper,
                       int minOverlap, BitVector* pMarkedReads, BWTIndexSet& motherIndices,
                       BWTIndexSet& fatherIndices, ContributingParent contribParent);
    
    ~FMMergeTrioProcess();
    
    FMMergeResult process(const SequenceWorkItem& item);
    
private:
    
    // Add the edges of X as candidates to the graph
    void addCandidates(StringGraph* pGraph,
                       const Vertex* pX,
                       const Edge* pEdgeToX,
                       OverlapBlockList* pBlockList,
                       FMMergeQueue& candidateQueue);
    
    // Check whether the candidate can be merged into the current graph
    bool checkCandidate(Edge* pXY, const OverlapBlockList* pBlockList) const;
    bool checkCandidateAndTestHet(const FMMergeCandidate& candidate, OverlapBlockList* pBlockList);
    
    // Check is a branch in a graph could be caused by a heterozygous site
    UniqueIfPhased checkIfHet(VariantInfoVectors& vi, const OverlapBlockList* pBlockList, const Vertex* pX, const EdgeDir& direction);
    void checkPresenceInParents(const std::string& seq, ParentCheckSequenceResult& result, const int kmerLength = 31, const int startBase = 0);
    
    // Extend path from a particular vertex in only one direction
    StringGraph* buildVariantGraph(const std::string& rootId, const std::string& rootSequence, const bool rootIsRC, EdgeDir direction);
    void extendForwardOnly(StringGraph* pGraph, Vertex* pXstart, EdgeDir& rootOrientation);
    void extendBackwardOnly(StringGraph* pGraph, Vertex* pXstart, EdgeDir& rootOrientation);
    
    //
    std::string makeVertexID(BWTInterval interval);
    
    // Just some printing for debugging
    void printVariantPathSummary(const std::string& keepRootString, const std::string& variantRootString, std::string& mergedSeqKeep, std::string& mergedSeqVariant, const Vertex* pX, int num_v_keep, int num_v_variant,const EdgeDir& direction);
    
    const OverlapAlgorithm* m_pOverlapper;
    const int m_minOverlap;
    BitVector* m_pMarkedReads;
    BWTIndexSet m_motherIndices;
    BWTIndexSet m_fatherIndices;
    ContributingParent m_contribParent;
};

#endif /* defined(__sga_git__FMMergeTrioProcess__) */
