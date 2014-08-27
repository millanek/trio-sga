//
//  FMMergeTrioProcess.cpp
//  sga_git
//
//  Created by Milan Malinsky on 18/04/2014.
//  Copyright (c) 2014 University of Cambridge. All rights reserved.
//

#include "FMMergeTrioProcess.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"


//
FMMergeTrioProcess::FMMergeTrioProcess(const OverlapAlgorithm* pOverlapper, int minOverlap, BitVector* pMarkedReads, BWTIndexSet& motherIndices, BWTIndexSet& fatherIndices, ContributingParent contribParent) :
m_pOverlapper(pOverlapper),
m_minOverlap(minOverlap),
m_pMarkedReads(pMarkedReads),
m_motherIndices(motherIndices),
m_fatherIndices(fatherIndices),
m_contribParent(contribParent)
{
    
}

//
FMMergeTrioProcess::~FMMergeTrioProcess()
{
    
    
}

//
FMMergeResult FMMergeTrioProcess::process(const SequenceWorkItem& item)
{
    // Calculate the intervals in the forward FM-index for this read
    const BWT* pBWT = m_pOverlapper->getBWT();
    
    // Find the interval in the fm-index containing the read
    std::string readString = item.read.seq.toString();
    BWTInterval readInterval = BWTAlgorithms::findInterval(pBWT, readString);
    
    // Update the interval by looking for the $ to map the interval into indices in the bit array
    BWTAlgorithms::updateInterval(readInterval, '$', pBWT);
    
    // The read must be present in the index
    assert(readInterval.isValid());
    
    // Check if this read has been used yet
    bool used = false;
    for(int64_t i = readInterval.lower; i <= readInterval.upper; ++i)
    {
        if(m_pMarkedReads->test(i))
        {
            used = true;
            break;
        }
    }
    
    FMMergeResult result;
    if(!used)
    {
        // Construct a new local graph around this read
        StringGraph* pGraph = new StringGraph;
        std::stringstream ssID;
        ssID << "IDX-" << readInterval.lower;
        std::string rootID = ssID.str();
        
        Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(rootID, readString);
        pGraph->addVertex(pVertex);
        
        // Add the root vertex to the result structure
        result.usedIntervals.push_back(readInterval);
        
        // Enqueue the read for overlap detection in both directions
        FMMergeQueue queue;
        
        // Construct the overlap block list for this node
        SeqRecord record;
        record.id = pVertex->getID();
        record.seq = pVertex->getSeq().toString();
        OverlapBlockList blockList;
        m_pOverlapper->overlapRead(record, m_minOverlap, &blockList);
        
        removeContainmentBlocks((int)pVertex->getSeqLen(), &blockList);
        
        // Construct the initial list of candidates from the block list
        addCandidates(pGraph, pVertex, NULL, &blockList, queue);
        
       // std::cout << "Root vertex:\t\t" << pVertex->getSeq().toString() << std::endl;
       // std::cout << std::endl;
        
        while(!queue.empty())
        {
            FMMergeCandidate currCandidate = queue.front();
            queue.pop();
            
            // Determine whether this is a valid vertex to merge or not.
            // It is valid if it has a single edge in the direction of the vertex
            // that added it to the candidate list
            SeqRecord record;
            record.id = currCandidate.pVertex->getID();
            record.seq = currCandidate.pVertex->getSeq().toString();
            
            OverlapBlockList candidateBlockList;
            m_pOverlapper->overlapRead(record, m_minOverlap, &candidateBlockList);
            removeContainmentBlocks((int)currCandidate.pVertex->getSeqLen(), &candidateBlockList);
            
            bool validMergeNode = checkCandidateAndTestHet(currCandidate, &candidateBlockList);
            if(validMergeNode)
            {
                addCandidates(pGraph, currCandidate.pVertex, currCandidate.pEdge, &candidateBlockList, queue);
                result.usedIntervals.push_back(currCandidate.interval);
            }
            else
            {
                // Mark this vertex for later removal
                currCandidate.pVertex->setColor(GC_RED);
            }
        }
        
        // The graph has now been constructed. Remove all the nodes that are marked invalid for merging
        pGraph->sweepVertices(GC_RED);
        
        SGDuplicateVisitor dupVisit(true);
        pGraph->visit(dupVisit);
        
        // Merge nodes
        pGraph->simplify();
        
        // If there was a cycle in the graph, it is possible that more than 1 vertex
        // remains in the graph. Copy the vertex sequences into the result object.
        pGraph->getVertexSequences(result.mergedSequences);
        result.isMerged = true;
        delete pGraph;
    }
    else
    {
        result.isMerged = false;
    }
    
    if(result.isMerged)
    {
        // If some work was performed, update the bitvector so other threads do not try to merge the same set of reads.
        // This uses compare-and-swap instructions to ensure the uppdate is atomic.
        // If some other thread has merged this set (and updated
        // the bitvector), we discard all the merged data.
        
        // As a given set of reads should all be merged together, we only need to make sure we atomically update
        // the bit for the read with the lowest index in the set.
        
        // Sort the intervals into ascending order and remove any duplicate intervals (which can occur
        // if the subgraph has a simple cycle)
        std::sort(result.usedIntervals.begin(), result.usedIntervals.end(), BWTInterval::compare);
        std::vector<BWTInterval>::iterator newEnd = std::unique(result.usedIntervals.begin(),
                                                                result.usedIntervals.end(),
                                                                BWTInterval::equal);
        result.usedIntervals.erase(newEnd, result.usedIntervals.end());
        
        // Check if the bit in the vector has already been set for the lowest read index
        // If it has some other thread has already output this set so we do nothing
        int64_t lowestIndex = result.usedIntervals.front().lower;
        bool currentValue = m_pMarkedReads->test(lowestIndex);
        bool updateSuccess = false;
        
        if(currentValue == false)
        {
            // Attempt to update the bit vector with an atomic CAS. If this returns false
            // the bit was set by some other thread
            updateSuccess = m_pMarkedReads->updateCAS(lowestIndex, currentValue, true);
        }
        
        if(updateSuccess)
        {
            // We successfully atomically set the bit for the first read in this set
            // to true. We can safely update the rest of the bits and keep the merged sequences
            // for output.
            std::vector<BWTInterval>::const_iterator iter = result.usedIntervals.begin();
            for(; iter != result.usedIntervals.end(); ++iter)
            {
                for(int64_t i = iter->lower; i <= iter->upper; ++i)
                {
                    if(i == lowestIndex) //already set
                        continue;
                    
                    currentValue = m_pMarkedReads->test(i);
                    if(currentValue)
                    {
                        // This value should not be true, emit a warning
                        std::cout << "Warning: Bit " << i << " was set outside of critical section\n";
                    }
                    else
                    {
                        m_pMarkedReads->updateCAS(i, currentValue, true);
                    }
                }
            }
        }
        else
        {
            // Some other thread merged these reads already, discard the intermediate
            // data and set the result to false
            result.mergedSequences.clear();
            result.usedIntervals.clear();
            result.isMerged = false;
        }
    }
    return result;
}



// Add the edges starting from pX as candidate vertices
// using the blockList.
// Precondition: pX is a valid vertex in the merge graph. In other words, there is a
// unique assembly that includes pX and the root vertex.
void FMMergeTrioProcess::addCandidates(StringGraph* pGraph, const Vertex* pX, const Edge* pEdgeToX, OverlapBlockList* pBlockList, FMMergeQueue& candidateQueue)
{
    // Count the number of edges in each direction.
    size_t numAnti = 0;
    size_t numSense = 0;
    
    for(OverlapBlockList::const_iterator iter = pBlockList->begin(); iter != pBlockList->end(); ++iter)
    {
        if(iter->getEdgeDir() == ED_SENSE)
            ++numSense;
        if(iter->getEdgeDir() == ED_ANTISENSE)
            ++numAnti;
    }
    
    if (numAnti == 2) {
        VariantInfoVectors vi;
        UniqueIfPhased keep = checkIfHet(vi, pBlockList, pX, ED_ANTISENSE);
        if (keep == KEEP_FIRST) {
            pBlockList->remove(vi.overlapVec[1]);
            StringGraph* pKeepGraph = buildVariantGraph(vi.vertexIDVec[0], vi.vertexStringVec[0], vi.isRCVec[0],ED_ANTISENSE);
            StringGraph* pVariantGraph = buildVariantGraph(vi.vertexIDVec[1], vi.vertexStringVec[1], vi.isRCVec[1],ED_ANTISENSE);
            int num_v = (int)pVariantGraph->getNumVertices() - 1; int num_v_keep = (int)pKeepGraph->getNumVertices() - 1;
            pKeepGraph->simplify(); pVariantGraph->simplify(); // Merge nodes
            std::vector<std::string> mergedSequences; std::vector<std::string> mergedSequencesKeep;
            pVariantGraph->getVertexSequences(mergedSequences); pKeepGraph->getVertexSequences(mergedSequencesKeep);
            
            printVariantPathSummary(vi.fullStringVec[0],vi.fullStringVec[1], mergedSequencesKeep[0], mergedSequences[0], pX, num_v_keep, num_v,ED_ANTISENSE);
        }
        if (keep == KEEP_SECOND) {
            pBlockList->remove(vi.overlapVec[0]);
            StringGraph* pKeepGraph = buildVariantGraph(vi.vertexIDVec[1], vi.vertexStringVec[1], vi.isRCVec[1],ED_ANTISENSE);
            StringGraph* pVariantGraph = buildVariantGraph(vi.vertexIDVec[0], vi.vertexStringVec[0], vi.isRCVec[0],ED_ANTISENSE);
            int num_v = (int)pVariantGraph->getNumVertices() - 1; int num_v_keep = (int)pKeepGraph->getNumVertices() - 1;
            pKeepGraph->simplify(); pVariantGraph->simplify(); // Merge nodes
            std::vector<std::string> mergedSequences; std::vector<std::string> mergedSequencesKeep;
            pVariantGraph->getVertexSequences(mergedSequences); pKeepGraph->getVertexSequences(mergedSequencesKeep);
            
            printVariantPathSummary(vi.fullStringVec[1],vi.fullStringVec[0], mergedSequencesKeep[0], mergedSequences[0], pX, num_v_keep, num_v,ED_ANTISENSE);
        }
    }
    
    if (numSense == 2) {
        VariantInfoVectors vi;
        UniqueIfPhased keep = checkIfHet(vi, pBlockList, pX, ED_SENSE);
        
        if (keep == KEEP_FIRST) {
            pBlockList->remove(vi.overlapVec[1]);
            StringGraph* pKeepGraph = buildVariantGraph(vi.vertexIDVec[0], vi.vertexStringVec[0], vi.isRCVec[0],ED_SENSE);
            StringGraph* pVariantGraph = buildVariantGraph(vi.vertexIDVec[1], vi.vertexStringVec[1], vi.isRCVec[1],ED_SENSE);
            int num_v = (int)pVariantGraph->getNumVertices() - 1; int num_v_keep = (int)pKeepGraph->getNumVertices() - 1;
            pKeepGraph->simplify(); pVariantGraph->simplify(); // Merge nodes
            std::vector<std::string> mergedSequences; std::vector<std::string> mergedSequencesKeep;
            pVariantGraph->getVertexSequences(mergedSequences); pKeepGraph->getVertexSequences(mergedSequencesKeep);
            
            printVariantPathSummary(vi.fullStringVec[0],vi.fullStringVec[1], mergedSequencesKeep[0], mergedSequences[0], pX, num_v_keep, num_v,ED_SENSE);
        }
        if (keep == KEEP_SECOND) {
            pBlockList->remove(vi.overlapVec[0]);
            StringGraph* pKeepGraph = buildVariantGraph(vi.vertexIDVec[1], vi.vertexStringVec[1], vi.isRCVec[1],ED_SENSE);
            StringGraph* pVariantGraph = buildVariantGraph(vi.vertexIDVec[0], vi.vertexStringVec[0], vi.isRCVec[0],ED_SENSE);
            int num_v = (int)pVariantGraph->getNumVertices() - 1; int num_v_keep = (int)pKeepGraph->getNumVertices() - 1;
            pKeepGraph->simplify(); pVariantGraph->simplify(); // Merge nodes
            std::vector<std::string> mergedSequences; std::vector<std::string> mergedSequencesKeep;
            pVariantGraph->getVertexSequences(mergedSequences); pKeepGraph->getVertexSequences(mergedSequencesKeep);
            
            printVariantPathSummary(vi.fullStringVec[1],vi.fullStringVec[0], mergedSequencesKeep[0], mergedSequences[0], pX, num_v_keep, num_v,ED_SENSE);
        }
    }
    
    // For each edge block, if it is unique for the direction add the vertex it describes as a candidate
    for(OverlapBlockList::const_iterator iter = pBlockList->begin(); iter != pBlockList->end(); ++iter)
    {
        EdgeDir currDir = iter->getEdgeDir();
        if((currDir == ED_SENSE && numSense == 1) ||
           (currDir == ED_ANTISENSE && numAnti == 1))
        {
            // Skip edges in the direction to X
            if(pEdgeToX != NULL && pEdgeToX->getTwinDir() == currDir)
                continue;
            
            // Construct new candidate vertices and add them to the graph
            std::string vertexID = iter->toCanonicalID();
            assert(vertexID != pX->getID());
            std::string vertexSeq = iter->getFullString(pX->getSeq().toString());
            Overlap ovrXY = iter->toOverlap(pX->getID(), vertexID, (int)pX->getSeqLen(), (int)vertexSeq.length());
            
            // The vertex may already exist in the graph if the graph contains a loop
            Vertex* pY = pGraph->getVertex(vertexID);
            
            // Generate the new vertex
            if(pY == NULL)
            {
                pY = new(pGraph->getVertexAllocator()) Vertex(vertexID, vertexSeq);
                pGraph->addVertex(pY);
            }
            
            // Construct a description of the edge based on the overlap
            EdgeDesc ed = SGAlgorithms::overlapToEdgeDesc(pY, ovrXY);
            
            // If an edge with the same description as XY exists for X do not add a new edge or candidate
            assert(!pX->hasEdge(ed));
            
            // Construct the found edge and add it to the graph
            Edge* pXY = SGAlgorithms::createEdgesFromOverlap(pGraph, ovrXY, false);
            
            // Add the new vertex as a candidate
            FMMergeCandidate candidate;
            candidate.pVertex = pY;
            candidate.pEdge = pXY;
            candidate.interval = iter->getCanonicalInterval();
            candidateQueue.push(candidate);
            
        }
    }
}

// Build a unipath in one direction only
StringGraph* FMMergeTrioProcess::buildVariantGraph(const std::string& rootId, const std::string& rootSequence, const bool rootIsRC, EdgeDir buildDirection) {
    StringGraph* pVariantGraph = new StringGraph;
    Vertex* pVertex = new(pVariantGraph->getVertexAllocator()) Vertex(rootId, rootSequence);
    pVariantGraph->addVertex(pVertex);
    EdgeDir rootOrientation;
   // if (buildDirection == ED_SENSE) {
        if (rootIsRC) rootOrientation = ED_ANTISENSE; else rootOrientation = ED_SENSE;
   // } else if (buildDirection == ED_ANTISENSE) {
   //     if (rootIsRC) rootOrientation = ED_SENSE; else rootOrientation = ED_ANTISENSE;
   // }
    if (buildDirection == ED_SENSE) {
        extendForwardOnly(pVariantGraph, pVertex, rootOrientation);
    } else if (buildDirection == ED_ANTISENSE) {
        extendBackwardOnly(pVariantGraph, pVertex, rootOrientation);
    }
    return pVariantGraph;
}

void FMMergeTrioProcess::extendForwardOnly(StringGraph* pGraph, Vertex* pXstart, EdgeDir& rootOrientation) {
    int numSense; int numAntisense; int numBackwardBranches = 0; Vertex* pX = pXstart; Vertex* pY; EdgeDir o = rootOrientation;
    do {
        SeqRecord record;
        record.id = pX->getID();
        record.seq = pX->getSeq().toString();
        OverlapBlockList blockList;
        m_pOverlapper->overlapRead(record, m_minOverlap, &blockList);
        removeContainmentBlocks((int)pX->getSeqLen(), &blockList);
        OverlapBlock senseBlock; OverlapBlock antiSenseBlock;
        numSense = 0; numAntisense = 0;
        for(OverlapBlockList::const_iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
        {
            if(iter->getEdgeDir() == ED_SENSE && o == ED_SENSE) { ++numSense; senseBlock = *iter; }
            else if(iter->getEdgeDir() == ED_ANTISENSE && o == ED_ANTISENSE) { ++numSense; senseBlock = *iter; }
            else { antiSenseBlock = *iter; numAntisense++; }
        }
        
       // std::cout << "We have: " << numSense << " overlaps" << std::endl;
        if (numSense == 1) {
            std::string vertexID = senseBlock.toCanonicalID();
            assert(vertexID != pX->getID());
            std::string vertexSeq = senseBlock.getFullString(pX->getSeq().toString());
            Overlap ovrXY = senseBlock.toOverlap(pX->getID(), vertexID, (int)pX->getSeqLen(), (int)vertexSeq.length());
            
            if (ovrXY.match.isRC() && o == ED_SENSE) { o = ED_ANTISENSE; }
            else if (ovrXY.match.isRC() && o == ED_ANTISENSE) { o = ED_SENSE; }
            
            // The vertex may already exist in the graph if the graph contains a loop
            pY = pGraph->getVertex(vertexID);
            
            // Generate the new vertex
            if(pY == NULL)
            {
                pY = new(pGraph->getVertexAllocator()) Vertex(vertexID, vertexSeq);
                pGraph->addVertex(pY);
            }
            
            // Construct a description of the edge based on the overlap
            EdgeDesc ed = SGAlgorithms::overlapToEdgeDesc(pY, ovrXY);
            
            // If an edge with the same description as XY exists for X do not add a new edge or candidate
            assert(!pX->hasEdge(ed));
            
            // Construct the found edge and add it to the graph
            Edge* pXY = SGAlgorithms::createEdgesFromOverlap(pGraph, ovrXY, false);
            
            // Check the other direction
            SeqRecord record;
            record.id = pY->getID();
            record.seq = pY->getSeq().toString();
            OverlapBlockList candidateBlockList;
            m_pOverlapper->overlapRead(record, m_minOverlap, &candidateBlockList);
            removeContainmentBlocks((int)pY->getSeqLen(), &candidateBlockList);
            /*
            EdgeDir mergeDir = pXY->getTwinDir();
            size_t mergeDirCount = 0;
            for(OverlapBlockList::const_iterator iter = candidateBlockList.begin(); iter != candidateBlockList.end(); ++iter)
            {
                if(iter->getEdgeDir() == mergeDir)
                    mergeDirCount += 1;
            }
            assert(mergeDirCount > 0);
            std::cerr << mergeDirCount; */
            bool validMergeNode = checkCandidate(pXY, &candidateBlockList);
            if (!validMergeNode) {
                //std::cout << "Branch in backward direction..." << std::endl;
                numBackwardBranches++;
                pY->setColor(GC_RED);
            }
            pX = pY;
        } else {
            std::cout << "Next we have: " << numSense << " overlaps" << std::endl;
            for(OverlapBlockList::const_iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
            {
                if(iter->getEdgeDir() == ED_SENSE && o == ED_SENSE) {
                    std::cout << iter->getFullString(pX->getSeq().toString()) << std::endl; }
                if(iter->getEdgeDir() == ED_ANTISENSE && o == ED_ANTISENSE) {
                    std::cout << iter->getFullString(pX->getSeq().toString()) << std::endl;
                }
            }
            std::cout << "Total branches in the other direction: " << numBackwardBranches << std::endl;
        }
    } while (numSense == 1);
}

void FMMergeTrioProcess::extendBackwardOnly(StringGraph* pGraph, Vertex* pXstart, EdgeDir& rootOrientation) {
    int numAntisense; int numBackwardBranches = 0; Vertex* pX = pXstart; Vertex* pY; EdgeDir o = rootOrientation;
    do {
        SeqRecord record;
        record.id = pX->getID();
        record.seq = pX->getSeq().toString();
        OverlapBlockList blockList;
        m_pOverlapper->overlapRead(record, m_minOverlap, &blockList);
        removeContainmentBlocks((int)pX->getSeqLen(), &blockList);
        OverlapBlock senseBlock; OverlapBlock antiSenseBlock;
        numAntisense = 0;
        for(OverlapBlockList::const_iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
        {
            if(iter->getEdgeDir() == ED_ANTISENSE && o == ED_SENSE) { ++numAntisense; senseBlock = *iter; }
            else if(iter->getEdgeDir() == ED_SENSE && o == ED_ANTISENSE) { ++numAntisense; senseBlock = *iter; }
            // else { antiSenseBlock = *iter; }
        }
        
        // std::cout << "We have: " << numSense << " overlaps" << std::endl;
        if (numAntisense == 1) {
            std::string vertexID = senseBlock.toCanonicalID();
            assert(vertexID != pX->getID());
            std::string vertexSeq = senseBlock.getFullString(pX->getSeq().toString());
            Overlap ovrXY = senseBlock.toOverlap(pX->getID(), vertexID, (int)pX->getSeqLen(), (int)vertexSeq.length());
            if (ovrXY.match.isRC() && o == ED_SENSE) { o = ED_ANTISENSE; }
            else if (ovrXY.match.isRC() && o == ED_ANTISENSE) { o = ED_SENSE; }
            
            // The vertex may already exist in the graph if the graph contains a loop
            pY = pGraph->getVertex(vertexID);
            
            // Generate the new vertex
            if(pY == NULL)
            {
                pY = new(pGraph->getVertexAllocator()) Vertex(vertexID, vertexSeq);
                pGraph->addVertex(pY);
            }
            
            // Construct a description of the edge based on the overlap
            EdgeDesc ed = SGAlgorithms::overlapToEdgeDesc(pY, ovrXY);
            
            // If an edge with the same description as XY exists for X do not add a new edge or candidate
            assert(!pX->hasEdge(ed));
            
            // Construct the found edge and add it to the graph
            Edge* pXY = SGAlgorithms::createEdgesFromOverlap(pGraph, ovrXY, false);
            
            // Check the other direction
            SeqRecord record;
            record.id = pY->getID();
            record.seq = pY->getSeq().toString();
            OverlapBlockList candidateBlockList;
            m_pOverlapper->overlapRead(record, m_minOverlap, &candidateBlockList);
            removeContainmentBlocks((int)pY->getSeqLen(), &candidateBlockList);
            /*
             EdgeDir mergeDir = pXY->getTwinDir();
             size_t mergeDirCount = 0;
             for(OverlapBlockList::const_iterator iter = candidateBlockList.begin(); iter != candidateBlockList.end(); ++iter)
             {
             if(iter->getEdgeDir() == mergeDir)
             mergeDirCount += 1;
             }
             assert(mergeDirCount > 0);
             std::cerr << mergeDirCount; */
            bool validMergeNode = checkCandidate(pXY, &candidateBlockList);
            if (!validMergeNode) {
                //std::cout << "Branch in backward direction..." << std::endl;
                numBackwardBranches++;
                pY->setColor(GC_RED);
            }
            pX = pY;
        } else {
            std::cout << "Next we have: " << numAntisense << " overlaps" << std::endl;
            for(OverlapBlockList::const_iterator iter = blockList.begin(); iter != blockList.end(); ++iter)
            {
                if(iter->getEdgeDir() == ED_ANTISENSE && o == ED_SENSE) {
                    std::cout << iter->getFullString(pX->getSeq().toString()) << std::endl; }
                if(iter->getEdgeDir() == ED_SENSE && o == ED_ANTISENSE) {
                    std::cout << iter->getFullString(pX->getSeq().toString()) << std::endl;
                }
            }
            std::cout << "Total branches in the other direction: " << numBackwardBranches << std::endl;
        }
    } while (numAntisense == 1);
}


// Check if the candidate node can be merged with the node it is linked to. Returns true if so
bool FMMergeTrioProcess::checkCandidate(Edge* pXY, const OverlapBlockList* pBlockList) const
{
    // Get the direction of the edge back to the node that generated this candidate
    // pEdge is the edge TO the candidate vertex so the twin direction is the direction away
    // from the candidate vertex
    EdgeDir mergeDir = pXY->getTwinDir();
    
    size_t mergeDirCount = 0;
    for(OverlapBlockList::const_iterator iter = pBlockList->begin(); iter != pBlockList->end(); ++iter)
    {
        if(iter->getEdgeDir() == mergeDir)
            mergeDirCount += 1;
    }
    
    assert(mergeDirCount > 0);

    //std::cerr << mergeDirCount;
    return mergeDirCount == 1;
}

// Check if the candidate node can be merged with the node it is linked to. Returns true if so
bool FMMergeTrioProcess::checkCandidateAndTestHet(const FMMergeCandidate& candidate, OverlapBlockList* pBlockList)
{
    // Get the direction of the edge back to the node that generated this candidate
    // pEdge is the edge TO the candidate vertex so the twin direction is the direction away
    // from the candidate vertex
    EdgeDir mergeDir = candidate.pEdge->getTwinDir();
    
    size_t mergeDirCount = 0;
    for(OverlapBlockList::const_iterator iter = pBlockList->begin(); iter != pBlockList->end(); ++iter)
    {
        if(iter->getEdgeDir() == mergeDir)
            mergeDirCount += 1;
    }
    
    if (mergeDirCount == 2) {
        VariantInfoVectors vi;
        UniqueIfPhased keep = checkIfHet(vi, pBlockList, candidate.pVertex, mergeDir);
        if (keep == KEEP_FIRST) {
            pBlockList->remove(vi.overlapVec[1]);
            mergeDirCount--;
        }
        if (keep == KEEP_SECOND) {
            pBlockList->remove(vi.overlapVec[0]);
            mergeDirCount--;
        }
    }
    
    assert(mergeDirCount > 0);
    
    //std::cerr << mergeDirCount;
    return mergeDirCount == 1;
}



UniqueIfPhased FMMergeTrioProcess::checkIfHet(VariantInfoVectors& vi, const OverlapBlockList* pBlockList, const Vertex* pX, const EdgeDir& direction) {
    ParentCheckSequenceResult first; ParentCheckSequenceResult second; int i = 1;
    for(OverlapBlockList::const_iterator iter = pBlockList->begin(); iter != pBlockList->end(); ++iter) {
        if (iter->getEdgeDir() == direction) {
            vi.overlapVec.push_back(*iter);
            std::string vertexID = iter->toCanonicalID(); vi.vertexIDVec.push_back(vertexID);
            assert(vertexID != pX->getID());
            std::string vertexSeq = iter->getFullString(pX->getSeq().toString()); vi.vertexStringVec.push_back(vertexSeq);
            Overlap ovrXY = iter->toOverlap(pX->getID(), vertexID, (int)pX->getSeqLen(), (int)vertexSeq.length());
            std::string fullStr;
            if (ovrXY.match.isRC()) {
                fullStr = reverseComplement(vertexSeq); vi.fullStringVec.push_back(fullStr); vi.isRCVec.push_back(true);
            } else {
                fullStr = vertexSeq; vi.fullStringVec.push_back(fullStr); vi.isRCVec.push_back(false);
            }
            if (direction == ED_SENSE) {
                if (i == 1) { checkPresenceInParents(fullStr, first, 31,m_minOverlap-1); i++; }
                else if (i == 2) { checkPresenceInParents(fullStr, second, 31,m_minOverlap-1); i++; }
            } else if (direction == ED_ANTISENSE) {
                if (i == 1) { checkPresenceInParents(reverseComplement(fullStr), first, 31,m_minOverlap-1); i++; }
                else if (i == 2) { checkPresenceInParents(reverseComplement(fullStr), second, 31,m_minOverlap-1); i++; }
            }
        }
    } assert(i <= 3);
    
    if (m_contribParent == MATERNAL) {
        if (first.bothParentsSummary == ABSENT_IN_MOTHER && second.bothParentsSummary == PRESENT_IN_BOTH) {
            return KEEP_SECOND;
        } else if (first.bothParentsSummary == PRESENT_IN_BOTH && second.bothParentsSummary == ABSENT_IN_MOTHER) {
            return KEEP_FIRST;
        }
    }
    if (m_contribParent == PATERNAL) {
        if (first.bothParentsSummary == ABSENT_IN_FATHER && second.bothParentsSummary == PRESENT_IN_BOTH) {
            return KEEP_SECOND;
        } else if (first.bothParentsSummary == PRESENT_IN_BOTH && second.bothParentsSummary == ABSENT_IN_FATHER) {
            return KEEP_FIRST;
        }
    }
    return INCONCLUSIVE;
}

void FMMergeTrioProcess::checkPresenceInParents(const std::string& seq, ParentCheckSequenceResult& result, const int kmerLength, const int startBase) {
    int n = (int)seq.size();
    int nk = n - kmerLength + 1;
    
    std::vector<int> solidVectorMother(n, 0);
    std::vector<int> solidVectorFather(n, 0);
    bool allSolidMother = true;
    bool allSolidFather = true;
    
    for(int i = startBase; i < nk; ++i)
    {
        std::string kmer = seq.substr(i, kmerLength);
        size_t countMother = BWTAlgorithms::countSequenceOccurrences(kmer, m_motherIndices);
        size_t countFather = BWTAlgorithms::countSequenceOccurrences(kmer, m_fatherIndices);
        
        if (countMother == 0) { allSolidMother = false; }
        if (countFather == 0) { allSolidFather = false; }
        // std::cout << "Count mother: " << countMother << " Count father: " << countFather << "\n";
    }
    if (allSolidMother) { result.setPresenceInMother(PRESENT_IN_PARENT); }
    else { result.setPresenceInMother(ABSENT_IN_PARENT); }
    if (allSolidFather) { result.setPresenceInFather(PRESENT_IN_PARENT); }
    else { result.setPresenceInFather(ABSENT_IN_PARENT); }
    
}

// Just some printing for dubugging
void FMMergeTrioProcess::printVariantPathSummary(const std::string& keepRootString, const std::string& variantRootString, std::string& mergedSeqKeep, std::string& mergedSeqVariant, const Vertex* pX, int num_v_keep, int num_v_variant, const EdgeDir& direction) {
    if (direction == ED_SENSE) std::cout << "Sense direction: " << std::endl;
    else if (direction == ED_ANTISENSE) std::cout << "Antisense direction: " << std::endl;
    std::cout << "Added by:\t\t" << pX->getSeq().toString() << std::endl;
    std::cout << "Keeping:\t\t" << keepRootString << std::endl;
    if (direction == ED_SENSE) {
        if (keepRootString != mergedSeqKeep.substr(0,100)) mergedSeqKeep = reverseComplement(mergedSeqKeep);
    } else {
        if (keepRootString != mergedSeqKeep.substr(mergedSeqKeep.length()-100)) mergedSeqKeep = reverseComplement(mergedSeqKeep);
    }
    
    std::cout << "Extended into:\t\t" << mergedSeqKeep << " num extensions:" << num_v_keep << std::endl;
    std::cout << "Not keeping:\t\t" << variantRootString << std::endl;
    if (direction == ED_SENSE) {
        if (variantRootString != mergedSeqVariant.substr(0,100)) mergedSeqVariant = reverseComplement(mergedSeqVariant);
    } else {
        if (variantRootString != mergedSeqVariant.substr(mergedSeqVariant.length()-100)) mergedSeqVariant = reverseComplement(mergedSeqVariant);
    }
    std::cout << "Extended into:\t\t" << mergedSeqVariant << " num extensions:" << num_v_variant << std::endl;
    std::cout << std::endl;
}

