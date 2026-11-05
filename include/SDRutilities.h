#ifndef SDR_UTILITIES_H
#define SDR_UTILITIES_H

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <utility>
#include <vector>
#include <set>

#include <array>
#include <string>
#include "SDRtypes.h"


namespace karyoke
{

// -------------------------------------------------- //
// --------- SDR MUTATION DETECTION HELPERS --------- //
// -------------------------------------------------- //

// oldStrandID always refers to one of the  original intact chromosomes (46 for humans), never to a
// newly-created strand ID produced by an earlier rearrangement. This is what allows every detect<Mutation>() 
// function below to group a cell's fragments by oldStrandID and treat that as "all the rearrangement activity 
// affecting a given original chromosome", rather than needing to trace strand lineage across events.


// Determines which original strand a record "belongs to" for karyogram slot placement. A record with exactly one centromere-bearing
// fragment unambiguously belongs to that fragment's strand, this is the correct rule for any BALANCED rearrangement, since by
// definition each derivative keeps exactly one parent's centromere.
inline int determineHomeStrandID(const SDRdataRecord& record)
{
    int centromereStrandID = -1;
    int centromereCount = 0;

    for (const SDRfragment& fragment : record.fragments)			// Find which fragment holds the centromere
    {
        if (fragment.hasCentromere)						// If the fragment has centromere
        {
            centromereStrandID = fragment.oldStrandID;				// Store the associated original strand ID of that fragment
            ++centromereCount;							// Increment the centromere count for that data record
        }
    }

    if (centromereCount == 1)							// Only supports single-centromere records for determining home strand
    {
        return centromereStrandID;
    }

    // Zero or 2+ centromeres, fall back to the original strand-ordering convention (first fragment old strand ID appearing in the data record is the home slot).
    return record.fragments.empty() ? -1 : record.fragments[0].oldStrandID;
}



















// Mutated strands are only those with new strand IDs > numOriginalStrands.
// EX: for the 46 human chromosomes with strand IDs (1-46), mutated strands begin with new Strand ID > 46. This will be consistent across SDR files.
inline bool isRearrangementCandidate(int newStrandID, int numOriginalStrands)
{
    return newStrandID > numOriginalStrands;
}


// Signed distance from position a to position b in bases (in Mbp, 1 bp = 0.000001 Mbp), rounded to a whole base. Positions are base
// coordinates: 50.000000 -> 50.000001 is the next base (touching, nothing missing), 50.000000 -> 50.000002 is a 1 bp gap. Rounding is what
// makes size boundaries reliable: 50.000100 - 50.000000 is not exactly 0.000100 in floating point, so plain Mbp comparisons can land on
// either side of it.
inline long long basesBetween(double a, double b)
{
    return std::llround((b - a) * 1.0e6);
}


// Threshold check used by every detector: a and b (positions in Mbp) are "the same" when they are at most `tolerance` apart, compared in whole
// bases (see basesBetween()) so a distance sitting exactly on the tolerance always counts as within it. Every structural-variation detector
// (deletion-type mutations, balanced inversions and translocations) uses svToleranceMbp below.
inline bool approxEqual(double a, double b, double tolerance)
{
    return std::llabs(basesBetween(a, b)) <= std::llround(tolerance * 1.0e6);
}


// Structural variations are >= 100 bp; anything smaller is an indel (Youk et al. 2024, chosen by Samy 2026-11-05). Two positions up to 100 bases
// apart (at most 99 bases missing between them, e.g. 50.000000 -> 50.000100) are therefore treated as touching by every structural-variation detector
// and the drawing; 50.000000 -> 50.000101 (100 bases missing) is a structural-variation gap. A balanced inversion or translocation may likewise lose
// up to 99 bases at its breakpoints and still count as balanced (Samy, 2026-11-05; it was 10 bp).
constexpr double svToleranceMbp = 0.000100;




// --------------------------------------------------------------------------- //
// INDELS: small insertions and deletions of < 100 bp
// --------------------------------------------------------------------------- //

// Indels are < 100 bp (at most 99 bases); 100 bp and larger is a structural variation, scored by the other detectors.
constexpr long long indelMaxBases = 99;

// Number of bases in a piece, counting both end bases (a piece with start == end is one base).
inline long long fragmentBases(const SDRfragment& fragment)
{
    return std::llabs(basesBetween(fragment.oldStartPosition, fragment.oldEndPosition)) + 1;
}

// A piece of < 100 bp, e.g. one inserted base (start == end).
inline bool isIndelSizedFragment(const SDRfragment& fragment)
{
    return fragmentBases(fragment) <= indelMaxBases;
}

// Returns a copy of the cell with every indel taken out of its rearranged linear records, and lists the indels removed when indels is non-null.
// Every structural-variation detector and the karyogram drawing work on this cleaned copy, so indels never reach them:
//  - Small insertion: a run of < 100 bp pieces (usually one-base pieces copied from random places in the genome) between two larger pieces, adding up to
//    < 100 bases, is dropped. Without this, the copied bases tied the receiving strand to every strand they were copied from (Miscellaneous on all of them,
//    or chromoplexy when two of those strands also had a translocation) and added to the chromothripsis fragment count.
//  - Small deletion: two consecutive pieces of the same strand and orientation, with 1-50 bases missing between them (or none when bases were inserted there),
//    e.g. 1/0/50.000000/0, 1/50.000005/263.2/1 is 4 bases deleted (50.000001-50.000004); 1/0/50.000000/0, 1/50.000100/263.2/1 is 99 (the largest indel),
//    are joined back into one piece, so the record reads as if the indel never happened (often a plain whole copy of the chromosome).
// Gaps of 100 bp or more, pieces of 100 bp or more, and pieces at either end of a record are left alone. Circular records and the baseline records (newStrandID <=
// numOriginalStrands) are not changed. A site with both a deletion and an insertion is ONE indel.
//  - Excised partner: as for long deletions, the deleted bases may be listed as their own record, or not at all. A single-piece record of the same strand
//    covering the deleted bases (each end within one base of the gap's, so 1/50.000001/50.000004/0 and 1/50.000000/50.000005/0 both match the gap above) is
//    part of the indel and is taken out too, so it is neither scored nor drawn; each such record is used for one indel only. Without this the leftover piece
//    was a change no other type explained, so it counted as a complex genomic rearrangement. The partner may be linear or circular and listed in either
//    orientation (Samy, 2026-11-03): an ecDNA of < 100 bp cut from an indel gap is scored as the indel deletion, not as ecDNA, and is not drawn, and so is a
//    reversed linear piece (unlike a long deletion, whose reversed excised piece counts as foreign material).
inline SDRsubHeader removeIndels(const SDRsubHeader& subHeader, int numOriginalStrands, std::vector<SDRindelEvent>* indels = nullptr)
{
    SDRsubHeader cleaned = subHeader;
    std::set<int> excisedRecordIDs;							// Excised indel partners, removed once every record is processed

    // The record listing the deleted bases [gapLow, gapHigh] of strandID on their own (see above), or -1.
    auto findExcisedPartner = [&](int strandID, double gapLow, double gapHigh) -> int
    {
        for (const SDRdataRecord& candidate : subHeader.dataRecords)
        {
            if (!isRearrangementCandidate(candidate.newStrandID, numOriginalStrands) || candidate.fragments.size() != 1 || excisedRecordIDs.count(candidate.newStrandID))
            {
                continue;
            }

            const SDRfragment& piece = candidate.fragments[0];
            const double low = std::min(piece.oldStartPosition, piece.oldEndPosition);
            const double high = std::max(piece.oldStartPosition, piece.oldEndPosition);

            if (piece.oldStrandID == strandID
                && std::llabs(basesBetween(gapLow, low)) <= 1 && std::llabs(basesBetween(high, gapHigh)) <= 1)
            {
                return candidate.newStrandID;
            }
        }

        return -1;
    };

    for (SDRdataRecord& record : cleaned.dataRecords)
    {
        if (!record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands) || record.fragments.size() < 2)
        {
            continue;
        }

        const std::vector<SDRfragment>& original = record.fragments;
        std::vector<SDRfragment> kept = {original[0]};
        std::size_t i = 1;

        while (i < original.size())
        {
            // Run of small pieces starting at i, and the bases they add (a piece with start == end is one base).
            std::size_t runEnd = i;
            long long insertedBases = 0;

            while (runEnd < original.size() && isIndelSizedFragment(original[runEnd]))
            {
                insertedBases += fragmentBases(original[runEnd]);
                ++runEnd;
            }

            const bool insertion = runEnd > i && runEnd < original.size() && !isIndelSizedFragment(kept.back()) && insertedBases <= indelMaxBases;

            if (runEnd > i && !insertion)						// Small pieces that are not a flanked indel insertion stay as they are.
            {
                kept.insert(kept.end(), original.begin() + static_cast<long>(i), original.begin() + static_cast<long>(runEnd));
                i = runEnd;
                continue;
            }

            if (!insertion)
            {
                insertedBases = 0;
            }

            SDRfragment& left = kept.back();
            const SDRfragment& right = original[runEnd];
            const bool leftReversed = left.oldStartPosition > left.oldEndPosition;
            const bool rightReversed = right.oldStartPosition > right.oldEndPosition;

            // Bases missing between the two pieces in reading direction; only meaningful for the same strand and orientation. The next base
            // (distance 1) is touching, and so is a shared position (distance 0, how the other mutation types list a clean break).
            const long long distance = leftReversed ? basesBetween(right.oldStartPosition, left.oldEndPosition) : basesBetween(left.oldEndPosition, right.oldStartPosition);
            const long long deletedBases = std::max(0LL, distance - 1);
            const bool joinable = left.oldStrandID == right.oldStrandID && leftReversed == rightReversed && !isIndelSizedFragment(left) && !isIndelSizedFragment(right)
                                  && distance >= 0 && deletedBases <= indelMaxBases && (deletedBases > 0 || insertion);

            if (joinable || insertion)
            {
                int excisedStrandID = -1;

                if (joinable && deletedBases > 0)
                {
                    excisedStrandID = findExcisedPartner(left.oldStrandID, std::min(left.oldEndPosition, right.oldStartPosition), std::max(left.oldEndPosition, right.oldStartPosition));

                    if (excisedStrandID != -1)
                    {
                        excisedRecordIDs.insert(excisedStrandID);
                    }
                }

                if (indels)
                {
                    indels->push_back({left.oldStrandID, record.newStrandID, left.oldEndPosition, joinable ? static_cast<int>(deletedBases) : 0,
                                       static_cast<int>(insertedBases), excisedStrandID});
                }
            }

            if (joinable)
            {
                left.oldEndPosition = right.oldEndPosition;
                left.hasCentromere = left.hasCentromere || right.hasCentromere;
            }
            else
            {
                kept.push_back(right);
            }

            i = runEnd + 1;
        }

        record.fragments = kept;
    }

    cleaned.dataRecords.erase(std::remove_if(cleaned.dataRecords.begin(), cleaned.dataRecords.end(),
                                             [&](const SDRdataRecord& record) { return excisedRecordIDs.count(record.newStrandID) > 0; }),
                              cleaned.dataRecords.end());

    return cleaned;
}

// --------------------------------------------------------------------------- //
// Detects INDELS (< 100 bp insertions/deletions) within a single cell.
// --------------------------------------------------------------------------- //

// One event per site (see removeIndels() for the shapes). Indels are never part of a composite mutation, so there is no hierarchy: they are counted
// on the receiving strand only, never on the strands the inserted bases were copied from, and are not drawn on the karyogram.
inline std::vector<SDRindelEvent> detectIndels(const SDRsubHeader& subHeader, int numOriginalStrands)
{
    std::vector<SDRindelEvent> indels;
    removeIndels(subHeader, numOriginalStrands, &indels);
    return indels;
}













// Finds record pairs representing a dicentric+acentric outcome, now generalized to records with any number of fragments (not just
// exactly 2), covering translocations exchanging segments from the middle of a chromosome, not just its ends.
//
// Within the dicentric candidate record, exactly two fragments must carry a centromere, from two DIFFERENT original strands, every
// other fragment in that record must be acentric and reference one of those same two strands. The first listed centromere-bearing
// fragment's strand becomes the dicentric record's own home strand; the second becomes the paired acentric record's home strand. This
// ordering is deliberate, whichever strand's fragment is listed first is the one that "owns" the dicentric drawing.
inline std::map<int, int> findDicentricAcentricHomeOverrides(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader)
{
    std::map<int, int> overrides;
    std::set<int> claimedAcentricRecords;				// newStrandIDs of acentric records already paired with a dicentric record, prevents the same acentric record being matched to more than one dicentric candidate.
    std::set<int> claimedExcisedRecords;				// newStrandIDs of single-fragment deletion records already used to fill a gap for some earlier pair, prevents the same deletion being reused by a different pair.

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;
    const double tolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

    // Given pooled fragments for ONE strand (from a dicentric+acentric pair), checks whether they (plus zero or more ADDITIONAL 
    // single-fragment excised records for that strand) exactly span its full declared length. This is what lets a deletion coexist
    // with the translocation that created the dicentric/acentric pair. Only single-fragment excised pieces are matched (one deletion
    // per gap).
    auto tryTileWithDeletions = [&](std::vector<SDRfragment> pooled, int strandID, std::set<std::size_t> excludeIndices, std::vector<std::size_t>& usedExcisedIndices, bool allowUnmatchedGaps) -> bool
    {
        const std::size_t sizeIndex = static_cast<std::size_t>(strandID);
        if (sizeIndex >= masterHeader.intactChromosomeSizes.size())			// strandID must index a valid ORIGINAL chromosome size. Out of range means its full length is unknown, so quit and return false rather than guess.
        {
            return false;
        }

        const double fullSize = masterHeader.intactChromosomeSizes[sizeIndex];		// Size of chromosome in Mbp
        if (fullSize <= 0.0)								// A zero/negative declared size can't be tiled against.
        {
            return false;
        }

	// Sort THIS STRAND's own pooled fragments (gathered from both the dicentric and acentric candidates, by the
	// caller) by position, so the loop below can walk them in order and find the gaps between them.
        std::sort(pooled.begin(), pooled.end(), [](const SDRfragment& a, const SDRfragment& b)
        {
            return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
        });


	// Find deletions that may affect clean tiling of fragments to the original chromosome's intact sizes
        struct Gap
        {
            double startPos;
            double endPos;
        };

        std::vector<Gap> gaps;


        const double firstStart = std::min(pooled.front().oldStartPosition, pooled.front().oldEndPosition);

        if (firstStart > 0.0 && !approxEqual(firstStart, 0.0, tolerance))		// Leading gap, the pooled fragments don't reach back to position 0.
        {
            gaps.push_back({0.0, firstStart});
        }

        for (std::size_t k = 0; k + 1 < pooled.size(); ++k)				// Walk consecutive pooled fragments looking for a gap between each pair.
        {
            const double higherPos = std::max(pooled[k].oldStartPosition, pooled[k].oldEndPosition);
            const double lowerPos = std::min(pooled[k + 1].oldStartPosition, pooled[k + 1].oldEndPosition);

            if (approxEqual(higherPos, lowerPos, tolerance))				// Fragments are contiguous within tolerance, no gap here.
            {
                continue;
            }

            if (lowerPos <= higherPos)							// Fragments overlap or are out of order, not a valid tiling, reject the whole strand.
            {
                return false;
            }

            gaps.push_back({higherPos, lowerPos});					// Genuine gap between these two fragments.
        }

        const double lastEnd = std::max(pooled.back().oldStartPosition, pooled.back().oldEndPosition);

        if (fullSize > lastEnd && !approxEqual(lastEnd, fullSize, tolerance))		// Trailing gap, the pooled fragments don't reach the chromosome's true end.
        {
            gaps.push_back({lastEnd, fullSize});
        }

        if (gaps.empty())								// Clean tiling already, no deletion involved.
        {
            return true;
        }


        std::vector<bool> gapMatched(gaps.size(), false);				// Tracks which gaps have been filled by a matching deletion record so far.

        for (std::size_t r = 0; r < records.size(); ++r)				// Search every other record in the cell for a single-fragment deletion filling one of the gaps above.
        {
            if (excludeIndices.count(r))						// Skip records already reserved: the dicentric/acentric records themselves, and any excised record already claimed elsewhere.
            {
                continue;
            }

            const SDRdataRecord& candidate = records[r];

            if (!candidate.linear || !isRearrangementCandidate(candidate.newStrandID, numOriginalStrands))	// A deletion's excised piece is its own plain LINEAR record (not a circular ecDNA fragment), and must be a genuine rearrangement.
            {
                continue;
            }

            if (candidate.fragments.size() != 1)					// Only a single-fragment record can fill exactly one gap, a multi-fragment excised piece isn't matched here.
            {
                continue;
            }

            const SDRfragment& fragment = candidate.fragments[0];

            if (fragment.oldStrandID != strandID || fragment.oldStartPosition > fragment.oldEndPosition)	// Must reference this same strand, and not be inverted (a plain deletion has no orientation).
            {
                continue;
            }

            for (std::size_t g = 0; g < gaps.size(); ++g)
            {
                if (gapMatched[g])							// Already filled, don't double-claim the same gap.
                {
                    continue;
                }

                if (approxEqual(fragment.oldStartPosition, gaps[g].startPos, tolerance) &&
                    approxEqual(fragment.oldEndPosition, gaps[g].endPos, tolerance))
                {
                    gapMatched[g] = true;
                    usedExcisedIndices.push_back(r);					// Remember which record filled this gap, so the caller can mark it claimed.
                    break;
                }
            }
        }

        for (bool matched : gapMatched)							// Every gap must be filled, otherwise the tiling of records to the original chromosomes invovled is incomplete.
        {
            if (!matched && !allowUnmatchedGaps)					// The excised record is optional: in the relaxed second pass a gap with no excised record is accepted as a lone deletion.
            {
                return false;
            }
        }

        return true;
    };



    // Two passes: first every dicentric pairs only when ALL gaps in the pooled strands are explained by excised records (the strongest proof), then the remaining dicentrics
    // pair with gaps allowed to lack an excised record, since the excised record is optional in the SDR file.
    for (int pass = 0; pass < 2; ++pass)
    {
    const bool allowUnmatchedGaps = (pass == 1);

    for (std::size_t i = 0; i < records.size(); ++i)					// Outer loop: try every record as a potential DICENTRIC candidate.
    {
        const SDRdataRecord& dicentricCandidate = records[i];

        if (overrides.count(dicentricCandidate.newStrandID))				// Already paired in the strict pass.
        {
            continue;
        }

        if (!dicentricCandidate.linear || !isRearrangementCandidate(dicentricCandidate.newStrandID, numOriginalStrands))	// Must be a linear rearranged record.
        {
            continue;
        }

        std::vector<int> centromereStrandsInOrder;					// Strand ID of every centromere-bearing fragment, in the order they appear on the record.
        std::set<int> allStrandsInRecord;						// Every distinct strand referenced anywhere in the record.

        for (const SDRfragment& fragment : dicentricCandidate.fragments)		// Store strand IDs of all centromere bearing strands.
        {
            allStrandsInRecord.insert(fragment.oldStrandID);

            if (fragment.hasCentromere)
            {
                centromereStrandsInOrder.push_back(fragment.oldStrandID);
            }
        }

        if (centromereStrandsInOrder.size() != 2)					// We are looking for dicentrics here, the record must have exactly 2 centromeres.
        {
            continue;
        }

        const int strandX = centromereStrandsInOrder[0];				// First listed centromere-bearing strand becomes home for the DICENTRIC record.
        const int strandY = centromereStrandsInOrder[1];				// And the second listed centromere-bearing strand becomes home for the paired ACENTRIC record, once a valid match is confirmed below.

        if (strandX == strandY || allStrandsInRecord.size() != 2)			// The two centromeres must come from two DIFFERENT strands, and no fragment here may reference a third strand.
        {
            continue;
        }

        for (std::size_t j = 0; j < records.size(); ++j)				// Inner loop: search for the ACENTRIC record that pairs with this dicentric candidate.
        {
            if (i == j)
            {
                continue;
            }

            const SDRdataRecord& acentricCandidate = records[j];

            if (claimedAcentricRecords.count(acentricCandidate.newStrandID))		// Already paired with a different dicentric candidate earlier in the outer loop, skip.
            {
                continue;
            }

            if (!acentricCandidate.linear || !isRearrangementCandidate(acentricCandidate.newStrandID, numOriginalStrands)) // Must be a rearranged record and linear, otherwise skip.
            {
                continue;
            }

            std::set<int> acentricStrands;
            bool hasCentromereFlag = false;

            for (const SDRfragment& fragment : acentricCandidate.fragments)
            {
                acentricStrands.insert(fragment.oldStrandID);

                if (fragment.hasCentromere)
                {
                    hasCentromereFlag = true;
                    break;
                }
            }

	    // Must be fully acentric, and reference EXACTLY strandX and strandY, the two centromeres the dicentric record took.
            if (hasCentromereFlag || acentricStrands.size() != 2 || !acentricStrands.count(strandX) || !acentricStrands.count(strandY))
            {
                continue;
            }

            std::vector<SDRfragment> pooledX;						// Every fragment referencing strandX, from BOTH the dicentric and acentric candidates, checked below to see if together they tile strandX's full length.
            std::vector<SDRfragment> pooledY;						// Same, for strandY.

            for (const SDRfragment& fragment : dicentricCandidate.fragments)
            {
                (fragment.oldStrandID == strandX ? pooledX : pooledY).push_back(fragment);
            }

            for (const SDRfragment& fragment : acentricCandidate.fragments)
            {
                (fragment.oldStrandID == strandX ? pooledX : pooledY).push_back(fragment);
            }

            std::set<std::size_t> excludeIndices = {i, j};				// The dicentric and acentric records themselves are never valid deletion candidates to fill their own gaps.

            for (std::size_t r = 0; r < records.size(); ++r)
            {
                if (claimedExcisedRecords.count(records[r].newStrandID))		// Also exclude deletion records already claimed by an earlier, different pair.
                {
                    excludeIndices.insert(r);
                }
            }

            std::vector<std::size_t> usedExcisedX;

            if (!tryTileWithDeletions(pooledX, strandX, excludeIndices, usedExcisedX, allowUnmatchedGaps))	// Confirm strandX's fragments (plus any matching deletion records) fully reconstruct strandX's original length.
            {
                continue;
            }

            std::set<std::size_t> excludeForY = excludeIndices;

            for (std::size_t idx : usedExcisedX)					// Whatever strandX's tiling used can't be reused by strandY's tiling, each deletion record fills at most one gap, once.
            {
                excludeForY.insert(idx);
            }

            std::vector<std::size_t> usedExcisedY;

            if (!tryTileWithDeletions(pooledY, strandY, excludeForY, usedExcisedY, allowUnmatchedGaps))	// Same check for strandY.
            {
                continue;
            }

            overrides[dicentricCandidate.newStrandID] = strandX;			// Confirmed match: the dicentric record draws under strandX's slot
            overrides[acentricCandidate.newStrandID] = strandY;				// and the acentric record draws under strandY's slot.
            claimedAcentricRecords.insert(acentricCandidate.newStrandID);

            for (std::size_t idx : usedExcisedX)					// Mark every deletion record used by either strand's tiling as claimed, so a later pair can't reuse them.
            {
                claimedExcisedRecords.insert(records[idx].newStrandID);
            }

            for (std::size_t idx : usedExcisedY)
            {
                claimedExcisedRecords.insert(records[idx].newStrandID);
            }

            break;									// Found this dicentric candidate's match, stop searching acentric candidates for it.
        }
    }
    }


    return overrides;
}











// --------------------------------------------------------------- //
// Function to detect DICENTRIC CHROMOSOMES in a single cell       //
// --------------------------------------------------------------- //

// A dicentric chromosome is a single LINEAR record carrying exactly two centromeres from two DIFFERENT original strands, with no third strand
// mixed in, the same shape findChromoplexyAcentricRecombinations(), searches for as a dicentric CANDIDATE, but here purely for
// COUNTING purposes: this does not require a matching acentric partner to be found, since that pairing is a drawing-placement concern, not a
// detection concern. 

// However, findChromoplexyAcentricRecombinations (the drawing-side dicentric recognition) still requires exactly 2 total strands. 
// It needs that to find a clean, position-tiling acentric partner, and a third strand's material would break that proof. 
// So a record this detector now counts as dicentric (centromere pair + a third strand's acentric fragment) won't be recognized by the drawing 
// code as a clean dicentric/acentric pair. It'll fall through to the generic chromothripsis-style drawing fallback instead.
inline std::vector<SDRdicentricEvent> detectDicentrics(const SDRsubHeader& subHeader, int numOriginalStrands)
{
    std::vector<SDRdicentricEvent> dicentrics;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (!record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands))	// Must be a linear rearranged record.
        {
            continue;
        }

        std::vector<int> centromereStrandsInOrder;
        std::vector<int> involvedOldStrandIDs;			// Every distinct strand referenced, in first-appearance order - not just the two centromere-bearing ones.

        for (const SDRfragment& fragment : record.fragments)
        {
            bool alreadySeen = false;

            for (int strandID : involvedOldStrandIDs)
            {
                if (strandID == fragment.oldStrandID)
                {
                    alreadySeen = true;
                    break;
                }
            }

            if (!alreadySeen)
            {
                involvedOldStrandIDs.push_back(fragment.oldStrandID);
            }

            if (fragment.hasCentromere)
            {
                centromereStrandsInOrder.push_back(fragment.oldStrandID);
            }
        }

        if (centromereStrandsInOrder.size() != 2)		// Exactly 2 centromeres required for a dicentric - any number of OTHER, acentric strands may also be fused in.
        {
            continue;
        }

        const int strandX = centromereStrandsInOrder[0];
        const int strandY = centromereStrandsInOrder[1];

        if (strandX == strandY)				// The two centromeres must come from two DIFFERENT strands.
        {
            continue;
        }

        SDRdicentricEvent event{};
        event.newStrandID = record.newStrandID;
        event.centromereStrandAid = strandX;
        event.centromereStrandBid = strandY;
        event.involvedOldStrandIDs = involvedOldStrandIDs;

        dicentrics.push_back(event);
    }

    return dicentrics;
}













// Detects the 3-strand chromoplexy shape that findDicentricAcentricHomeOverrides() does NOT cover: strand A's centromeric piece fuses with strand B's centromeric
// piece (dicentric), but instead of strand B's OTHER piece pairing back up with strand A (which is what findDicentricAcentricHomeOverrides looks for), it
// instead fuses onto a THIRD strand C's centromeric piece. The result is monocentric overall (A keeps its own centromere via the dicentric record, C
// keeps its own via the second record), and strand B has been fully consumed, both of its pieces now live on strands A and C, with nothing left representing
// B on its own.
//
// Strand A's and strand C's own remaining acentric fragments (whatever did NOT fuse with B) are the two "loose ends" of this same rearrangement and are
// re-paired with each other here, to be drawn together in their own slot.
//
// Two scope limits, compared to findDicentricAcentricHomeOverrides: the dicentric record here must have EXACTLY 2 fragments (both centromeric, one per
// strand), not "2 centromeres plus any number of extra acentric fragments", so a chromoplexy fusion combined with an additional deletion on strand A or
// strand B isn't picked up here. And strand A's/strand C's own leftover piece must be a single, clean fragment (or a single already-combined record) that
// exactly completes their length, there's no equivalent of tryTileWithDeletions() here to let extra deletion records fill in additional gaps.
inline std::vector<SDRchromoplexyAcentricPlacement> findChromoplexyAcentricRecombinations(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader)
{
    std::vector<SDRchromoplexyAcentricPlacement> placements;
    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;
    const double tolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

    // Checks that fragment1 and fragment2 together exactly tile strandID's whole declared length with no gap and no overlap, i.e. that strandID
    // really has been fully consumed by exactly these two pieces.
    auto fragmentsSpanFullLength = [&](const SDRfragment& fragment1, const SDRfragment& fragment2, int strandID) -> bool
    {
        const std::size_t sizeIndex = static_cast<std::size_t>(strandID);

        if (sizeIndex >= masterHeader.intactChromosomeSizes.size())				// Only care about intact chromosome sizes
        {
            return false;
        }

        const double fullSize = masterHeader.intactChromosomeSizes[sizeIndex];

        if (fullSize <= 0.0)
        {
            return false;
        }

        double lowerStart = std::min(fragment1.oldStartPosition, fragment1.oldEndPosition);
        double lowerEnd = std::max(fragment1.oldStartPosition, fragment1.oldEndPosition);
        double higherStart = std::min(fragment2.oldStartPosition, fragment2.oldEndPosition);
        double higherEnd = std::max(fragment2.oldStartPosition, fragment2.oldEndPosition);

        if (lowerStart > higherStart)
        {
            std::swap(lowerStart, higherStart);
            std::swap(lowerEnd, higherEnd);
        }

        return approxEqual(lowerStart, 0.0, tolerance) && approxEqual(lowerEnd, higherStart, tolerance) && approxEqual(higherEnd, fullSize, tolerance);
    };

    for (std::size_t i = 0; i < records.size(); ++i)					// Outer loop: try every record as a potential DICENTRIC candidate (the strand A + strand B fusion).
    {
        const SDRdataRecord& dicentric = records[i];

        if (!dicentric.linear || !isRearrangementCandidate(dicentric.newStrandID, numOriginalStrands) || dicentric.fragments.size() != 2)	// Exactly 2 fragments required, see the scope-limit note above.
        {
            continue;
        }

        std::vector<int> centromereStrands;						// Strand IDs of every centromere-bearing fragment in this record, in the order they appear.

        for (const SDRfragment& fragment : dicentric.fragments)
        {
            if (fragment.hasCentromere)
            {
                centromereStrands.push_back(fragment.oldStrandID);
            }
        }

        if (centromereStrands.size() != 2 || centromereStrands[0] == centromereStrands[1])
        {
            continue; 									// Not a two-different-strand dicentric fusion.
        }

        const int strandA = centromereStrands[0];					// First-listed centromere-bearing strand.
        const int strandB = centromereStrands[1];					// Second-listed centromere-bearing strand, the one that gets fully consumed.


        // This function assigns no explicit home-strand override for the dicentric record itself, unlike  findDicentricAcentricHomeOverrides(). It relies on 
	// determineHomeStrandID()'s own fallback rule (first fragment in the record, used when there are 2+ centromere flags) landing on strandA. That
        // holds here specifically because dicentric.fragments has EXACTLY 2 entries and BOTH are centromere-flagged (enforced above), so "first fragment positionally" 
	// and "first centromere-bearing fragment" (=strandA) are the same thing.
        const SDRfragment& fragmentA = (dicentric.fragments[0].oldStrandID == strandA) ? dicentric.fragments[0] : dicentric.fragments[1];
        const SDRfragment& fragmentB = (dicentric.fragments[0].oldStrandID == strandB) ? dicentric.fragments[0] : dicentric.fragments[1];

        // Look for a second record fusing strand B's OTHER (acentric) piece onto a third strand C's centromere.
        for (std::size_t j = 0; j < records.size(); ++j)
        {
            if (i == j)
            {
                continue;
            }

            const SDRdataRecord& monocentric = records[j];

            if (!monocentric.linear || !isRearrangementCandidate(monocentric.newStrandID, numOriginalStrands) || monocentric.fragments.size() != 2)
            {
                continue;
            }

            int strandC = -1;								// Resolved below; unlike strandA/strandB, a record with exactly ONE centromere flag (this one) is unambiguous, so determineHomeStrandID() lands on strandC directly without relying on the fallback rule.
            int centromereCountHere = 0;
            const SDRfragment* fragmentBOther = nullptr;
            const SDRfragment* fragmentC = nullptr;

            for (const SDRfragment& fragment : monocentric.fragments)
            {
                if (fragment.hasCentromere)
                {
                    ++centromereCountHere;

                    if (fragment.oldStrandID != strandB)
                    {
                        strandC = fragment.oldStrandID;
                        fragmentC = &fragment;
                    }
                }

                if (fragment.oldStrandID == strandB)
                {
                    fragmentBOther = &fragment;
                }
            }

            // Must be monocentric overall (strand C keeps exactly one centromere), must actually reference strand B, that B fragment must be the acentric half, and C must be a genuinely
            // third strand, not A or B again.
            if (centromereCountHere != 1 || strandC == -1 || strandC == strandA || strandC == strandB || fragmentBOther == nullptr || fragmentBOther->hasCentromere || fragmentC == nullptr)
            {
                continue;
            }

            // Confirm strand B really is fully consumed: its two pieces (one from each record) tile its whole length.
            if (!fragmentsSpanFullLength(fragmentB, *fragmentBOther, strandB))
            {
                continue;
            }

            // Find strand A's and strand C's own leftover acentric fragments, whatever did NOT fuse with strand B, and confirm each, together with its partner's centromeric fragment from
            // above, fully tiles its own strand's length. The SDR file may express these two leftover pieces either as ONE record that already combines both (fragments.size() == 2, one
            // fragment per strand) or as two separate single-fragment records, one per strand.
            std::vector<int> leftoverNewStrandIDs;

            for (const SDRdataRecord& candidate : records)
            {
                if (candidate.newStrandID == dicentric.newStrandID || candidate.newStrandID == monocentric.newStrandID)
                {
                    continue;
                }

                if (!candidate.linear || !isRearrangementCandidate(candidate.newStrandID, numOriginalStrands))
                {
                    continue;
                }

                if (candidate.fragments.size() == 2)
                {
                    // Already-combined shape: one fragment must be strand A's leftover, the other strand C's leftover, neither carrying a centromere.
                    const SDRfragment* candidateA = nullptr;
                    const SDRfragment* candidateC = nullptr;

                    for (const SDRfragment& fragment : candidate.fragments)
                    {
                        if (fragment.hasCentromere)
                        {
                            candidateA = nullptr;
                            candidateC = nullptr;
                            break;
                        }

                        if (fragment.oldStrandID == strandA)
                        {
                            candidateA = &fragment;
                        }
                        else if (fragment.oldStrandID == strandC)
                        {
                            candidateC = &fragment;
                        }
                    }

                    if (candidateA != nullptr && candidateC != nullptr && fragmentsSpanFullLength(fragmentA, *candidateA, strandA) && fragmentsSpanFullLength(*fragmentC, *candidateC, strandC))
                    {
                        leftoverNewStrandIDs = { candidate.newStrandID };
                        break;							// First qualifying combined record wins, the data isn't expected to contain more than one valid match here.
                    }
                }
            }

            if (leftoverNewStrandIDs.empty())
            {
                // Fall back to the two-separate-records shape.
                const SDRdataRecord* leftoverA = nullptr;
                const SDRdataRecord* leftoverC = nullptr;

                for (const SDRdataRecord& candidate : records)
                {
                    if (candidate.newStrandID == dicentric.newStrandID || candidate.newStrandID == monocentric.newStrandID)
                    {
                        continue;
                    }

                    if (!candidate.linear || !isRearrangementCandidate(candidate.newStrandID, numOriginalStrands) || candidate.fragments.size() != 1)
                    {
                        continue;
                    }

                    const SDRfragment& fragment = candidate.fragments[0];

                    if (fragment.hasCentromere)
                    {
                        continue;
                    }

                    if (leftoverA == nullptr && fragment.oldStrandID == strandA && fragmentsSpanFullLength(fragmentA, fragment, strandA))
                    {
                        leftoverA = &candidate;
                    }
                    else if (leftoverC == nullptr && fragment.oldStrandID == strandC && fragmentsSpanFullLength(*fragmentC, fragment, strandC))
                    {
                        leftoverC = &candidate;
                    }
                }

                if (leftoverA != nullptr && leftoverC != nullptr)
                {
                    leftoverNewStrandIDs = { leftoverA->newStrandID, leftoverC->newStrandID };
                }
            }

            if (leftoverNewStrandIDs.empty())						// Neither shape matched. Strand A's and strand C's own leftovers don't cleanly account for the rest of their length, so this isn't a valid chromoplexy placement.
            {
                continue;
            }

            SDRchromoplexyAcentricPlacement placement{};
            placement.consumedStrandID = strandB;					// Strand B is hidden entirely, it has no material of its own left.
            placement.leftoverNewStrandIDs = leftoverNewStrandIDs;			// Strand A's and strand C's own leftover piece(s), re-paired into one combined drawing slot.
            placements.push_back(placement);

            break; 									// Found this dicentric record's match, move on to the next candidate dicentric record.
        }
    }

    return placements;
}


















// =========================================================================== //
// SHARED HELPERS FOR THE MUTATION HIERARCHY (simple < more complex < complex genomic rearrangements)
// =========================================================================== //

// Forward declarations (defaults live here only, NOT on the later definitions).
inline std::vector<SDRchromoplexyEvent> detectChromoplexy(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader, bool collectAccompanyingDeletions = true);
inline std::vector<SDRchromothripsisEvent> detectChromothripsis(const SDRsubHeader& subHeader, int numOriginalStrands, int minFragmentCount = 10);
inline std::vector<SDRcomplexRearrangementEvent> detectComplexGenomicRearrangements(const SDRsubHeader& subHeader, int numOriginalStrands, int minFragmentCount = 10);
inline std::vector<SDRdeletionTranslocationEvent> detectDeletionTranslocations(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader, bool applyHierarchy = true);
inline std::vector<SDRtranslocationEvent> detectTranslocations(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader, bool applyHierarchy = true);
inline std::vector<SDRdeletionInsertionEvent> detectDeletionInsertions(const SDRsubHeader& subHeader, int numOriginalStrands, bool applyHierarchy = true);

// A fragment with its positions normalized so low <= high, remembering whether it was reversed (start > end).
struct SDRnormalizedFragment
{
    double low;
    double high;
    bool reversed;
};

// A gap between two neighbouring fragments of one old strand. leftIndex/rightIndex index the SORTED normalized fragments bordering the gap; -1 marks a terminal gap with no neighbour on that side.
struct SDRflankGap
{
    double start;
    double end;
    int leftIndex;
    int rightIndex;
};

// Normalizes + sorts the fragments of ONE old strand and finds every gap (leading/internal/trailing). Leading/trailing gaps are only searched when fullSizeMbp > 0.
// Returns false if two fragments overlap by more than tolerance (invalid shape).
inline bool analyzeFlankingFragments(const std::vector<SDRfragment>& fragments, double fullSizeMbp, double tolerance, std::vector<SDRnormalizedFragment>& sorted, std::vector<SDRflankGap>& gaps)
{
    sorted.clear();
    gaps.clear();

    if (fragments.empty())
    {
        return false;
    }

    for (const SDRfragment& fragment : fragments)
    {
        sorted.push_back({std::min(fragment.oldStartPosition, fragment.oldEndPosition), std::max(fragment.oldStartPosition, fragment.oldEndPosition), fragment.oldStartPosition > fragment.oldEndPosition});
    }

    std::sort(sorted.begin(), sorted.end(), [](const SDRnormalizedFragment& a, const SDRnormalizedFragment& b)
    {
        if (a.low != b.low)
        {
            return a.low < b.low;
        }
        return a.high < b.high;
    });

    const int count = static_cast<int>(sorted.size());

    if (fullSizeMbp > 0.0 && sorted.front().low > 0.0 && !approxEqual(sorted.front().low, 0.0, tolerance))
    {
        gaps.push_back({0.0, sorted.front().low, -1, 0});
    }

    for (int i = 0; i + 1 < count; ++i)
    {
        const double higher = sorted[i].high;
        const double lower = sorted[i + 1].low;

        if (approxEqual(higher, lower, tolerance))
        {
            continue;
        }

        if (lower < higher)								// Overlap = invalid shape
        {
            return false;
        }

        gaps.push_back({higher, lower, i, i + 1});
    }

    if (fullSizeMbp > 0.0 && fullSizeMbp > sorted.back().high && !approxEqual(sorted.back().high, fullSizeMbp, tolerance))
    {
        gaps.push_back({sorted.back().high, fullSizeMbp, count - 1, -1});
    }

    return true;
}

// True if the fragment on either side of the gap is reversed (a terminal gap only has one neighbour).
inline bool gapTouchesReversedFragment(const std::vector<SDRnormalizedFragment>& sorted, const SDRflankGap& gap)
{
    return (gap.leftIndex >= 0 && sorted[gap.leftIndex].reversed) || (gap.rightIndex >= 0 && sorted[gap.rightIndex].reversed);
}

// What (if anything) the material missing from a gap turned into elsewhere in the cell.
enum class SDRgapPartner
{
    None,			// Lone gap, material is simply gone
    LinearExcised,		// A single, non-reversed, linear, one-fragment record exactly filling the gap
    Circular,			// Reappears in a circular record -> ecDNA formation
    Foreign			// Reappears inside a linear multi-fragment record -> insertion / translocation / chromoplexy signature
};

// When a gap's material broke into SEVERAL circular pieces (e.g. a centric ring plus an acentric ecDNA cut from the same region),
// returns the pieces' [low, high] ranges in order, provided they tile the gap exactly (contiguous within tolerance, nothing left over).
// Returns an empty vector otherwise, including when a single piece fills the whole gap, which the exact-match checks already handle.
inline std::vector<std::pair<double, double>> circularPiecesTilingGap(const SDRsubHeader& subHeader, int numOriginalStrands, int oldStrandID, double gapStart, double gapEnd, double tolerance)
{
    std::vector<std::pair<double, double>> pieces;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands))
        {
            continue;
        }

        for (const SDRfragment& fragment : record.fragments)
        {
            const double low = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
            const double high = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

            if (fragment.oldStrandID == oldStrandID && low >= gapStart - tolerance && high <= gapEnd + tolerance)
            {
                pieces.push_back({low, high});
            }
        }
    }

    if (pieces.size() < 2)
    {
        return {};
    }

    std::sort(pieces.begin(), pieces.end());

    double reached = gapStart;

    for (const std::pair<double, double>& piece : pieces)
    {
        if (!approxEqual(piece.first, reached, tolerance))			// Overlap or uncovered stretch, not a clean split of this gap.
        {
            return {};
        }

        reached = piece.second;
    }

    return approxEqual(reached, gapEnd, tolerance) ? pieces : std::vector<std::pair<double, double>>{};
}

// Looks through every rearranged record (except ignoreNewStrandIDs) for a fragment of oldStrandID matching [gapStart, gapEnd] within tolerance.
// Priority: Circular > Foreign > LinearExcised > None. excisedStrandID is set only for LinearExcised (else -1).
// acceptReversedExcised: a single-fragment linear record filling the gap counts as the excised piece in either orientation (long deletions, deletion-inversions
// and deletion-translocations, Samy 2026-11-03: the cut-out piece is loose, so an inverted listing such as 5/60/40/0 is still just the deleted segment, as long as
// it is not inserted into another chromosome). Otherwise only a forward one does and a reversed one is Foreign.
inline SDRgapPartner classifyGapPartner(const SDRsubHeader& subHeader, int numOriginalStrands, int oldStrandID, double gapStart, double gapEnd, double tolerance, const std::set<int>& ignoreNewStrandIDs, int& excisedStrandID,
                                        bool acceptReversedExcised = false)
{
    excisedStrandID = -1;

    bool circular = false;
    bool foreign = false;
    int linearID = -1;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands) || ignoreNewStrandIDs.count(record.newStrandID))
        {
            continue;
        }

        for (const SDRfragment& fragment : record.fragments)
        {
            if (fragment.oldStrandID != oldStrandID)
            {
                continue;
            }

            const double low = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
            const double high = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

            if (!approxEqual(low, gapStart, tolerance) || !approxEqual(high, gapEnd, tolerance))
            {
                continue;
            }

            if (!record.linear)
            {
                circular = true;
            }
            else if (record.fragments.size() == 1 && (acceptReversedExcised || fragment.oldStartPosition <= fragment.oldEndPosition))
            {
                if (linearID == -1)
                {
                    linearID = record.newStrandID;
                }
            }
            else
            {
                foreign = true;
            }
        }
    }

    // A gap whose material became several circular pieces (a ring plus ecDNA, or several ecDNA) is a circular formation too, not a deletion.
    if (circular || !circularPiecesTilingGap(subHeader, numOriginalStrands, oldStrandID, gapStart, gapEnd, tolerance).empty())
    {
        return SDRgapPartner::Circular;
    }

    if (foreign)
    {
        return SDRgapPartner::Foreign;
    }

    if (linearID != -1)
    {
        excisedStrandID = linearID;
        return SDRgapPartner::LinearExcised;
    }

    return SDRgapPartner::None;
}

// For every old strand: the linear, rearranged, PURE (all fragments from that strand), 2+ fragment records ("flanking/remaining" records).
inline std::map<int, std::vector<const SDRdataRecord*>> collectFlankingRecords(const SDRsubHeader& subHeader, int numOriginalStrands)
{
    std::map<int, std::vector<const SDRdataRecord*>> result;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (!record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands) || record.fragments.size() < 2)
        {
            continue;
        }

        bool pure = true;
        for (const SDRfragment& fragment : record.fragments)
        {
            if (fragment.oldStrandID != record.fragments[0].oldStrandID)
            {
                pure = false;
                break;
            }
        }

        if (pure)
        {
            result[record.fragments[0].oldStrandID].push_back(&record);
        }
    }

    // TERMINAL deletions leave a SINGLE-fragment remaining record (no flanking pair). Such a record is only accepted for a strand that has no multi-fragment
    // remaining record, when exactly one of that strand's pure single-fragment linear records carries the centromere and is not reversed (the excised tip is acentric).
    std::map<int, std::vector<const SDRdataRecord*>> singleFragmentCentric;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (!record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands) || record.fragments.size() != 1)
        {
            continue;
        }

        const SDRfragment& fragment = record.fragments[0];

        if (fragment.hasCentromere && fragment.oldStartPosition <= fragment.oldEndPosition)
        {
            singleFragmentCentric[fragment.oldStrandID].push_back(&record);
        }
    }

    for (const auto& [strandID, candidates] : singleFragmentCentric)
    {
        if (candidates.size() == 1 && result.find(strandID) == result.end())
        {
            result[strandID].push_back(candidates[0]);
        }
    }

    return result;
}

// Chromosome's declared full size in Mbp (0 if unknown); intactChromosomeSizes index == oldStrandID.
inline double fullSizeForStrand(const SDRmasterHeader& masterHeader, int oldStrandID)
{
    const std::size_t index = static_cast<std::size_t>(oldStrandID);
    return (index < masterHeader.intactChromosomeSizes.size()) ? masterHeader.intactChromosomeSizes[index] : 0.0;
}

// Strands claimed by the complex genomic rearrangements that silence everything simpler: chromothripsis (10+ fragments in a 1-2 strand cluster), complex genomic
// rearrangements (4+ pieces of one chromosome out of order, or shattered exchanged material) and chromoplexy (3+ strand cluster).
// ecDNA is handled per-gap through the Circular partner check, so it is not part of this set.
inline std::set<int> collectHigherOrderStrands(const SDRsubHeader& subHeader, int numOriginalStrands)
{
    std::set<int> strands;

    for (const SDRchromothripsisEvent& event : detectChromothripsis(subHeader, numOriginalStrands))
    {
        strands.insert(event.involvedStrandIDs.begin(), event.involvedStrandIDs.end());
    }

    for (const SDRcomplexRearrangementEvent& event : detectComplexGenomicRearrangements(subHeader, numOriginalStrands))	// Silences the same simpler events chromothripsis did
    {
        strands.insert(event.involvedStrandIDs.begin(), event.involvedStrandIDs.end());
    }

    const SDRmasterHeader emptyHeader{};
    for (const SDRchromoplexyEvent& event : detectChromoplexy(subHeader, numOriginalStrands, emptyHeader, false))	// false avoids recursing back into detectDeletions
    {
        strands.insert(event.involvedStrandIDs.begin(), event.involvedStrandIDs.end());
    }

    return strands;
}


// --------------------------------------------------------------------------- //
// Detects LONG DELETION events within a single cell.
// --------------------------------------------------------------------------- //

// Modified detectDeletions function to detect multiple deletions within a single strand. A strand with N deletions will contain N + 1 fragments in SDR data field 3, and will have
// N + 1 data entries, 1 entry being the original strand with all the gaps, and N entries representing each deletion causing the gaps.
// Function also handles terminal chromosome deletions at the beginning or end of chromosomes, so the deleted segment will lack two flanking fragments.
inline std::vector<SDRdeletionEvent> detectDeletions(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader, bool applyHierarchy = true)
{
    const double delTolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

    std::vector<SDRdeletionEvent> deletions;

    for (const auto& [oldStrandID, records] : collectFlankingRecords(subHeader, numOriginalStrands))
    {
        if (records.size() != 1)							// Exactly one remaining record per strand, otherwise ambiguous.
        {
            continue;
        }

        const SDRdataRecord& flanking = *records[0];

        std::vector<SDRnormalizedFragment> sorted;
        std::vector<SDRflankGap> gaps;

        if (!analyzeFlankingFragments(flanking.fragments, fullSizeForStrand(masterHeader, oldStrandID), delTolerance, sorted, gaps))
        {
            continue;
        }

        for (const SDRflankGap& gap : gaps)
        {
            if (gapTouchesReversedFragment(sorted, gap))				// Belongs to a deletion-inversion (or is skipped as too complex).
            {
                continue;
            }

            int excisedStrandID = -1;
            const SDRgapPartner partner = classifyGapPartner(subHeader, numOriginalStrands, oldStrandID, gap.start, gap.end, delTolerance, {flanking.newStrandID}, excisedStrandID, true);

            if (partner == SDRgapPartner::Circular || partner == SDRgapPartner::Foreign)	// ecDNA / deletion-insertion / translocation signature, not a plain deletion.
            {
                continue;
            }

            SDRdeletionEvent event{};
            event.oldStrandID = oldStrandID;
            event.deletionStart = gap.start;
            event.deletionEnd = gap.end;
            event.remainingStrandID = flanking.newStrandID;
            event.excisedStrandID = excisedStrandID;					// -1 when there is no excised record.

            deletions.push_back(event);
        }
    }

    // Strands taking part in a deletion-translocation or deletion-insertion keep their remaining material in MIXED records (foreign fragments next to their own),
    // which the pure-record pass above cannot see. Pool each such strand's fragments over every linear multi-fragment rearranged record and report the gaps that
    // belong to no larger event as ordinary deletions: not the translocation's own deletion, not an insertion's donor gap (a Foreign partner), not an inversion's flank.
    {
        const std::map<int, std::vector<const SDRdataRecord*>> pureRecords = collectFlankingRecords(subHeader, numOriginalStrands);

        std::set<int> mixedStrands;
        std::vector<std::array<double, 3>> claimedGaps;				// {strandID, start, end} of every deletion-translocation's own deletion
        const SDRmasterHeader& header = masterHeader;

        for (const SDRdeletionTranslocationEvent& event : detectDeletionTranslocations(subHeader, numOriginalStrands, header, false))
        {
            mixedStrands.insert(event.oldStrandAid);
            mixedStrands.insert(event.oldStrandBid);
            claimedGaps.push_back({static_cast<double>(event.deletedOldStrandID), event.deletionStart, event.deletionEnd});
        }

        for (const SDRtranslocationEvent& event : detectTranslocations(subHeader, numOriginalStrands, header, false))
        {
            mixedStrands.insert(event.oldStrandAid);
            mixedStrands.insert(event.oldStrandBid);
        }

        for (const SDRdeletionInsertionEvent& event : detectDeletionInsertions(subHeader, numOriginalStrands, false))
        {
            mixedStrands.insert(event.recipientOldStrandID);
        }

        for (int oldStrandID : mixedStrands)
        {
            if (pureRecords.count(oldStrandID))					// Already handled by the pure-record pass.
            {
                continue;
            }

            std::vector<SDRfragment> pooled;
            std::vector<int> pooledRecordIDs;

            for (const SDRdataRecord& record : subHeader.dataRecords)
            {
                if (!record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands) || record.fragments.size() < 2)
                {
                    continue;
                }

                for (const SDRfragment& fragment : record.fragments)
                {
                    if (fragment.oldStrandID == oldStrandID)
                    {
                        pooled.push_back(fragment);
                        pooledRecordIDs.push_back(record.newStrandID);
                    }
                }
            }

            std::vector<SDRnormalizedFragment> sorted;
            std::vector<SDRflankGap> gaps;

            if (pooled.empty() || !analyzeFlankingFragments(pooled, fullSizeForStrand(masterHeader, oldStrandID), delTolerance, sorted, gaps))
            {
                continue;
            }

            for (const SDRflankGap& gap : gaps)
            {
                bool claimed = false;
                for (const std::array<double, 3>& claim : claimedGaps)
                {
                    claimed = claimed || (static_cast<int>(claim[0]) == oldStrandID && approxEqual(claim[1], gap.start, delTolerance) && approxEqual(claim[2], gap.end, delTolerance));
                }

                if (claimed || gapTouchesReversedFragment(sorted, gap))
                {
                    continue;
                }

                // The record holding the neighbouring fragment becomes the deletion's "remaining" record.
                int remainingID = -1;
                for (std::size_t k = 0; k < pooled.size() && remainingID == -1; ++k)
                {
                    const double low = std::min(pooled[k].oldStartPosition, pooled[k].oldEndPosition);
                    const double high = std::max(pooled[k].oldStartPosition, pooled[k].oldEndPosition);

                    if ((gap.leftIndex >= 0 && approxEqual(high, gap.start, delTolerance)) || (gap.leftIndex < 0 && approxEqual(low, gap.end, delTolerance)))
                    {
                        remainingID = pooledRecordIDs[k];
                    }
                }

                int excisedStrandID = -1;
                const SDRgapPartner partner = classifyGapPartner(subHeader, numOriginalStrands, oldStrandID, gap.start, gap.end, delTolerance, {}, excisedStrandID, true);

                if (partner == SDRgapPartner::Circular || partner == SDRgapPartner::Foreign)
                {
                    continue;
                }

                SDRdeletionEvent event{};
                event.oldStrandID = oldStrandID;
                event.deletionStart = gap.start;
                event.deletionEnd = gap.end;
                event.remainingStrandID = remainingID;
                event.excisedStrandID = excisedStrandID;

                deletions.push_back(event);
            }
        }
    }

    if (applyHierarchy)									// Chromothripsis / chromoplexy silence every simple deletion on their strands.
    {
        const std::set<int> silenced = collectHigherOrderStrands(subHeader, numOriginalStrands);
        deletions.erase(std::remove_if(deletions.begin(), deletions.end(), [&](const SDRdeletionEvent& e) { return silenced.count(e.oldStrandID) > 0; }), deletions.end());
    }

    return deletions;
}













// Helper function to check for candidate balanced inversion fragments, where the start and end positions of a fragment are flipped.
inline bool isReversedFragment(const SDRfragment& fragment)
{
    return fragment.oldStartPosition > fragment.oldEndPosition;
}


// A record can list a chromosome from its other end: exactly three pieces of one chromosome, the two outer ones reversed and reaching opposite chromosome ends
// (the start and end of the chromosome swap places), the middle one forward, e.g. 24/132.6/52.9, 24/21.7/52.9, 24/21.7/0. Read from the other end that is
// 0-21.7, 52.9-21.7 (inverted), 52.9-132.6: only the middle piece is inverted. Returns the fragments in that reading (order reversed, each start/end swapped),
// or the record's own fragments when it is not listed upside down, so the detectors score both listings the same way.
inline std::vector<SDRfragment> fragmentsInReadingOrder(const SDRdataRecord& record, const SDRmasterHeader& masterHeader)
{
    const double tolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)
    const std::vector<SDRfragment>& fragments = record.fragments;

    if (!record.linear || fragments.size() != 3 || fragments[0].oldStrandID != fragments[1].oldStrandID || fragments[1].oldStrandID != fragments[2].oldStrandID)
    {
        return fragments;
    }

    const double strandSize = fullSizeForStrand(masterHeader, fragments[0].oldStrandID);
    const bool upsideDown = strandSize > 0.0
                            && isReversedFragment(fragments[0]) && !isReversedFragment(fragments[1]) && isReversedFragment(fragments[2])
                            && approxEqual(fragments[0].oldStartPosition, strandSize, tolerance) && approxEqual(fragments[2].oldEndPosition, 0.0, tolerance);

    if (!upsideDown)
    {
        return fragments;
    }

    std::vector<SDRfragment> flipped(fragments.rbegin(), fragments.rend());

    for (SDRfragment& fragment : flipped)
    {
        std::swap(fragment.oldStartPosition, fragment.oldEndPosition);
    }

    return flipped;
}


// ---------------------------------------------------------------------------- //
// Detects BALANCED INVERSIONS within a single cell.
// ---------------------------------------------------------------------------- //

// A balanced inversion is recognized within a single new-strand record: fragments referencing the same old strand ID must form an
// intact region with no gaps (balanced) and can have more than one reverse fragment.
inline std::vector<SDRinversionEvent> detectInversions(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader, bool applyHierarchy = true)
{
    const double balInvTolerance = svToleranceMbp;				// Up to 99 bases lost at the breakpoints is still balanced

    std::vector<SDRinversionEvent> inversions;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (!record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands))
        {
            continue;
        }

        std::map<int, std::vector<SDRfragment>> groupOldStrand;

        for (const SDRfragment& fragment : fragmentsInReadingOrder(record, masterHeader))	// A pericentric inversion listed upside down reads as one inner inversion
        {
            groupOldStrand[fragment.oldStrandID].push_back(fragment);
        }

        for (auto& [oldStrandID, fragments] : groupOldStrand)
        {
            if (fragments.size() < 2)
            {
                continue;
            }

            bool anyReversed = false;
            for (const SDRfragment& fragment : fragments)
            {
                anyReversed = anyReversed || isReversedFragment(fragment);
            }

            if (!anyReversed)
            {
                continue;
            }

            // Terminal gaps are only meaningful when the group is the whole record (a pure record).
            const double fullSize = (fragments.size() == record.fragments.size()) ? fullSizeForStrand(masterHeader, oldStrandID) : 0.0;

            std::vector<SDRnormalizedFragment> sorted;
            std::vector<SDRflankGap> gaps;

            if (!analyzeFlankingFragments(fragments, fullSize, balInvTolerance, sorted, gaps))
            {
                continue;
            }

            for (int index = 0; index < static_cast<int>(sorted.size()); ++index)
            {
                if (!sorted[index].reversed)
                {
                    continue;
                }

                bool touchesGap = false;
                for (const SDRflankGap& gap : gaps)
                {
                    if (gap.leftIndex == index || gap.rightIndex == index)
                    {
                        touchesGap = true;
                        break;
                    }
                }

                if (touchesGap)								// Deleted material next to the inversion -> deletion-inversion or skipped, not balanced.
                {
                    continue;
                }

                // A balanced inversion needs a break on BOTH sides of the inverted segment. A reversed segment that reaches either end of
                // the chromosome (a terminal inversion) is not scored as one; countMutationsPerOriginalStrand() counts it as a Complex Genomic Rearrangement.
                const double strandSize = fullSizeForStrand(masterHeader, oldStrandID);

                if (approxEqual(sorted[index].low, 0.0, balInvTolerance) || (strandSize > 0.0 && approxEqual(sorted[index].high, strandSize, balInvTolerance)))
                {
                    continue;
                }

                SDRinversionEvent event{};
                event.oldStrandID = oldStrandID;
                event.inversionStart = sorted[index].low;
                event.inversionEnd = sorted[index].high;
                event.newStrandID = record.newStrandID;

                inversions.push_back(event);
            }
        }
    }

    if (applyHierarchy)
    {
        const std::set<int> silenced = collectHigherOrderStrands(subHeader, numOriginalStrands);
        inversions.erase(std::remove_if(inversions.begin(), inversions.end(), [&](const SDRinversionEvent& e) { return silenced.count(e.oldStrandID) > 0; }), inversions.end());
    }

    return inversions;
}












// A candidate "half" of a balanced translocation: a record with exactly two fragments, referencing two distinct old strand IDs,
// both fragments in their original (non-reversed) orientation.
inline bool getTranslocationCandidateFragments(const SDRdataRecord& record, int numOriginalStrands, SDRfragment& fragmentA, SDRfragment& fragmentB)
{

    if (!record.linear)									// Fragments must be linear in balanced translocations
    {
	return false;
    }

    if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))		// Make sure we are skipping the intact strands for mutation detection
    {
	return false;
    }

    if (record.fragments.size() != 2)							// Each resulting new data record must contain exactly two fragments each, otherwise it is not a balanced translocation.
    {
        return false;
    }

    const SDRfragment& first = record.fragments[0];					// Track first balanced translocation fragment
    const SDRfragment& second = record.fragments[1];					// Track second balanced translocation fragment

    if (first.oldStrandID == second.oldStrandID)					// Balanced translocations occur on two different intact/old strand IDs, cannot have them on the same strand ID.
    {
        return false;
    }

    if (isReversedFragment(first) || isReversedFragment(second))			// Check if fragments are inverted, if so, not a balanced translocation mutation
    {
        return false;
    }

    fragmentA = first;
    fragmentB = second;

    return true;
}











// Checks whether two fragments referencing the same old strand ID meet exactly at a shared boundary within tolerance, i.e. they
// are the two pieces resulting from a single break in that strand. If so, populates breakpoint with the position of that break.
inline bool fragmentsShareBreakpoint(const SDRfragment& fragment1, const SDRfragment& fragment2, double& breakpoint)
{
    double balTransTolerance = svToleranceMbp; 					// Balanced translocations can lose up to 99 bases at the breakpoint (see svToleranceMbp)

    if (fragment1.oldStrandID != fragment2.oldStrandID)					// Both fragments must reference the same original strand. Checks whether they're the two pieces of a single break in that strand.
    {
        return false;
    }

    if (isReversedFragment(fragment1) || isReversedFragment(fragment2))			// Balanced translocation fragments cannot be reversed.
    {
        return false;
    }

    if (approxEqual(fragment1.oldEndPosition, fragment2.oldStartPosition, balTransTolerance)) // If fragment start and ends are within tolerance, record the location of the break for the first strands.
    {
        breakpoint = fragment1.oldEndPosition;
        return true;
    }

    if (approxEqual(fragment2.oldEndPosition, fragment1.oldStartPosition, balTransTolerance)) // If fragment start and ends are within tolerance, record the location of the break for the second strand.
    {
        breakpoint = fragment2.oldEndPosition;
        return true;
    }

    return false;
}













// ------------------------------------------------------------------------------------	//
// Function to detect BALANCED TRANSLOCATION entries in SDR file
// ------------------------------------------------------------------------------------ //


// Generalized to detect translocations anywhere on a chromosome (not just at its ends), and to correctly handle an exchanged piece that re-attaches in reversed orientation
// (which looks locally like an inversion but isn't one, since it spans two strands). Pools every fragment referencing a given old strand across BOTH candidate records,
// normalizes each to [min, max] (so reversed fragments sort correctly).
inline std::vector<SDRtranslocationEvent> detectTranslocations(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader, bool applyHierarchy)
{
    std::vector<SDRtranslocationEvent> translocations;					// Translocations vector to store translocation events

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;			// Records vector to store the SDR data records/entries

    const double balTraTolerance = svToleranceMbp;					// Up to 99 bases lost at a breakpoint is still balanced (see svToleranceMbp)


    // Sub function to check if the strand being investigated is valid (linear for BalTrans, and a mutation record). Need this since we are
    // checking across multiple data records for two distinct old strand IDs sharing DNA information in a balanced manner (no information loss).
    // Groups a given record's fragments by each old Strand ID that fragment came from
    auto getRecordStrandFragments = [&](const SDRdataRecord& record, std::map<int, std::vector<SDRfragment>>& groupStrand) -> bool
    {
        if (!record.linear)								// Data record must be linear to be considered
        {
            return false;
        }

        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))		// Data record must have a new strand ID > numOriginalStrands
        {										// Otherwise not a mutation
            return false;
        }

        for (const SDRfragment& fragment : record.fragments)				// If initial checks pass, store the fragment associated with
        {										// A given oldStrandID
            groupStrand[fragment.oldStrandID].push_back(fragment);
        }

        // A record holding 3+ chromosomes can carry a translocation AND an unrelated insertion, e.g. 38/0/119.4, 30/0.8/2.2, 38/119.4/130.2, 18/11.4/132.6:
        // chr 30's piece is only inserted into chr 38 (a deletion-insertion, found by detectDeletionInsertions()), while chr 38 and chr 18 exchange ends.
        // Such inserted pieces (flanked on both sides by one other chromosome) are set aside so the translocation between the other two can still be checked.
        // Safeguard: only done for 3+ chromosomes, so a two-chromosome central-segment exchange (A, B, A) keeps its inserted-looking piece.
        if (groupStrand.size() > 2)
        {
            groupStrand.clear();

            for (std::size_t k = 0; k < record.fragments.size(); ++k)
            {
                const bool inserted = k > 0 && k + 1 < record.fragments.size()
                                      && record.fragments[k - 1].oldStrandID == record.fragments[k + 1].oldStrandID
                                      && record.fragments[k - 1].oldStrandID != record.fragments[k].oldStrandID;

                if (!inserted)
                {
                    groupStrand[record.fragments[k].oldStrandID].push_back(record.fragments[k]);
                }
            }
        }

        return !groupStrand.empty();							// Return true if the groupStrand is non-empty
    };



    for (std::size_t i = 0; i < records.size(); ++i)				// Loop through all SDR data records, looking for a first candidate
    {
        std::map<int, std::vector<SDRfragment>> firstStrandGroup;		// This variable stores the first group of two original strands contributing to the 
										// balanced translocation.
        if (!getRecordStrandFragments(records[i], firstStrandGroup))		// If the first record does not fit the required balTrans shape, skip record
        {
            continue;
        }

        if (firstStrandGroup.size() != 2)					// Simple balanced translocation mutations only involve two chromosomes
        {
            continue;
        }

	for (std::size_t j = i + 1; j < records.size(); ++j)			// Loop through all following records after the i-th data record,
        {									// looking for a second candidate for the balTrans
            std::map<int, std::vector<SDRfragment>> secondStrandGroup;		// Store the second group of two strands that contribute to balanced translocations.

            if (!getRecordStrandFragments(records[j], secondStrandGroup))	// If secondStrandGroup does not fit balTrans shape, skip record
            {
                continue;
            }

            if (secondStrandGroup.size() != 2)					// Simple balTrans mutations only involve two chromosomes
            {
                continue;
            }


	    // Store fragment original strand IDs for the first candidate
            std::vector<int> firstStrands;

	    for (auto& [strandID, fragments] : firstStrandGroup)
            {
                firstStrands.push_back(strandID);
            }

	    // Store the fragment original strand IDs for the second candidate
            std::vector<int> secondStrands;

            for (auto& [strandID, fragments] : secondStrandGroup)
            {
                secondStrands.push_back(strandID);
            }

	    // If the two original strands referenced by each candidate record are not equivalent, then this is not a balanced translocation
            if (firstStrands != secondStrands)
            {
                continue;
            }

	    // Once the strand IDs are confirmed to be identical, set strand A and strand B IDs
	    const int strandIDs[2] = {firstStrands[0], firstStrands[1]};	// equivalent to secondStrands

	    // Pools every fragment referencing each strand across BOTH records (this is what allows central-segment translocations, not only end swaps), remembering which
	    // record (0 = first, 1 = second) each fragment came from. The pooled fragments must tile each strand with every JUNCTION between the two records contiguous.
	    // Gaps INSIDE one record (or at a chromosome end) are independent deletions, not part of the translocation, and are reported by detectDeletions().
	    std::vector<double> breakpoints[2];
	    bool valid = true;

	    for (int s = 0; s < 2 && valid; ++s)
	    {
	        std::vector<SDRfragment> pooled;
	        std::vector<int> owners;

	        for (const SDRfragment& fragment : firstStrandGroup[strandIDs[s]])
	        {
	            pooled.push_back(fragment);
	            owners.push_back(0);
	        }

	        for (const SDRfragment& fragment : secondStrandGroup[strandIDs[s]])
	        {
	            pooled.push_back(fragment);
	            owners.push_back(1);
	        }

	        const double fullSize = fullSizeForStrand(masterHeader, strandIDs[s]);
	        std::vector<SDRnormalizedFragment> sorted;
	        std::vector<SDRflankGap> gaps;

	        if (fullSize <= 0.0 || !analyzeFlankingFragments(pooled, fullSize, balTraTolerance, sorted, gaps))
	        {
	            valid = false;
	            break;
	        }

	        // Which records (bit 0 = first, bit 1 = second) hold a fragment spanning exactly [low, high].
	        auto ownersOf = [&](double low, double high) -> int
	        {
	            int mask = 0;

	            for (std::size_t k = 0; k < pooled.size(); ++k)
	            {
	                if (approxEqual(std::min(pooled[k].oldStartPosition, pooled[k].oldEndPosition), low, balTraTolerance) && approxEqual(std::max(pooled[k].oldStartPosition, pooled[k].oldEndPosition), high, balTraTolerance))
	                {
	                    mask |= (1 << owners[k]);
	                }
	            }

	            return mask;
	        };

	        for (const SDRflankGap& gap : gaps)						// A gap at the junction (neighbours from the two different records) means material was lost there: deletion-translocation.
	        {
	            if (gap.leftIndex < 0 || gap.rightIndex < 0)
	            {
	                continue;
	            }

	            const int leftOwners = ownersOf(sorted[gap.leftIndex].low, sorted[gap.leftIndex].high);
	            const int rightOwners = ownersOf(sorted[gap.rightIndex].low, sorted[gap.rightIndex].high);

	            if ((leftOwners == 1 && rightOwners == 2) || (leftOwners == 2 && rightOwners == 1))
	            {
	                valid = false;
	            }
	        }

	        for (std::size_t k = 0; k + 1 < sorted.size(); ++k)				// Contiguous junctions between the two records are the translocation breakpoints.
	        {
	            if (!approxEqual(sorted[k].high, sorted[k + 1].low, balTraTolerance))
	            {
	                continue;
	            }

	            const int leftOwners = ownersOf(sorted[k].low, sorted[k].high);
	            const int rightOwners = ownersOf(sorted[k + 1].low, sorted[k + 1].high);

	            if ((leftOwners == 1 && rightOwners == 2) || (leftOwners == 2 && rightOwners == 1))
	            {
	                breakpoints[s].push_back(sorted[k].high);
	            }
	        }

	        valid = valid && !breakpoints[s].empty();					// Each strand needs at least one real break where the two records exchange material.
	    }

	    if (!valid)
	    {
	        continue;
	    }

	    // All balanced translocation mutation checks passed, store the strand and fragment information, as well as the break locations
	    SDRtranslocationEvent event{};
	    event.oldStrandAid = strandIDs[0];
	    event.oldStrandBid = strandIDs[1];
	    event.breakpointsA = breakpoints[0];
	    event.breakpointsB = breakpoints[1];
	    event.newStrandID1 = records[i].newStrandID;
	    event.newStrandID2 = records[j].newStrandID;

	    translocations.push_back(event);
	}
    }

    if (applyHierarchy)									// Chromothripsis / chromoplexy silence every balanced translocation on their strands.
    {
        const std::set<int> silenced = collectHigherOrderStrands(subHeader, numOriginalStrands);
        translocations.erase(std::remove_if(translocations.begin(), translocations.end(), [&](const SDRtranslocationEvent& e) { return silenced.count(e.oldStrandAid) > 0 || silenced.count(e.oldStrandBid) > 0; }), translocations.end());
    }

    return translocations;
}













// ----------------------------------------------------------------- //
// Detects EXTRACHROMOSOMAL DNA (ecDNA) events within a single cell. //
// ----------------------------------------------------------------- //

// Identical fragment signature to the long deletion events with the added indicator of data field 4 isLinear = 0 (for circular fragments).
// Generalized to N independent ecDNA mutations, and for ecDNA made up of multiple deleted fragments from the same original strand. Can have either
// One remaining strand with N deletions and N+1 fragments, followed by N records of ecDNA depicting the N ecDNA mutations, or followed by one or
// multiple records with multiple fragments forming ecDNA.
//
// Also covers ecDNA assembled from excised fragments of TWO OR MORE DIFFERENT original strands: each contributing strand still has its own
// separate flanking/remaining record (the same deletion-shaped gap(s) as the single-origin case), but some or all of their excised gap fragments
// are combined onto ONE shared circular new strand instead of each getting its own. A contributing strand may supply more than one of its own gap
// fragments to that shared strand (oldStrandIDs/remainingStrandIDs then repeat that strand once per matched segment, same as the same-origin
// case), alongside fragments from other strands.
inline std::vector<SDRecDNAevent> detectECDNA(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader, bool applyHierarchy = true)
{
    std::vector<SDRecDNAevent> ecDNAevents;					// Vector to store ecDNA event information (strand IDs, fragment lengths)

    // Group every fragment in this cell by the old strand ID it references, remembering which new-strand record it came from.
    // Unlike detectDeletions, the excised piece here is EXPECTED to be circular.
    std::map<int, std::vector<std::pair<const SDRdataRecord*, SDRfragment>>> groupOldStrand;


    // Tracks each newStrandID's TRUE total fragment count across its whole record, used below to reject a candidate whose fragments are only
    // PARTLY from the old strand being examined.
    std::map<int, std::size_t> newStrandIDtotalFragments;


    for (const SDRdataRecord& record : subHeader.dataRecords)			// Loop through each data record per cell subheader
    {
	if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))	// Check if a data record consists of a rearranged fragment, if not, skip
	{
	    continue;
	}

        newStrandIDtotalFragments[record.newStrandID] = record.fragments.size();


	for (const SDRfragment& fragment : record.fragments)			// Check the fragments the new strand is composed of
	{
	    groupOldStrand[fragment.oldStrandID].push_back({&record, fragment});// Store the fragments corresponding to each original (old) strand ID
	}
    }


    struct Gap								// Gap struct to store the start and end positions of multiple gaps caused by multiple ecDNA events in a single fragment
    {
        double startPos;
        double endPos;
    };

    // Populated as a side effect of the per-old-strand loop below, for every strand that has a valid remaining/flanking record, regardless of whether 
    // that strand also has any same-origin excised candidates of its own, and regardless of how many gaps it has. Used
    // afterward to match multi-origin circular records, whose fragments reference several different old strands at once (and may reference the SAME old 
    // strand more than once, if that strand has multiple gaps).
    std::map<int, std::pair<int, std::vector<Gap>>> gapsByOldStrand;		// oldStrandID -> {remainingNewStrandID, gaps}


    // This regroups the same three fragments, but now keyed by newStrandID instead of oldStrandID,
    // collapsing "which record" down to just an ID, and separately keeping a lookup back to the actual record object:
    for (auto& [oldStrandID, entries] : groupOldStrand)				// Loop through all the entries/fragments associated with a given old strand ID
    {
        std::map<int, std::vector<SDRfragment>> groupNewStrand;			// Vector to store the new strand ID and the fragments the new strand is composed of -> use to tell how many fragments each new strand ID is made of
        std::map<int, const SDRdataRecord*> recordByNewStrand;			// Need to get back from newStrandID the actual SDR data record to check if strand is linear or circular

        for (const auto& [record, fragment] : entries)
        {
            groupNewStrand[record->newStrandID].push_back(fragment);
            recordByNewStrand[record->newStrandID] = record;
        }


	// recordByNewStrand looks like this:

	// recordByNewStrand[newStrandID] = &recordnewStrandID

	// where &recordNewStrandID is the entire SDR data record to check if the strand is linear or circular


        // An ecDNA event touches at least two new-strand records: one linear record with the two or more flanking fragments, one
        // or more CIRCULAR records with the excised fragment(s).
        if (groupNewStrand.size() < 2)
        {
            continue;
        }

	bool strandHasMultiFragmentLinear = false;				// Decides whether a single-fragment linear record may act as the remaining record (terminal excision).
	for (const auto& [candidateID, candidateFragments] : groupNewStrand)
	{
	    if (candidateFragments.size() >= 2 && candidateFragments.size() == newStrandIDtotalFragments[candidateID] && recordByNewStrand[candidateID]->linear)
	    {
		strandHasMultiFragmentLinear = true;
	    }
	}

	int remainingStrandID = -1;						// newStrandID for the strand missing an excised fragment
        std::vector<SDRfragment>* flankingFragments = nullptr;			// Vector containing the pointer to the two fragments missing a central excised segment
        bool multipleRemainingRecords = false;
        std::vector<std::pair<int, std::vector<SDRfragment>*>> excisedCandidates; // {newStrandID, fragments}


        for (auto& [newStrandID, fragments] : groupNewStrand)
        {

	    if (fragments.size() != newStrandIDtotalFragments[newStrandID]) // Reject if this record has fragments from OTHER old strands too.
            {
                continue;
            }

            const SDRdataRecord* record = recordByNewStrand[newStrandID];	// Assign the SDR data record to a given newStrandID

            if (record->linear && (fragments.size() >= 2 || (!strandHasMultiFragmentLinear && fragments.size() == 1 && fragments[0].hasCentromere && fragments[0].oldStartPosition <= fragments[0].oldEndPosition)))	// Linear remaining record: 2+ flanking fragments, or (only if the strand has none) a single centromere-bearing fragment left by a terminal excision
            {
		if (flankingFragments != nullptr)
                {
                    multipleRemainingRecords = true; 				// More than one candidate "remaining" record - invalid shape.
                    break;
                }

                remainingStrandID = newStrandID;
                flankingFragments = &fragments;
            }
            else if (!record->linear)						// Necessary format is one entry that is non-linear with one or more fragments, all referencing the same oldStrandID
            {
		// Every fragment in an ecDNA (excised) record must be acentric.
                bool hasCentromereFlag = false;

		for (const SDRfragment& fragment : fragments)
                {
                    if (fragment.hasCentromere)
                    {
                        hasCentromereFlag = true;
                        break;
                    }
                }

                if (hasCentromereFlag)
                {
                    continue; // Disqualified - not added as an excised candidate.
                }

                excisedCandidates.push_back({newStrandID, &fragments});
            }
        }


	if (multipleRemainingRecords || flankingFragments == nullptr)		// Without a valid remaining/flanking record, there's nothing to compute a gap from at all.
        {
            continue;
        }

	const double ecDNAtolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

	// Gaps (leading, internal and trailing) from NORMALIZED fragments, so an inverted fragment next to the gap doesn't invalidate the shape (an ecDNA gap touching
	// an inversion is still an ecDNA formation).
        const std::size_t sizeIndex = static_cast<std::size_t>(oldStrandID);
        const double chromosomeFullSizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size()) ? masterHeader.intactChromosomeSizes[sizeIndex] : 0.0;

        std::vector<Gap> gaps;
        std::vector<std::size_t> gapOrigin;					// Index of the deletion (gap on the remaining record) each entry of gaps was cut from.
        {
            std::vector<SDRnormalizedFragment> sortedFlanking;
            std::vector<SDRflankGap> flankGaps;

            if (analyzeFlankingFragments(*flankingFragments, chromosomeFullSizeMbp, ecDNAtolerance, sortedFlanking, flankGaps))
            {
                for (const SDRflankGap& flankGap : flankGaps)
                {
                    // A gap split into several circular pieces (e.g. a centric ring next to an acentric ecDNA) is matched piece by piece,
                    // so each acentric piece can claim its own part. Centric pieces are never candidates, so their part stays unclaimed.
                    const std::vector<std::pair<double, double>> pieces = circularPiecesTilingGap(subHeader, numOriginalStrands, oldStrandID, flankGap.start, flankGap.end, ecDNAtolerance);

                    const std::size_t origin = &flankGap - flankGaps.data();

                    if (pieces.empty())
                    {
                        gaps.push_back({flankGap.start, flankGap.end});
                        gapOrigin.push_back(origin);
                    }
                    else
                    {
                        for (const std::pair<double, double>& piece : pieces)
                        {
                            gaps.push_back({piece.first, piece.second});
                            gapOrigin.push_back(origin);
                        }
                    }
                }
            }
        }


	// Remember this strand's gap shape, however many gaps it has, for the multi-origin matching pass below, regardless of whether it also has 
	// any same-origin excised candidates.
	if (!gaps.empty())
	{
	    gapsByOldStrand[oldStrandID] = { remainingStrandID, gaps };
	}


	if (excisedCandidates.empty())						// Nothing of this strand's own origin to match against its gaps, the strand may still be
        {									// picked up by the multi-origin pass below via gapsByOldStrand.
            continue;
        }

	// Each circular candidate is matched on its own against the strand's still-unclaimed gaps. A gap with no circular partner (a lone gap, or one with a linear
	// excised record) is simply not ecDNA and is left to the deletion detectors, so gaps and circular pieces no longer need to match up one-to-one overall.
	std::vector<bool> gapClaimed(gaps.size(), false);

	// One deletion whose missing segment shattered into several acentric circles is ONE ecDNA formation (one mutational event), so circles cut from the
	// same gap are merged into a single event. A centric ring cut from that gap is still scored separately, as a ring chromosome.
	const std::size_t strandEventsStart = ecDNAevents.size();
	std::vector<std::set<std::size_t>> eventOrigins;			// Gaps each of this strand's events came from, index-aligned with ecDNAevents from strandEventsStart.

	for (const auto& [newStrandID, fragments] : excisedCandidates)
        {
            std::vector<bool> trialClaimed = gapClaimed;			// Rolled back if this candidate fails to match completely.
            std::vector<std::pair<double, double>> matchedSegments;
            std::set<std::size_t> matchedOrigins;
            bool allMatched = true;

            for (const SDRfragment& fragment : *fragments)
            {
                const double fragLow = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
                const double fragHigh = std::max(fragment.oldStartPosition, fragment.oldEndPosition);
                bool matched = false;

                for (std::size_t i = 0; i < gaps.size(); ++i)
                {
                    if (trialClaimed[i])
                    {
                        continue;
                    }

                    if (approxEqual(fragLow, gaps[i].startPos, ecDNAtolerance) && approxEqual(fragHigh, gaps[i].endPos, ecDNAtolerance))
                    {
                        trialClaimed[i] = true;
                        matched = true;
                        matchedSegments.push_back({gaps[i].startPos, gaps[i].endPos});
                        matchedOrigins.insert(gapOrigin[i]);
                        break;
                    }
                }

                if (!matched)
                {
                    allMatched = false;
                    break;
                }
            }

            if (!allMatched || matchedSegments.empty())
            {
                continue;
            }

            gapClaimed = trialClaimed;

            std::size_t sameDeletion = eventOrigins.size();			// An earlier event of this strand sharing a gap with this circle, if any.
            for (std::size_t k = 0; k < eventOrigins.size() && sameDeletion == eventOrigins.size(); ++k)
            {
                for (std::size_t origin : matchedOrigins)
                {
                    if (eventOrigins[k].count(origin))
                    {
                        sameDeletion = k;
                        break;
                    }
                }
            }

            if (sameDeletion < eventOrigins.size())				// Another circle of the same shattered segment, extend that event (keeps the first circle's excisedStrandID).
            {
                SDRecDNAevent& event = ecDNAevents[strandEventsStart + sameDeletion];
                event.oldStrandIDs.insert(event.oldStrandIDs.end(), matchedSegments.size(), oldStrandID);
                event.ecDNAsegments.insert(event.ecDNAsegments.end(), matchedSegments.begin(), matchedSegments.end());
                event.remainingStrandIDs.insert(event.remainingStrandIDs.end(), matchedSegments.size(), remainingStrandID);
                eventOrigins[sameDeletion].insert(matchedOrigins.begin(), matchedOrigins.end());
                continue;
            }

            SDRecDNAevent event{};
            event.oldStrandIDs.assign(matchedSegments.size(), oldStrandID);
            event.ecDNAsegments = matchedSegments;
            event.remainingStrandIDs.assign(matchedSegments.size(), remainingStrandID);
            event.excisedStrandID = newStrandID;

            ecDNAevents.push_back(event);
            eventOrigins.push_back(matchedOrigins);
        }

    }

    // ---------------------------------------------------------------------------- //
    // -------------- Multi-foreign-fragment ecDNA detection section -------------- //
    // ---------------------------------------------------------------------------- //

    // Circular records combining fragments from TWO OR MORE different original strands. Each contributing strand may supply ONE OR MORE of its own fragments (matched
    // against that strand's own gap set computed above), not just a single fragment each. So a strand with several gaps can have more than one of them recombine into the same circular
    // piece, alongside fragments from other strands.
    //
    // A record can never be picked up by both this section and the same-origin loop above: the same-origin excisedCandidates check above requires 
    // fragments.size() == newStrandIDtotalFragments[newStrandID], i.e. EVERY fragment on the record belongs to the one oldStrandID being grouped.
    // A record spanning 2+ distinct strands fails that check under every oldStrandID grouping it touches, so it never ends up in excisedCandidates 
    // and only ever gets evaluated here.
    for (const SDRdataRecord& record : subHeader.dataRecords)			// Every record in the cell is a candidate circular ecDNA segment, checked independently of the per-strand loop above.
    {
        if (record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands) || record.fragments.size() < 2)	// Must be circular (this IS the excised piece, not a flanking/remaining record), a genuine rearrangement, and have 2+ fragments to even possibly span multiple strands.
        {
            continue;
        }

        std::set<int> distinctStrands;						// Every distinct original strand this circular record draws a fragment from.

        for (const SDRfragment& fragment : record.fragments)
        {
            distinctStrands.insert(fragment.oldStrandID);
        }

        if (distinctStrands.size() < 2)
        {
            continue; 								// Only one contributing strand, already handled above for single-origin fragment ecDNA.
        }

        // Track, per contributing strand, which of that strand's own gaps this record has already claimed, so two fragments from the same 
        // strand can't both match the same gap.
        std::map<int, std::vector<bool>> claimedByStrand;

        for (int strandID : distinctStrands)
        {
            const auto strandIt = gapsByOldStrand.find(strandID);
	    // No known gap shape for this strand at all yields an EMPTY claim vector, so any fragment from it will fail to match below rather than crash.
            claimedByStrand[strandID] = std::vector<bool>(strandIt != gapsByOldStrand.end() ? strandIt->second.second.size() : 0, false);
        }

        SDRecDNAevent event{};							// Built up fragment-by-fragment below; oldStrandIDs/ecDNAsegments/remainingStrandIDs stay index-aligned with each other.
        event.excisedStrandID = record.newStrandID;
        bool allMatched = true;							// Every fragment on this record must match a real gap, or the whole record is rejected.

        for (const SDRfragment& fragment : record.fragments)			// Keep the record's own fragment order in the assembled event.
        {
            if (fragment.hasCentromere)						// Every fragment in an ecDNA circular record must be acentric, a centromere here means this isn't a genuine excised gap fragment.
            {
                allMatched = false;
                break;
            }

            const auto strandIt = gapsByOldStrand.find(fragment.oldStrandID);

            if (strandIt == gapsByOldStrand.end())				// This fragment's strand has no known gap shape at all (no valid remaining/flanking record was found for it above), so there's nothing to verify it against.
            {
                allMatched = false;
                break;
            }

            const int remainingStrandIDForFragment = strandIt->second.first;	// newStrandID of this fragment's own strand's remaining/flanking record.
            const std::vector<Gap>& strandGaps = strandIt->second.second;	// Every gap previously computed for this fragment's own strand.
            std::vector<bool>& claimed = claimedByStrand[fragment.oldStrandID];

            bool matched = false;

            for (std::size_t i = 0; i < strandGaps.size(); ++i)			// Look for one of this fragment's OWN strand's gaps that matches its position and hasn't been claimed yet.
            {
                if (claimed[i])							// Already filled by an earlier fragment on this same record, don't double-claim it.
                {
                    continue;
                }

                if (approxEqual(fragment.oldStartPosition, strandGaps[i].startPos, 0.001) && approxEqual(fragment.oldEndPosition, strandGaps[i].endPos, 0.001))
                {
                    claimed[i] = true;
                    matched = true;
                    event.oldStrandIDs.push_back(fragment.oldStrandID);
                    event.ecDNAsegments.push_back({strandGaps[i].startPos, strandGaps[i].endPos});
                    event.remainingStrandIDs.push_back(remainingStrandIDForFragment);
                    break;
                }
            }

            if (!matched)							// This fragment doesn't correspond to any of its own strand's known gaps, reject the whole record.
            {
                allMatched = false;
                break;
            }
        }

        if (!allMatched)							// Some fragment failed to match, this record isn't a valid multi-origin ecDNA shape.
        {
            continue;
        }

        ecDNAevents.push_back(event);						// Every fragment matched a distinct gap of its own strand, this is a confirmed multi-origin ecDNA event.
    }


    if (applyHierarchy)								// ecDNA pieces on a strand of a chromothripsis or chromoplexy cluster are scored as that event, not also as ecDNA.
    {
        const std::set<int> higherOrderStrands = collectHigherOrderStrands(subHeader, numOriginalStrands);	// Chromothripsis and chromoplexy strands

        ecDNAevents.erase(std::remove_if(ecDNAevents.begin(), ecDNAevents.end(), [&](const SDRecDNAevent& e)
        {
            for (int strandID : e.oldStrandIDs)
            {
                if (higherOrderStrands.count(strandID))
                {
                    return true;
                }
            }
            return false;
        }), ecDNAevents.end());
    }

    return ecDNAevents;
}














// --------------------------------------------------------------- //
// Function to detect RING CHROMOSOMES in a single cell            //
// --------------------------------------------------------------- //

// A ring chromosome is any circular (non-linear) rearranged record that STILL carries a centromere, the mirror image of ecDNA, whose defining
// trait is the opposite (circular, but with NO centromere at all, see detectECDNA()). A ring chromosome's fragments may
// originate from multiple different original strands (a fused, multi-chromosome ring), so no restriction is placed on how many
// distinct old strand IDs contribute to it.
inline std::vector<SDRchromosomeRingEvent> detectChromosomeRings(const SDRsubHeader& subHeader, int numOriginalStrands)
{
    std::vector<SDRchromosomeRingEvent> rings;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands))	// Must be circular and a genuine rearrangement, not an intact/baseline strand.
        {
            continue;
        }

        bool hasCentromere = false;
        int centromereOldStrandID = -1;
        std::vector<int> involvedOldStrandIDs;

        for (const SDRfragment& fragment : record.fragments)
        {
            bool alreadySeen = false;

            for (int strandID : involvedOldStrandIDs)
            {
                if (strandID == fragment.oldStrandID)
                {
                    alreadySeen = true;
                    break;
                }
            }

            if (!alreadySeen)
            {
                involvedOldStrandIDs.push_back(fragment.oldStrandID);
            }

            if (fragment.hasCentromere && !hasCentromere)		// Record the FIRST centromere-bearing fragment's original strand.
            {
                hasCentromere = true;
                centromereOldStrandID = fragment.oldStrandID;
            }
        }

        if (!hasCentromere)						// No centromere present - this is an ecDNA-shaped fragment instead, not a ring chromosome.
        {
            continue;
        }

        SDRchromosomeRingEvent event{};
        event.newStrandID = record.newStrandID;
        event.involvedOldStrandIDs = involvedOldStrandIDs;
        event.centromereOldStrandID = centromereOldStrandID;

        rings.push_back(event);
    }

    return rings;
}




















// ---------------------------------------------------------------- //
// Function to detect DELETION-INVERSION mutations in a single cell //
// ---------------------------------------------------------------- //

// One new-strand record contains three fragments from the same old strand: two in normal orientation and one reversed. A second
// new-strand record contains a single fragment from the same old strand = the deleted segment. Unlike a plain balanced inversion,
// the three fragments are NOT fully contiguous: there is exactly one gap among them (where the deletion occurred), and that gap must
// exactly match the second record's excised fragment (within tolerance).
inline std::vector<SDRdeletionInversionEvent> detectDeletionInversions(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader, bool applyHierarchy = true)
{
    const double delInvTolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

    std::vector<SDRdeletionInversionEvent> delInvs;

    for (const auto& [oldStrandID, records] : collectFlankingRecords(subHeader, numOriginalStrands))
    {
        if (records.size() != 1)
        {
            continue;
        }

        const SDRdataRecord& flanking = *records[0];

        std::vector<SDRnormalizedFragment> sorted;
        std::vector<SDRflankGap> gaps;

        // A deletion-inversion listed from the other chromosome end (outer pieces reversed and swapped, middle forward) is read the normal way first.
        if (!analyzeFlankingFragments(fragmentsInReadingOrder(flanking, masterHeader), fullSizeForStrand(masterHeader, oldStrandID), delInvTolerance, sorted, gaps))
        {
            continue;
        }

        for (int index = 0; index < static_cast<int>(sorted.size()); ++index)
        {
            if (!sorted[index].reversed)
            {
                continue;
            }

            std::vector<const SDRflankGap*> touching;				// Gaps bordering this inverted fragment (terminal ones included)
            for (const SDRflankGap& gap : gaps)
            {
                if (gap.leftIndex == index || gap.rightIndex == index)
                {
                    touching.push_back(&gap);
                }
            }

            if (touching.size() != 1)						// Gaps on both sides: skipped for now. No gap: balanced inversion, handled elsewhere.
            {
                continue;
            }

            const SDRflankGap& gap = *touching[0];
            const int otherIndex = (gap.leftIndex == index) ? gap.rightIndex : gap.leftIndex;

            if (otherIndex >= 0 && sorted[otherIndex].reversed)			// Two inversions meeting the gap: not a simple deletion-inversion.
            {
                continue;
            }

            int excisedStrandID = -1;
            const SDRgapPartner partner = classifyGapPartner(subHeader, numOriginalStrands, oldStrandID, gap.start, gap.end, delInvTolerance, {flanking.newStrandID}, excisedStrandID, true);	// The loose excised piece may be listed inverted

            if (partner == SDRgapPartner::Circular || partner == SDRgapPartner::Foreign)	// ecDNA formation (do not count the inversion) / more complex.
            {
                continue;
            }

            SDRdeletionInversionEvent event{};
            event.oldStrandID = oldStrandID;
            event.inversionStart = sorted[index].low;
            event.inversionEnd = sorted[index].high;
            event.deletionStart = gap.start;
            event.deletionEnd = gap.end;
            event.remainingStrandID = flanking.newStrandID;
            event.excisedStrandID = excisedStrandID;

            delInvs.push_back(event);
        }
    }

    if (applyHierarchy)
    {
        const std::set<int> silenced = collectHigherOrderStrands(subHeader, numOriginalStrands);
        delInvs.erase(std::remove_if(delInvs.begin(), delInvs.end(), [&](const SDRdeletionInversionEvent& e) { return silenced.count(e.oldStrandID) > 0; }), delInvs.end());
    }

    return delInvs;
}











// Checks whether two fragments referencing the same old strand ID have a genuine gap between them (rather than meeting cleanly at a
// breakpoint), i.e. material between them is missing. If so, populates gapStart/gapEnd with the gap's bounds.
inline bool fragmentsHaveGap(const SDRfragment& fragment1, const SDRfragment& fragment2, double& gapStart, double& gapEnd)
{
    const double delTolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

    if (fragment1.oldStrandID != fragment2.oldStrandID)		// Both fragments being compared must have the same old strand ID
    {
        return false;
    }

    if (isReversedFragment(fragment1) || isReversedFragment(fragment2))	// Neither fragment being compared should be inverted
    {
        return false;
    }

    // If fragment 1 end is less than fragment 2 start and the positions are not approximately equal, then a genuine gap exists and compute the gap
    if (fragment1.oldEndPosition < fragment2.oldStartPosition && !approxEqual(fragment1.oldEndPosition, fragment2.oldStartPosition, delTolerance))
    {
        gapStart = fragment1.oldEndPosition;
        gapEnd = fragment2.oldStartPosition;
        return true;
    }

    // If fragment 2 end is less than fragment 1 start and the positions are not approximately equal, then a genuine gap exists and compute the gap
    if (fragment2.oldEndPosition < fragment1.oldStartPosition && !approxEqual(fragment2.oldEndPosition, fragment1.oldStartPosition, delTolerance))
    {
        gapStart = fragment2.oldEndPosition;
        gapEnd = fragment1.oldStartPosition;
        return true;
    }

    // Otherwise, no gap
    return false;
}











// --------------------------------------------------------------- //
// Function to detect DELETION-TRANSLOCATIONS in a single cell     //
// --------------------------------------------------------------- //

// Two candidate records, each with exactly 2 fragments from 2 different old strands (linear, any orientation, a reversed fragment represents an inverted exchange, e.g.
// one forming a dicentric/acentric pair, and is not itself a reason to reject). For each of the two strands involved, the two records' fragments for that strand are compared
// using NORMALIZED [min,max] bounds (so a reversed fragment's boundaries are still checked correctly, unlike a raw oldStart/oldEnd comparison): if they touch exactly,
// that strand is "clean" (no deletion); if BOTH strands are clean, this is a plain balanced translocation, handled elsewhere. If exactly one strand has a genuine gap
// instead, a third record (single fragment, from the gapped strand, non-reversed, exactly filling the gap) confirms a deletion-translocation.

inline std::vector<SDRdeletionTranslocationEvent> detectDeletionTranslocations(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader, bool applyHierarchy)
{
    std::vector<SDRdeletionTranslocationEvent> delTras;					// Store deletion-translocation events

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;
    const double tolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

    // Single-use function to determine which fragments in a record should be considered for deletion-translocation.
    auto getRecordStrandFragments = [&](const SDRdataRecord& record, std::map<int, std::vector<SDRfragment>>& byStrand) -> bool
    {
        if (!record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands))
        {
            return false;
        }

        for (const SDRfragment& fragment : record.fragments)
        {
            byStrand[fragment.oldStrandID].push_back(fragment);
        }

        return !byStrand.empty();
    };

    // A fragment of the pooled strand together with which of the two translocation records (0 = first, 1 = second) it came from.
    struct PooledFragment
    {
        SDRfragment fragment;
        int owner;
    };

    // Finds every gap of ONE strand's pooled fragments (both records together) against the strand's full length. Any number of gaps is allowed: the one
    // sitting at the translocation junction is the deletion-translocation, every other gap is an independent deletion (reported by detectDeletions()).
    // Returns false for an invalid shape (unknown size, overlapping fragments).
    auto analyzeStrand = [&](const std::vector<PooledFragment>& pooled, int strandID, std::vector<SDRnormalizedFragment>& sorted, std::vector<SDRflankGap>& gaps) -> bool
    {
        const double fullSize = fullSizeForStrand(masterHeader, strandID);

        if (fullSize <= 0.0)
        {
            return false;
        }

        std::vector<SDRfragment> fragments;
        for (const PooledFragment& entry : pooled)
        {
            fragments.push_back(entry.fragment);
        }

        return analyzeFlankingFragments(fragments, fullSize, tolerance, sorted, gaps);
    };

    // Which translocation records (bit 0 = first, bit 1 = second) hold a pooled fragment spanning exactly [low, high].
    auto ownersOf = [&](const std::vector<PooledFragment>& pooled, double low, double high) -> int
    {
        int mask = 0;

        for (const PooledFragment& entry : pooled)
        {
            if (approxEqual(std::min(entry.fragment.oldStartPosition, entry.fragment.oldEndPosition), low, tolerance) && approxEqual(std::max(entry.fragment.oldStartPosition, entry.fragment.oldEndPosition), high, tolerance))
            {
                mask |= (1 << entry.owner);
            }
        }

        return mask;
    };


    for (std::size_t i = 0; i < records.size(); ++i)					// Loop through all SDR data records
    {
        std::map<int, std::vector<SDRfragment>> firstStrandGroup;

        // The record must be a linear rearrangement made of fragments from exactly two different original strands.
        if (!getRecordStrandFragments(records[i], firstStrandGroup) || firstStrandGroup.size() != 2)
        {
            continue;
        }

        for (std::size_t j = i + 1; j < records.size(); ++j)
        {
            std::map<int, std::vector<SDRfragment>> secondStrandGroup;

            if (!getRecordStrandFragments(records[j], secondStrandGroup) || secondStrandGroup.size() != 2)
            {
                continue;
            }

            std::vector<int> firstStrands;
            for (auto& [strandID, fragments] : firstStrandGroup)
            {
                firstStrands.push_back(strandID);
            }

            std::vector<int> secondStrands;
            for (auto& [strandID, fragments] : secondStrandGroup)
            {
                secondStrands.push_back(strandID);
            }

            if (firstStrands != secondStrands)						// Both records must draw on the same two strands.
            {
                continue;
            }

            const int strandIDs[2] = {firstStrands[0], firstStrands[1]};

            std::vector<PooledFragment> pooled[2];
            std::vector<SDRnormalizedFragment> sorted[2];
            std::vector<SDRflankGap> gaps[2];
            bool valid = true;

            for (int s = 0; s < 2; ++s)
            {
                for (const SDRfragment& fragment : firstStrandGroup[strandIDs[s]])
                {
                    pooled[s].push_back({fragment, 0});
                }

                for (const SDRfragment& fragment : secondStrandGroup[strandIDs[s]])
                {
                    pooled[s].push_back({fragment, 1});
                }

                valid = valid && analyzeStrand(pooled[s], strandIDs[s], sorted[s], gaps[s]);
            }

            if (!valid)
            {
                continue;
            }

            if (gaps[0].empty() && gaps[1].empty())					// No gap anywhere: plain balanced translocation, handled elsewhere.
            {
                continue;
            }

            // The deletion-translocation gap sits AT the junction: the fragment ending at the gap start and the fragment starting at the gap end come from the two
            // DIFFERENT translocation records. Terminal gaps have no neighbour on one side, so they never qualify. Gaps inside one record are independent deletions.
            int junctionStrandIndex = -1;
            const SDRflankGap* junctionGap = nullptr;
            int junctionCount = 0;

            for (int s = 0; s < 2; ++s)
            {
                for (const SDRflankGap& gap : gaps[s])
                {
                    if (gap.leftIndex < 0 || gap.rightIndex < 0)
                    {
                        continue;
                    }

                    const int leftOwners = ownersOf(pooled[s], sorted[s][gap.leftIndex].low, sorted[s][gap.leftIndex].high);
                    const int rightOwners = ownersOf(pooled[s], sorted[s][gap.rightIndex].low, sorted[s][gap.rightIndex].high);

                    if ((leftOwners == 1 && rightOwners == 2) || (leftOwners == 2 && rightOwners == 1))
                    {
                        ++junctionCount;
                        junctionStrandIndex = s;
                        junctionGap = &gap;
                    }
                }
            }

            if (junctionCount != 1)							// None: just a translocation with lone deletions elsewhere. Several: too complex for this shape.
            {
                continue;
            }

            const int deletedOldStrandID = strandIDs[junctionStrandIndex];
            const double gapStart = junctionGap->start;
            const double gapEnd = junctionGap->end;

            // Breakpoints of the OTHER strand: where its two records' fragments meet.
            const int otherIndex = 1 - junctionStrandIndex;
            std::vector<double> cleanBreakPositions;

            for (std::size_t k = 0; k + 1 < sorted[otherIndex].size(); ++k)
            {
                if (!approxEqual(sorted[otherIndex][k].high, sorted[otherIndex][k + 1].low, tolerance))
                {
                    continue;
                }

                const int leftOwners = ownersOf(pooled[otherIndex], sorted[otherIndex][k].low, sorted[otherIndex][k].high);
                const int rightOwners = ownersOf(pooled[otherIndex], sorted[otherIndex][k + 1].low, sorted[otherIndex][k + 1].high);

                if ((leftOwners == 1 && rightOwners == 2) || (leftOwners == 2 && rightOwners == 1))
                {
                    cleanBreakPositions.push_back(sorted[otherIndex][k].high);
                }
            }

            int excisedStrandID = -1;
            const SDRgapPartner partner = classifyGapPartner(subHeader, numOriginalStrands, deletedOldStrandID, gapStart, gapEnd, tolerance, {records[i].newStrandID, records[j].newStrandID}, excisedStrandID, true);	// The loose excised piece may be listed inverted

            if (partner == SDRgapPartner::Circular || partner == SDRgapPartner::Foreign)	// ecDNA / more complex: not a deletion-translocation
            {
                continue;
            }

            SDRdeletionTranslocationEvent event{};
            event.oldStrandAid = strandIDs[0];
            event.oldStrandBid = strandIDs[1];
            event.deletedOldStrandID = deletedOldStrandID;
            event.cleanBreakPositions = cleanBreakPositions;
            event.deletionStart = gapStart;
            event.deletionEnd = gapEnd;
            event.newStrandID1 = records[i].newStrandID;
            event.newStrandID2 = records[j].newStrandID;
            event.excisedStrandID = excisedStrandID;					// -1 when there is no excised record

            delTras.push_back(event);
        }
    }

    if (applyHierarchy)									// Chromothripsis / chromoplexy silence every deletion-translocation on their strands.
    {
        const std::set<int> silenced = collectHigherOrderStrands(subHeader, numOriginalStrands);
        delTras.erase(std::remove_if(delTras.begin(), delTras.end(), [&](const SDRdeletionTranslocationEvent& e) { return silenced.count(e.oldStrandAid) > 0 || silenced.count(e.oldStrandBid) > 0; }), delTras.end());
    }

    return delTras;
}















// ------------------------------------------------------------------------------------ //
// Function to detect DELETION-INSERTION entries in SDR file
// ------------------------------------------------------------------------------------ //

// A deletion-insertion event is characterized by two new-strand records:
// 1) A "donor" record: two fragments from the same old strand ID, with a real gap between them (identical shape to a plain long deletion's "remaining" piece).
// 2) A "recipient" record: three fragments - two from a DIFFERENT old strand ID that are themselves contiguous (the recipient's own material, split by the insertion),
// and one FOREIGN fragment (from the donor's old strand ID) sitting between them, whose position exactly matches the donor's gap.
// Unlike a balanced translocation, only ONE segment moves, the donor loses material with nothing coming back, and the recipient's own material is fully retained, just
// split by the inserted piece.
inline std::vector<SDRdeletionInsertionEvent> detectDeletionInsertions(const SDRsubHeader& subHeader, int numOriginalStrands, bool applyHierarchy)
{
    const double tolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

    std::vector<SDRdeletionInsertionEvent> delInsEvents;

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;

    for (const SDRdataRecord& donorRecord : records)				// Candidate DONOR: any pure, linear, rearranged record with 2+ fragments
    {
        if (!donorRecord.linear || !isRearrangementCandidate(donorRecord.newStrandID, numOriginalStrands) || donorRecord.fragments.size() < 2)
        {
            continue;
        }

        const int donorOldStrandID = donorRecord.fragments[0].oldStrandID;

        bool pure = true;
        for (const SDRfragment& fragment : donorRecord.fragments)
        {
            pure = pure && fragment.oldStrandID == donorOldStrandID;
        }

        if (!pure)
        {
            continue;
        }

        std::vector<SDRnormalizedFragment> sorted;
        std::vector<SDRflankGap> gaps;

        if (!analyzeFlankingFragments(donorRecord.fragments, 0.0, tolerance, sorted, gaps))	// Internal gaps only: the segment must have flanking material on both sides
        {
            continue;
        }

        for (const SDRflankGap& gap : gaps)
        {
            if (gapTouchesReversedFragment(sorted, gap))
            {
                continue;
            }

            const double segmentStart = gap.start;
            const double segmentEnd = gap.end;

            for (const SDRdataRecord& recipientRecord : records)		// Look for the matching RECIPIENT record
            {
                if (recipientRecord.newStrandID == donorRecord.newStrandID || !recipientRecord.linear || !isRearrangementCandidate(recipientRecord.newStrandID, numOriginalStrands) || recipientRecord.fragments.size() < 3)
                {
                    continue;
                }

                // The recipient holds ONE foreign fragment (from the donor strand), flanked on both sides by the recipient strand. The recipient's own
                // material may have gaps elsewhere (independent deletions), only the insertion site itself has to be a clean break. Other chromosomes may also
                // be attached to the recipient (an unrelated translocation on the same chromosome, e.g. 38/0/119.4, 30/0.8/2.2, 38/119.4/130.2, 18/11.4/132.6),
                // they are scored by their own detector.
                int foreignIndex = -1;
                bool valid = true;

                for (std::size_t k = 0; k < recipientRecord.fragments.size() && valid; ++k)
                {
                    if (recipientRecord.fragments[k].oldStrandID == donorOldStrandID)
                    {
                        valid = foreignIndex == -1;
                        foreignIndex = static_cast<int>(k);
                    }
                }

                if (!valid || foreignIndex <= 0 || foreignIndex + 1 >= static_cast<int>(recipientRecord.fragments.size()))
                {
                    continue;
                }

                const int recipientOldStrandID = recipientRecord.fragments[foreignIndex - 1].oldStrandID;

                if (recipientRecord.fragments[foreignIndex + 1].oldStrandID != recipientOldStrandID)	// Both neighbours must be the same recipient chromosome
                {
                    continue;
                }

                // The inserted piece may go in either orientation (e.g. 29/0/21.5, 20/27.4/25.0, 29/21.5/132.6 is chr 20's piece inserted inverted), so it is
                // compared by its low/high bounds. Every other fragment must keep its original orientation; a reversed one is a different event.
                bool otherReversed = false;
                for (std::size_t k = 0; k < recipientRecord.fragments.size(); ++k)
                {
                    otherReversed = otherReversed || (static_cast<int>(k) != foreignIndex && isReversedFragment(recipientRecord.fragments[k]));
                }

                if (otherReversed)
                {
                    continue;
                }

                const SDRfragment& foreignFragment = recipientRecord.fragments[foreignIndex];
                const double foreignLow = std::min(foreignFragment.oldStartPosition, foreignFragment.oldEndPosition);
                const double foreignHigh = std::max(foreignFragment.oldStartPosition, foreignFragment.oldEndPosition);

                if (!approxEqual(foreignLow, segmentStart, tolerance) || !approxEqual(foreignHigh, segmentEnd, tolerance))
                {
                    continue;
                }

                const SDRfragment& lowerNativeFrag = recipientRecord.fragments[foreignIndex - 1];		// Neighbours of the inserted piece, in record order
                const SDRfragment& upperNativeFrag = recipientRecord.fragments[foreignIndex + 1];

                if (!approxEqual(lowerNativeFrag.oldEndPosition, upperNativeFrag.oldStartPosition, tolerance))	// Recipient's own material must be cleanly split
                {
                    continue;
                }

                SDRdeletionInsertionEvent event{};
                event.donorOldStrandID = donorOldStrandID;
                event.recipientOldStrandID = recipientOldStrandID;
                event.segmentStart = segmentStart;
                event.segmentEnd = segmentEnd;
                event.insertionPoint = lowerNativeFrag.oldEndPosition;
                event.donorRemainingNewStrandID = donorRecord.newStrandID;
                event.recipientNewStrandID = recipientRecord.newStrandID;

                delInsEvents.push_back(event);
            }
        }
    }

    if (applyHierarchy)								// Chromothripsis / chromoplexy silence every deletion-insertion on their strands.
    {
        const std::set<int> silenced = collectHigherOrderStrands(subHeader, numOriginalStrands);
        delInsEvents.erase(std::remove_if(delInsEvents.begin(), delInsEvents.end(), [&](const SDRdeletionInsertionEvent& e) { return silenced.count(e.donorOldStrandID) > 0 || silenced.count(e.recipientOldStrandID) > 0; }), delInsEvents.end());
    }

    return delInsEvents;
}











// ------------------------------------------------------------------------------------ //
// Function to detect CHROMOPLEXY events in SDR file
// ------------------------------------------------------------------------------------ //

// Builds a graph where each old strand ID is a node. For every rearranged record, every PAIR of distinct old strand IDs referenced by its fragments gets an edge, regardless
// of whether that record matches any specific named mutation shape (translocation, deletion-translocation, etc.). Any connected component spanning 3 or more distinct
// strands is reported as one chromoplexy event.
// Only strands that RECEIVE foreign material are linked: in real chromoplexy every chromosome involved gets material from another one. A strand whose material only sits
// inside another chromosome (flanked by that chromosome's own fragments, e.g. 38/0/119.4, 30/0.8/2.2, 38/119.4/130.2) only donated it (a deletion-insertion), so it is
// left out of the graph and does not merge two unrelated events into one chromoplexy.
inline std::vector<SDRchromoplexyEvent> detectChromoplexy(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader, bool collectAccompanyingDeletions)
{
    std::vector<SDRchromoplexyEvent> chromoplexyEvents;

    const std::vector<SDRdeletionEvent> deletions = collectAccompanyingDeletions ? detectDeletions(subHeader, numOriginalStrands, masterHeader, false) : std::vector<SDRdeletionEvent>{};


    // Initial State: Every chromosome is isolated in its own group (parent[x] = x).
    // Mutation 1 (Chromosome 1 and 2 swap material): You call unite(1, 2). Chromosome 1 and 2 are now linked. parent[1] becomes 2.
    // Mutation 2 (Chromosome 2 and 3 swap material): You call unite(2, 3). Now, Chromosome 1, 2, and 3 are all part of the same group led by 3.
    // If you later ask findRoot(1), the code will climb up (1 -> 2 -> 3), optimize the path, and return 3.

    // Groups chromosomes based on their shared mutations
    std::map<int, int> parent;

    auto findRoot = [&](int x) -> int		// x are the individual chromosome IDs
    {
        if (!parent.count(x))			// If chromosome x was found in the parent collection or not, if not initialize x as the parent
        {
            parent[x] = x;
        }

        while (parent[x] != x)			// Find root element
        {
            parent[x] = parent[parent[x]];	// path compression
            x = parent[x];			// Move pointer to newly assigned parent
        }

        return x;				// Return the root representative of the set of chromosomes
    };


    // Merge two chromosome mutation groups into a single combined group
    auto unite = [&](int a, int b)
    {
        const int rootA = findRoot(a);
        const int rootB = findRoot(b);

        if (rootA != rootB)			// Check if a and b belong to different groups
        {
            parent[rootA] = rootB;		// Unite the separate groups if a and b are not joined already
        }
    };


    // A strand receives foreign material when it forms an end of a linear multi-strand record (its own arm with foreign material attached), or is part of a
    // circular multi-strand record (no ends, so no strand is just inserted into another). Strands only found in the middle of other chromosomes are donors only.
    std::set<int> receivingStrands;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands) || record.fragments.empty())
        {
            continue;
        }

        std::set<int> strandsInRecord;

        for (const SDRfragment& fragment : record.fragments)
        {
            strandsInRecord.insert(fragment.oldStrandID);
        }

        if (strandsInRecord.size() < 2)
        {
            continue;
        }

        if (record.linear)
        {
            receivingStrands.insert(record.fragments.front().oldStrandID);
            receivingStrands.insert(record.fragments.back().oldStrandID);
        }
        else
        {
            receivingStrands.insert(strandsInRecord.begin(), strandsInRecord.end());
        }
    }


    // For every rearranged record, track which distinct old strand IDs it references and union every pair of receiving strands found together.
    std::vector<std::pair<int, std::set<int>>> recordStrandSets; 			// {newStrandID, {distinct old strand IDs}}

    for (const SDRdataRecord& record : subHeader.dataRecords)				// Loop through all SDR data records per cell
    {
        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))		// If the record is not a mutation, but an intact record, skip
        {
            continue;
        }

	std::set<int> strandsInRecord;							// Track which chromosomes make up the fragments in the mutation

        for (const SDRfragment& fragment : record.fragments)
        {
            strandsInRecord.insert(fragment.oldStrandID);
        }

        if (strandsInRecord.size() < 2)							// Purely intra-strand record, contributes no edges (only one chromosome is involved).
        {
            continue;
        }

        recordStrandSets.push_back({record.newStrandID, strandsInRecord});		// Couple the new strand ID with the old strand IDs of its fragments

        std::vector<int> strandVec;							// The record's receiving strands; donor-only strands get no edges
        for (int strandID : strandsInRecord)
        {
            if (receivingStrands.count(strandID))
            {
                strandVec.push_back(strandID);
            }
        }

	// Loop through all original chromosome IDs, in the current record and following records
	for (std::size_t a = 0; a < strandVec.size(); ++a)
        {
            for (std::size_t b = a + 1; b < strandVec.size(); ++b)
            {
                unite(strandVec[a], strandVec[b]);					// Unite two old strand IDs from two records
            }
        }
    }


    // as you iterate all strands, every strand sharing the same root ends up collected together in one set, and using a set rather than a vector means 
    // duplicates (from path-compressed lookups) just don't matter.
    std::map<int, std::set<int>> componentsByRoot;					// maps one arbitrary representative strand the full set of every strand connected to it.

    // for every strand that's ever been registered as a node: call findRoot(strandID) to get its component's representative 
    // ("root") strand ID, and record that strandID belongs under that root.
    for (auto& [strandID, strandParent] : parent)
    {
        componentsByRoot[findRoot(strandID)].insert(strandID);
    }

    // Loops through every discovered component and checks if it has >= 3 strands
    for (auto& [root, strandSet] : componentsByRoot)
    {
        if (strandSet.size() < 3)							// Chromoplexy involved 3 or more chromosomes
        {
            continue;
        }

	// Build the list of records by newStrandID since the component qualifies with 3+ strands
        std::vector<int> contributingIDs;						// Will hold all records touching at least one of the componentsByRoot's 3+ strands.

	// The nested loop checks for each record if any old strand this record touches belong to the current 3+-strand component.
        for (auto& [newStrandID, strands] : recordStrandSets)
        {
	    // Check each strand reference in the given record
            for (int s : strands)
            {
                if (strandSet.count(s))							// If that strand is in strandSet, add newStrandID to the contributingIDs
                {
                    contributingIDs.push_back(newStrandID);				// Store the record ID and break, one match is proof the record belongs to the chromoplexy event, stop searching further.
                    break;
                }
            }
        }


	// All checks passed for chromoplexy, store the mutation event information
        SDRchromoplexyEvent event{};
        event.involvedStrandIDs = std::vector<int>(strandSet.begin(), strandSet.end());
        event.contributingNewStrandIDs = contributingIDs;


	// Extra non-essential check to assess whether one of the involved strands in chromoplexy also has a separate associated deletion, and stores that deletions information
        for (const SDRdeletionEvent& deletionEvent : deletions)
        {
            if (strandSet.count(deletionEvent.oldStrandID))
            {
                event.strandsWithAccompanyingDeletions.push_back(deletionEvent.oldStrandID);
            }
        }

        chromoplexyEvents.push_back(event);
    }

    return chromoplexyEvents;

}













// ------------------------------------------------------------------------------------ //
// Function to detect CHROMOTHRIPSIS events in SDR file
// ------------------------------------------------------------------------------------ //

// True when pieces of ONE chromosome, read in the order given, keep their original order along it: forward pieces follow each other up the chromosome, and a run of
// consecutive reversed pieces is one inverted block whose pieces follow each other down the chromosome. Gaps (deletions) are allowed. Simple deletions and inversions,
// however many, keep this order. Also tries the pieces read from the other end, since a chromosome can be listed from either end.
inline bool piecesInOriginalOrder(const std::vector<SDRfragment>& pieces)
{
    const double tolerance = svToleranceMbp;					// Fewer than 100 bases missing = touching (see svToleranceMbp)

    auto inOrder = [tolerance](const std::vector<SDRfragment>& fragments) -> bool
    {
        double previousHigh = -1.0;							// Highest position reached so far along the chromosome

        for (std::size_t i = 0; i < fragments.size(); )
        {
            if (!isReversedFragment(fragments[i]))
            {
                if (fragments[i].oldStartPosition < previousHigh - tolerance)
                {
                    return false;
                }

                previousHigh = fragments[i].oldEndPosition;
                ++i;
                continue;
            }

            std::size_t j = i + 1;							// Run of reversed pieces i..j-1 = one inverted block

            while (j < fragments.size() && isReversedFragment(fragments[j]))
            {
                // A reversed piece lying wholly above this block (e.g. 7/80/50 then 7/120/80) is a second inversion next to the first, each piece flipped
                // in place, not part of the same inverted block: it starts its own block, which the check below keeps in order.
                if (fragments[j].oldEndPosition >= fragments[i].oldStartPosition - tolerance)
                {
                    break;
                }

                if (fragments[j].oldStartPosition > fragments[j - 1].oldEndPosition + tolerance)	// Inside the block, each piece lies below the previous one
                {
                    return false;
                }

                ++j;
            }

            if (fragments[j - 1].oldEndPosition < previousHigh - tolerance)		// The whole block lies above what came before it
            {
                return false;
            }

            previousHigh = fragments[i].oldStartPosition;
            i = j;
        }

        return true;
    };

    std::vector<SDRfragment> flipped(pieces.rbegin(), pieces.rend());		// The same pieces read from the chromosome's other end

    for (SDRfragment& fragment : flipped)
    {
        std::swap(fragment.oldStartPosition, fragment.oldEndPosition);
    }

    return inOrder(pieces) || inOrder(flipped);
}


// True when a record of ONE chromosome with at least minFragments pieces has them out of their original order, e.g. 40/132.6/127.0, 40/115.2/127.0, 40/96.4/115.2,
// 40/96.4/0 (see piecesInOriginalOrder()).
inline bool isNonSequentialRearrangement(const SDRdataRecord& record, std::size_t minFragments = 4)
{
    if (!record.linear || record.fragments.size() < minFragments)
    {
        return false;
    }

    for (const SDRfragment& fragment : record.fragments)
    {
        if (fragment.oldStrandID != record.fragments.front().oldStrandID)	// Multi-chromosome records are left to the translocation-type detectors.
        {
            return false;
        }
    }

    return !piecesInOriginalOrder(record.fragments);
}


// True when a record joining two or more chromosomes carries exchanged material that was itself shattered: one chromosome's consecutive pieces (2 or more)
// are out of their original order, or the piece meeting the other chromosome is inverted on its own, e.g. 13/132.6/14.4, 13/1.9/14.4, 34/121.8/132.6: chr 13
// broke at 14.4 and the piece next to the junction flipped before joining chr 34. A plain translocation exchanges each chromosome's material
// in one piece (or one reversed piece), and an inversion with breaks inside the exchanged material keeps its pieces in order, so neither matches.
inline bool hasShatteredExchangedMaterial(const SDRdataRecord& record)
{
    if (!record.linear)
    {
        return false;
    }

    const std::vector<SDRfragment>& fragments = record.fragments;

    for (std::size_t runStart = 0; runStart < fragments.size(); )
    {
        std::size_t runEnd = runStart + 1;						// Run of consecutive pieces of one chromosome: runStart..runEnd-1

        while (runEnd < fragments.size() && fragments[runEnd].oldStrandID == fragments[runStart].oldStrandID)
        {
            ++runEnd;
        }

        const bool bordersOtherChromosome = runStart > 0 || runEnd < fragments.size();

        if (bordersOtherChromosome && runEnd - runStart >= 2)
        {
            const std::vector<SDRfragment> run(fragments.begin() + runStart, fragments.begin() + runEnd);

            // The piece meeting the other chromosome is flipped on its own: its orientation differs from its neighbour in the run (either way round, since a
            // chromosome can be listed from either end, e.g. 13/132.6/14.4 reversed then 13/1.9/14.4 forward next to chr 34).
            const bool flippedAtStartJunction = runStart > 0 && isReversedFragment(run.front()) != isReversedFragment(run[1]);
            const bool flippedAtEndJunction = runEnd < fragments.size() && isReversedFragment(run.back()) != isReversedFragment(run[run.size() - 2]);

            if (!piecesInOriginalOrder(run) || flippedAtStartJunction || flippedAtEndJunction)
            {
                return true;
            }
        }

        runStart = runEnd;
    }

    return false;
}


// Same underlying strand-connectivity graph as detectChromoplexy(), built from raw fragment data, so clustering doesn't depend on any named mutation shape. Components of
// size 3+ are skipped (chromoplexy's territory). For components of size 1 or 2, rather than counting named detected mutation events , this counts the TOTAL 
// number of data-field-3 fragments across every rearranged record touching the cluster's strand(s), a direct measure of how many
// pieces the strand(s) were actually broken into, regardless of shape. Only counts the remaining fragments in the original strand, not the excised pieces, as a
// clearer metric for the number of rearrangements that occurred in a strand.
// A cluster below the fragment threshold is a complex genomic rearrangement instead (see detectComplexGenomicRearrangements()) when one of its records
// holds 4 or more pieces of one chromosome rearranged out of their original order (see isNonSequentialRearrangement()), which simple deletions and
// inversions cannot produce, or joins two chromosomes with exchanged material that broke and was rearranged on the way (see
// hasShatteredExchangedMaterial()). Each cluster is scored as one or the other, never both.
// complexRearrangements selects which of the two kinds of cluster is returned.
inline std::vector<SDRchromothripsisEvent> detectShatteredClusters(const SDRsubHeader& subHeader, int numOriginalStrands, int minFragmentCount, bool complexRearrangements)
{
    std::vector<SDRchromothripsisEvent> chromothripsisEvents;			// Store chromothripsis (or complex genomic rearrangement) events

    // Groups chromosomes based on their shared mutations
    std::map<int, int> parent;

    auto findRoot = [&](int x) -> int		// x are the individual chromosome IDs
    {
        if (!parent.count(x))			// If chromosome x was found in the parent collection or not, if not initialize x as the parent
        {
            parent[x] = x;
        }

        while (parent[x] != x)			// Find root element
        {
            parent[x] = parent[parent[x]];	// path compression
            x = parent[x];			// Move pointer to newly assigned parent
        }

        return x;				// Return the root representative of the set of chromosomes
    };


    // Merge two chromosome mutation groups into a single combined group
    auto unite = [&](int a, int b)
    {
        const int rootA = findRoot(a);
        const int rootB = findRoot(b);

        if (rootA != rootB)			// Check if a and b belong to different groups
        {
            parent[rootA] = rootB;		// Unite the separate groups if a and b are not joined already
        }
    };


    // Map from old strand ID -> set of newStrandIDs (records) whose fragments touch it, whether intra-strand-only or inter-strand.
    std::map<int, std::set<int>> recordsByStrand;

    // Map from newStrandID -> the record's own total fragment count, so each record's fragments are only ever tallied once even if it
    // touches multiple strands in a cluster.
    std::map<int, std::size_t> fragmentCountByNewStrandID;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
	// If record is not a mutated record, skip
        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))
        {
            continue;
        }

	// Old strand IDs of all fragments in a given record
        std::set<int> strandsInRecord;

        for (const SDRfragment& fragment : record.fragments)
        {
            strandsInRecord.insert(fragment.oldStrandID);
        }

	// Count how many fragments there are across all involved new Strand IDs in a chromothripsis event (either 1 or 2 new Strand ID, i.e. data records)
        fragmentCountByNewStrandID[record.newStrandID] = record.fragments.size();

        for (int strandID : strandsInRecord)
        {
            findRoot(strandID); 				// Register the node even if this record touches only one strand.
            recordsByStrand[strandID].insert(record.newStrandID);
        }

        if (strandsInRecord.size() < 2)				// Intra-strand record, no edge to add.
        {
            continue;
        }

	// Store all the strand IDs involved across all chromothripsis-involved records
        std::vector<int> strandVec(strandsInRecord.begin(), strandsInRecord.end());

	// Unite a and b data records that could involve two old strand IDs in a chromothripsis event
        for (std::size_t a = 0; a < strandVec.size(); a++)
        {
            for (std::size_t b = a + 1; b < strandVec.size(); b++)
            {
                unite(strandVec[a], strandVec[b]);
            }
        }
    }

    // Maps one arbitrary representative strand to the full set of every strand connected to it.
    std::map<int, std::set<int>> componentsByRoot;

    for (auto& [strandID, strandParent] : parent)
    {
        componentsByRoot[findRoot(strandID)].insert(strandID);
    }


    for (auto& [root, strandSet] : componentsByRoot)
    {
        if (strandSet.size() > 2)
        {
            continue; // 3+ strands, chromoplexy's territory, not chromothripsis.
        }

        // Any record with 2+ fragments is a "major" surviving piece and counts in full, whether it's a single-strand remaining piece
        // (like a deletion's) or a recombined multi-strand derivative (like a translocation/chromothripsis product). A record with
        // exactly 1 fragment is always a trivial excised byproduct (true for every excised/circular shape our detectors define)
        // and gets excluded.
        std::set<int> majorNewStrandIDs;

        for (int strandID : strandSet)
        {
            for (int newStrandID : recordsByStrand[strandID])
            {
                if (fragmentCountByNewStrandID[newStrandID] >= 2)		// Ignores single excised fragments, only strands made up of multiple fragments count
                {								// towards chromothripsis
                    majorNewStrandIDs.insert(newStrandID);
                }
            }
        }

        int totalFragmentCount = 0;						// Count how many contributing fragments there are to be considered as chromothripsis

        for (int newStrandID : majorNewStrandIDs)				// Count all non-single-excised fragment pieces across the involved parent strands
        {
            totalFragmentCount += static_cast<int>(fragmentCountByNewStrandID[newStrandID]);
        }

        bool hasShatteredExchange = false;					// Exchanged material between two chromosomes broke and was rearranged
        bool hasNonSequentialRecord = false;					// 4+ pieces of one chromosome out of their original order

        for (const SDRdataRecord& record : subHeader.dataRecords)
        {
            if (!majorNewStrandIDs.count(record.newStrandID))
            {
                continue;
            }

            hasShatteredExchange = hasShatteredExchange || hasShatteredExchangedMaterial(record);
            hasNonSequentialRecord = hasNonSequentialRecord || isNonSequentialRearrangement(record);
        }

	// Min fragment count = 10 by default, function argument parameter can be changed. Only the fragment count makes chromothripsis.
        const bool isChromothripsis = totalFragmentCount >= minFragmentCount;

	// Safeguard: a cluster that is already chromothripsis is never also a complex genomic rearrangement, so the two counts never overlap.
        const bool isComplexRearrangement = !isChromothripsis && (hasNonSequentialRecord || hasShatteredExchange);

        if (complexRearrangements ? !isComplexRearrangement : !isChromothripsis)
        {
            continue;
        }

	// Checks passed, store the chromothripsis (or complex genomic rearrangement) event
        SDRchromothripsisEvent event{};
        event.involvedStrandIDs = std::vector<int>(strandSet.begin(), strandSet.end());
        event.contributingNewStrandIDs = std::vector<int>(majorNewStrandIDs.begin(), majorNewStrandIDs.end());
        event.totalFragmentCount = totalFragmentCount;

        chromothripsisEvents.push_back(event);
    }

    return chromothripsisEvents;
}


// Chromothripsis: a 1-2 strand cluster with 10+ fragments (see detectShatteredClusters()).
inline std::vector<SDRchromothripsisEvent> detectChromothripsis(const SDRsubHeader& subHeader, int numOriginalStrands, int minFragmentCount)
{
    return detectShatteredClusters(subHeader, numOriginalStrands, minFragmentCount, false);
}


// Complex genomic rearrangement: a 1-2 strand cluster below the chromothripsis fragment threshold with a record of 4+ pieces of one chromosome
// out of their original order, or with shattered exchanged material between chromosomes (see detectShatteredClusters()). It silences simpler events on its strands, like chromothripsis.
inline std::vector<SDRcomplexRearrangementEvent> detectComplexGenomicRearrangements(const SDRsubHeader& subHeader, int numOriginalStrands, int minFragmentCount)
{
    return detectShatteredClusters(subHeader, numOriginalStrands, minFragmentCount, true);
}


// Counts the detected mutations of one cell per original strand, for the web viewer's per-chromosome popups.
// Uses the same detection functions (and hierarchy) as the SDR summary file, so the numbers agree with it.
// A mutation involving several strands (e.g. a translocation between chromosomes 5 and 12) is counted once on each strand involved,
// so the per-strand counts can add up to more than the cell total. Only mutation types that occur are listed, in summary-file order.
// Exception: a dicentric is counted only on the chromosome slot it is drawn in, not on the partner left acentric.
// A chromosome with a change that none of the other mutation types accounts for gets one "Complex Genomic Rearrangements" entry, even if it has other
// mutations. Those chromosomes are also returned in unexplainedStrandsOut when it is given (see countComplexGenomicRearrangements()).
// Indels (< 100 bp) are counted once per site on the strand that received them. Every other type, including the CGR catch-all, is detected on the
// cell with its indels removed (see removeIndels()), so an indel never makes a chromosome count as a complex genomic rearrangement.
inline std::map<int, std::vector<std::pair<std::string, int>>> countMutationsPerOriginalStrand(const SDRsubHeader& rawSubHeader, const SDRmasterHeader& masterHeader,
                                                                                              std::set<int>* unexplainedStrandsOut = nullptr)
{
    const int numOriginalStrands = masterHeader.intactChromosomeSizes.empty() ? 0 : static_cast<int>(masterHeader.intactChromosomeSizes[0]);

    std::vector<SDRindelEvent> indels;
    const SDRsubHeader subHeader = removeIndels(rawSubHeader, numOriginalStrands, &indels);

    std::map<int, std::map<std::string, int>> counts;			// strand ID -> mutation type -> count

    // Count one mutation of this type on every distinct strand in strandIDs.
    auto addEvent = [&counts](const std::string& type, std::vector<int> strandIDs)
    {
        std::sort(strandIDs.begin(), strandIDs.end());
        strandIDs.erase(std::unique(strandIDs.begin(), strandIDs.end()), strandIDs.end());

        for (int strandID : strandIDs)
        {
            counts[strandID][type]++;
        }
    };

    const std::vector<SDRdeletionEvent> deletions = detectDeletions(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRinversionEvent> inversions = detectInversions(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRtranslocationEvent> translocations = detectTranslocations(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRecDNAevent> ecDNAs = detectECDNA(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRdeletionInversionEvent> delInvs = detectDeletionInversions(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRdeletionTranslocationEvent> delTras = detectDeletionTranslocations(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRdeletionInsertionEvent> delIns = detectDeletionInsertions(subHeader, numOriginalStrands);
    const std::vector<SDRchromoplexyEvent> chromoplexy = detectChromoplexy(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRchromothripsisEvent> chromothripsis = detectChromothripsis(subHeader, numOriginalStrands);
    const std::vector<SDRcomplexRearrangementEvent> complexRearrangements = detectComplexGenomicRearrangements(subHeader, numOriginalStrands);
    const std::vector<SDRchromosomeRingEvent> rings = detectChromosomeRings(subHeader, numOriginalStrands);
    const std::vector<SDRdicentricEvent> dicentrics = detectDicentrics(subHeader, numOriginalStrands);

    for (const auto& e : deletions)      addEvent("Long Deletions", {e.oldStrandID});
    for (const auto& e : inversions)     addEvent("Balanced Inversions", {e.oldStrandID});
    for (const auto& e : translocations) addEvent("Balanced Translocations", {e.oldStrandAid, e.oldStrandBid});
    for (const auto& e : ecDNAs)         addEvent("Extrachromosomal DNA (ecDNA)", e.oldStrandIDs);
    for (const auto& e : delInvs)        addEvent("Deletion-Inversions", {e.oldStrandID});
    for (const auto& e : delTras)        addEvent("Deletion-Translocations", {e.oldStrandAid, e.oldStrandBid});
    for (const auto& e : delIns)         addEvent("Deletion-Insertions", {e.donorOldStrandID, e.recipientOldStrandID});
    for (const auto& e : chromoplexy)    addEvent("Chromoplexy", e.involvedStrandIDs);
    for (const auto& e : chromothripsis) addEvent("Chromothripsis", e.involvedStrandIDs);
    for (const auto& e : complexRearrangements) addEvent("Complex Genomic Rearrangements", e.involvedStrandIDs);
    for (const auto& e : rings)          addEvent("Ring Chromosomes", e.involvedOldStrandIDs);
    for (const auto& e : indels)         addEvent("Indels", {e.oldStrandID});

    // Dicentrics are the one exception to counting every strand involved: a dicentric is counted only on the chromosome
    // slot where the karyogram draws it (the one showing two centromeres). The partner chromosome, which is left acentric
    // (e.g. the donor in a deletion-insertion), is not counted. The slot is found with the same rules the drawing uses:
    // the dicentric/acentric pairing override when there is one, otherwise determineHomeStrandID().
    const std::map<int, int> drawnSlots = findDicentricAcentricHomeOverrides(subHeader, numOriginalStrands, masterHeader);

    for (const SDRdicentricEvent& e : dicentrics)
    {
        const auto overrideIt = drawnSlots.find(e.newStrandID);

        if (overrideIt != drawnSlots.end())
        {
            addEvent("Dicentrics", {overrideIt->second});
            continue;
        }

        for (const SDRdataRecord& record : subHeader.dataRecords)
        {
            if (record.newStrandID == e.newStrandID)
            {
                addEvent("Dicentrics", {determineHomeStrandID(record)});
                break;
            }
        }
    }

    // Unexplained changes: a change to a strand that none of the detectors above accounts for, e.g. a terminal inversion, scored whether or not the strand
    // also has detected mutations. Each such strand counts as one Complex Genomic Rearrangement, in the cell total and the karyogram's box as well.
    // Each strand's material in each rearranged record (other than a plain, whole, forward copy) must be accounted for by an event, and so must each of
    // its reversed fragments and each of its gaps. A strand with anything left over gets one Complex Genomic Rearrangements entry.
    {
        const double tolerance = svToleranceMbp;				// 99 bases missing at most, the loosest detector tolerance
        typedef std::pair<int, int> RecordStrand;				// {newStrandID, oldStrandID}
        typedef std::array<double, 3> StrandRange;				// {oldStrandID, low, high}

        std::set<RecordStrand> explainedPieces;				// This strand's material in this record belongs to an event
        std::set<RecordStrand> fullyExplained;				// ...including any reversals in it (multi-strand events)
        std::vector<StrandRange> explainedReversals;			// Inverted segments accounted for
        std::vector<StrandRange> explainedGaps;				// Missing segments accounted for
        std::set<int> wholeStrandEvents;					// Strands whose every change is part of one large event (chromoplexy, chromothripsis)

        auto explainRecords = [&](const std::vector<int>& newStrandIDs, const std::vector<int>& oldStrandIDs, bool full)
        {
            for (int newID : newStrandIDs)
            {
                for (int oldID : oldStrandIDs)
                {
                    explainedPieces.insert({newID, oldID});
                    if (full)
                    {
                        fullyExplained.insert({newID, oldID});
                    }
                }
            }
        };
        auto range = [](int strandID, double a, double b) -> StrandRange
        {
            return {static_cast<double>(strandID), std::min(a, b), std::max(a, b)};
        };

        for (const auto& e : deletions)
        {
            explainRecords({e.remainingStrandID}, {e.oldStrandID}, false);
            explainRecords({e.excisedStrandID}, {e.oldStrandID}, true);	// The excised record is only the deleted piece, so an inverted listing of it is explained too
            explainedGaps.push_back(range(e.oldStrandID, e.deletionStart, e.deletionEnd));
        }
        for (const auto& e : inversions)
        {
            explainRecords({e.newStrandID}, {e.oldStrandID}, false);
            explainedReversals.push_back(range(e.oldStrandID, e.inversionStart, e.inversionEnd));
        }
        for (const auto& e : translocations)
        {
            explainRecords({e.newStrandID1, e.newStrandID2}, {e.oldStrandAid, e.oldStrandBid}, true);
        }
        std::vector<StrandRange> ecDNAranges;
        for (const auto& e : ecDNAs)
        {
            for (std::size_t i = 0; i < e.oldStrandIDs.size(); ++i)
            {
                explainRecords({e.remainingStrandIDs[i]}, {e.oldStrandIDs[i]}, false);
                ecDNAranges.push_back(range(e.oldStrandIDs[i], e.ecDNAsegments[i].first, e.ecDNAsegments[i].second));

                // A deletion-inversion whose deleted piece closed into a circle is scored as ecDNA only (the higher tier absorbs it), so an
                // inverted fragment right next to the excised segment on the remaining strand belongs to this ecDNA as well.
                const double segLow = std::min(e.ecDNAsegments[i].first, e.ecDNAsegments[i].second);
                const double segHigh = std::max(e.ecDNAsegments[i].first, e.ecDNAsegments[i].second);

                for (const SDRdataRecord& record : subHeader.dataRecords)
                {
                    if (record.newStrandID != e.remainingStrandIDs[i])
                    {
                        continue;
                    }

                    for (const SDRfragment& fragment : fragmentsInReadingOrder(record, masterHeader))
                    {
                        const double low = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
                        const double high = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

                        if (fragment.oldStrandID == e.oldStrandIDs[i] && isReversedFragment(fragment)
                            && (approxEqual(low, segHigh, tolerance) || approxEqual(high, segLow, tolerance)))
                        {
                            explainedReversals.push_back(range(fragment.oldStrandID, low, high));
                        }
                    }
                }
            }
        }
        explainedGaps.insert(explainedGaps.end(), ecDNAranges.begin(), ecDNAranges.end());
        for (const auto& e : delInvs)
        {
            explainRecords({e.remainingStrandID}, {e.oldStrandID}, false);
            explainRecords({e.excisedStrandID}, {e.oldStrandID}, true);	// Only the deleted piece, so an inverted listing of it is explained too
            explainedGaps.push_back(range(e.oldStrandID, e.deletionStart, e.deletionEnd));
            explainedReversals.push_back(range(e.oldStrandID, e.inversionStart, e.inversionEnd));
        }
        for (const auto& e : delTras)
        {
            explainRecords({e.newStrandID1, e.newStrandID2}, {e.oldStrandAid, e.oldStrandBid}, true);
            explainRecords({e.excisedStrandID}, {e.deletedOldStrandID}, true);	// Only the deleted piece, so an inverted listing of it is explained too
            explainedGaps.push_back(range(e.deletedOldStrandID, e.deletionStart, e.deletionEnd));
        }
        for (const auto& e : delIns)
        {
            explainRecords({e.donorRemainingNewStrandID}, {e.donorOldStrandID}, false);
            explainRecords({e.recipientNewStrandID}, {e.donorOldStrandID, e.recipientOldStrandID}, true);
            explainedGaps.push_back(range(e.donorOldStrandID, e.segmentStart, e.segmentEnd));
        }
        for (const auto& e : chromoplexy)
        {
            explainRecords(e.contributingNewStrandIDs, e.involvedStrandIDs, true);
            wholeStrandEvents.insert(e.involvedStrandIDs.begin(), e.involvedStrandIDs.end());
        }
        for (const auto& e : chromothripsis)
        {
            explainRecords(e.contributingNewStrandIDs, e.involvedStrandIDs, true);
            wholeStrandEvents.insert(e.involvedStrandIDs.begin(), e.involvedStrandIDs.end());
        }
        for (const auto& e : complexRearrangements)
        {
            explainRecords(e.contributingNewStrandIDs, e.involvedStrandIDs, true);
            wholeStrandEvents.insert(e.involvedStrandIDs.begin(), e.involvedStrandIDs.end());
        }

        // Rings and dicentrics account for their own record, and for the acentric leftovers of the same strands (the ends cut off when a ring
        // closes, the acentric partner of a dicentric), so those leftovers are not scored again. Their material also fills the strands' gaps.
        auto explainWithLeftovers = [&](int newStrandID, const std::vector<int>& strandIDs)
        {
            explainRecords({newStrandID}, strandIDs, true);
            const std::set<int> strands(strandIDs.begin(), strandIDs.end());

            for (const SDRdataRecord& record : subHeader.dataRecords)
            {
                if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))
                {
                    continue;
                }

                bool onlyTheseStrands = !record.fragments.empty();
                bool acentric = true;

                for (const SDRfragment& fragment : record.fragments)
                {
                    onlyTheseStrands = onlyTheseStrands && strands.count(fragment.oldStrandID);
                    acentric = acentric && !fragment.hasCentromere;
                }

                if (onlyTheseStrands && (acentric || record.newStrandID == newStrandID))
                {
                    explainRecords({record.newStrandID}, strandIDs, record.newStrandID == newStrandID);
                }

                if (record.newStrandID == newStrandID)
                {
                    for (const SDRfragment& fragment : record.fragments)
                    {
                        explainedGaps.push_back(range(fragment.oldStrandID, fragment.oldStartPosition, fragment.oldEndPosition));
                    }
                }
            }
        };
        for (const auto& e : rings)
        {
            explainWithLeftovers(e.newStrandID, e.involvedOldStrandIDs);
        }
        for (const auto& e : dicentrics)
        {
            explainWithLeftovers(e.newStrandID, e.involvedOldStrandIDs);
        }

        // A reversed fragment or gap counts as accounted for when an event's segment overlaps it (detectors trim slightly differently).
        auto overlapsAny = [&](const std::vector<StrandRange>& ranges, int strandID, double low, double high)
        {
            for (const StrandRange& r : ranges)
            {
                if (static_cast<int>(r[0]) == strandID && r[1] < high - tolerance && low < r[2] - tolerance)
                {
                    return true;
                }
            }
            return false;
        };

        // The circles themselves: every circular fragment lying in an ecDNA segment (one formation can have several circles).
        for (const SDRdataRecord& record : subHeader.dataRecords)
        {
            if (record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands))
            {
                continue;
            }

            for (const SDRfragment& fragment : record.fragments)
            {
                const double low = std::min(fragment.oldStartPosition, fragment.oldEndPosition);
                const double high = std::max(fragment.oldStartPosition, fragment.oldEndPosition);

                for (const StrandRange& r : ecDNAranges)
                {
                    if (static_cast<int>(r[0]) == fragment.oldStrandID && low >= r[1] - tolerance && high <= r[2] + tolerance)
                    {
                        explainRecords({record.newStrandID}, {fragment.oldStrandID}, false);

                        // An ecDNA can be closed from forward and reversed pieces (e.g. 17/111.9/112.8 + 17/118.9/112.8), and a circle has no
                        // orientation, so a reversed piece inside the ecDNA is part of it and not an inversion left unexplained.
                        if (isReversedFragment(fragment))
                        {
                            explainedReversals.push_back(range(fragment.oldStrandID, low, high));
                        }
                    }
                }
            }
        }

        std::set<int> unexplainedStrands;
        std::map<int, std::vector<SDRfragment>> pooledByStrand;		// Every rearranged fragment per strand, for the gap check below
        std::set<int> wholeCopyStrands;						// Strands with an unchanged whole copy among the rearranged records

        for (const SDRdataRecord& record : subHeader.dataRecords)
        {
            if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))
            {
                continue;
            }

            if (record.linear && record.fragments.size() == 1)
            {
                const SDRfragment& fragment = record.fragments[0];
                const double strandSize = fullSizeForStrand(masterHeader, fragment.oldStrandID);

                if (approxEqual(fragment.oldStartPosition, 0.0, tolerance) && strandSize > 0.0 && approxEqual(fragment.oldEndPosition, strandSize, tolerance))
                {
                    wholeCopyStrands.insert(fragment.oldStrandID);
                    continue;							// An unchanged copy of the whole chromosome, nothing altered.
                }
            }

            for (const SDRfragment& fragment : fragmentsInReadingOrder(record, masterHeader))	// Upside-down listings are judged as read the normal way
            {
                const int strandID = fragment.oldStrandID;
                pooledByStrand[strandID].push_back(fragment);

                if (wholeStrandEvents.count(strandID))
                {
                    continue;
                }

                const RecordStrand key{record.newStrandID, strandID};

                if (!explainedPieces.count(key))
                {
                    unexplainedStrands.insert(strandID);
                }
                else if (isReversedFragment(fragment) && !fullyExplained.count(key)
                         && !overlapsAny(explainedReversals, strandID, fragment.oldEndPosition, fragment.oldStartPosition))
                {
                    unexplainedStrands.insert(strandID);			// An inversion no event accounts for, e.g. a terminal one next to a balanced inversion.
                }
            }
        }

        // Missing material: every gap left in a strand's pooled rearranged fragments must belong to a deletion-type event (or be filled by a ring).
        for (const auto& [strandID, fragments] : pooledByStrand)
        {
            if (wholeStrandEvents.count(strandID) || wholeCopyStrands.count(strandID))
            {
                continue;
            }

            std::vector<SDRnormalizedFragment> sorted;
            std::vector<SDRflankGap> gaps;

            if (!analyzeFlankingFragments(fragments, fullSizeForStrand(masterHeader, strandID), tolerance, sorted, gaps))
            {
                continue;								// Overlapping material (e.g. duplicated); not judged here.
            }

            for (const SDRflankGap& gap : gaps)
            {
                if (!overlapsAny(explainedGaps, strandID, gap.start, gap.end))
                {
                    unexplainedStrands.insert(strandID);
                }
            }
        }

        for (int strandID : unexplainedStrands)
        {
            addEvent("Complex Genomic Rearrangements", {strandID});	// Safeguard: never doubles a CGR cluster's entry, its strands are skipped above (wholeStrandEvents)
        }

        if (unexplainedStrandsOut)
        {
            *unexplainedStrandsOut = unexplainedStrands;
        }
    }

    // Same order as the "Mutation Summary" in the SDR summary file.
    const std::vector<std::string> typeOrder = {
        "Long Deletions", "Balanced Inversions", "Balanced Translocations", "Extrachromosomal DNA (ecDNA)", "Deletion-Inversions",
        "Deletion-Translocations", "Deletion-Insertions", "Chromoplexy", "Chromothripsis", "Complex Genomic Rearrangements", "Ring Chromosomes", "Dicentrics",
        "Indels"};

    std::map<int, std::vector<std::pair<std::string, int>>> result;

    for (const auto& [strandID, typeCounts] : counts)
    {
        for (const std::string& type : typeOrder)
        {
            const auto it = typeCounts.find(type);

            if (it != typeCounts.end())
            {
                result[strandID].emplace_back(type, it->second);
            }
        }
    }

    return result;
}


// Number of complex genomic rearrangements in one cell, for the Mutation Summary and the karyogram's box: the 1-2 strand clusters with 4+ pieces out
// of order or shattered exchanged material (detectComplexGenomicRearrangements()), plus one for every chromosome with a change no other mutation type accounts for.
// Indels are taken out first (see removeIndels()), so they never count here.
inline std::size_t countComplexGenomicRearrangements(const SDRsubHeader& rawSubHeader, const SDRmasterHeader& masterHeader)
{
    const int numOriginalStrands = masterHeader.intactChromosomeSizes.empty() ? 0 : static_cast<int>(masterHeader.intactChromosomeSizes[0]);
    const SDRsubHeader subHeader = removeIndels(rawSubHeader, numOriginalStrands);

    std::set<int> unexplainedStrands;
    countMutationsPerOriginalStrand(subHeader, masterHeader, &unexplainedStrands);

    return detectComplexGenomicRearrangements(subHeader, numOriginalStrands).size() + unexplainedStrands.size();
}

} // namespace karyoke








#endif
