#include <cmath>
#include <algorithm>
#include <map>
#include <string>
#include <set>
#include <iostream>

#include "Karyogram.h"
#include "SDRutilities.h"



namespace karyoke
{

// ------------------------------------------------------- //
// Helper functions to determine shape of mutations and
// how they should be drawn
// ------------------------------------------------------- //




// Recognizes a simple chromosome ring formation: exactly 3 records in one slot: one CIRCULAR record carrying a centromere (the ring
// itself), and two LINEAR, single-fragment, acentric records (the two excised chromosome ends) that, together with the ring, exactly tile
// the home strand's full declared length, with the ring genuinely in the middle. Scoped to this single-ring, two-end case only.
bool Karyogram::isChromosomeRingShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, const SDRmasterHeader& masterHeader)
{
    if (records.size() != 3)
    {
        return false;
    }

    const SDRdataRecord* ringRecord = nullptr;
    std::vector<const SDRdataRecord*> excisedRecords;

    for (const SDRdataRecord* record : records)
    {
        if (!record->linear)
        {
            if (ringRecord != nullptr || record->fragments.size() != 1)
            {
                return false; // Only ONE circular record expected, and scoped to a single-fragment ring for now.
            }

            const SDRfragment& fragment = record->fragments[0];

            if (fragment.oldStrandID != homeOldStrandID || !fragment.hasCentromere)
            {
                return false; // The ring itself must carry the centromere and originate from the home strand.
            }

            ringRecord = record;
        }
        else
        {
            if (record->fragments.size() != 1)
            {
                return false; // Scoped to simple, single-fragment excised ends.
            }

            const SDRfragment& fragment = record->fragments[0];

            if (fragment.oldStrandID != homeOldStrandID || fragment.hasCentromere || isReversedFragment(fragment))
            {
                return false; // Excised ends must be plain, acentric, non-reversed home-strand material.
            }

            excisedRecords.push_back(record);
        }
    }

    if (ringRecord == nullptr || excisedRecords.size() != 2)
    {
        return false;
    }

    const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);

    if (sizeIndex >= masterHeader.intactChromosomeSizes.size())
    {
        return false;
    }

    const double fullSize = masterHeader.intactChromosomeSizes[sizeIndex];

    if (fullSize <= 0.0)
    {
        return false;
    }

    // Tag each piece so we can confirm the ring lands genuinely in the MIDDLE once sorted, both ends must flank it, not sit on one side.
    struct TaggedFragment
    {
        SDRfragment fragment;
        bool isRing;
    };

    std::vector<TaggedFragment> allPieces = {
        {ringRecord->fragments[0], true},
        {excisedRecords[0]->fragments[0], false},
        {excisedRecords[1]->fragments[0], false}
    };

    std::sort(allPieces.begin(), allPieces.end(), [](const TaggedFragment& a, const TaggedFragment& b)
    {
        return std::min(a.fragment.oldStartPosition, a.fragment.oldEndPosition) < std::min(b.fragment.oldStartPosition, b.fragment.oldEndPosition);
    });

    if (!allPieces[1].isRing)
    {
        return false;
    }

    const double delTolerance = 0.000050;					// 50 bp = 0.000050 Mbp. Minimum tolerance for a long deletion SV to occur.

    const double firstLo = std::min(allPieces.front().fragment.oldStartPosition, allPieces.front().fragment.oldEndPosition);
    const double lastHi = std::max(allPieces.back().fragment.oldStartPosition, allPieces.back().fragment.oldEndPosition);

    if (!approxEqual(firstLo, 0.0, delTolerance) || !approxEqual(lastHi, fullSize, delTolerance))
    {
        return false; // Must span the whole chromosome, start to end.
    }

    for (std::size_t i = 0; i + 1 < allPieces.size(); ++i)
    {
        const double hi = std::max(allPieces[i].fragment.oldStartPosition, allPieces[i].fragment.oldEndPosition);
        const double lo = std::min(allPieces[i + 1].fragment.oldStartPosition, allPieces[i + 1].fragment.oldEndPosition);

        if (!approxEqual(hi, lo, delTolerance))
        {
            return false; // Must tile with no gaps or overlaps between the three pieces.
        }
    }

    return true;
}













// Recognizes the DONOR side of a deletion-translocation event: two
// records in one slot, a "recombined" record (fragments from the
// home strand PLUS exactly one foreign fragment from elsewhere) and
// a "excised" record (a single home-strand-only fragment). Unlike
// isDeletionShape(), the recombined record is allowed foreign
// material, since the swapped-in piece isn't a gap to fill in, it's
// already correctly part of the recombined molecule as drawn.
bool Karyogram::isDeletionTranslocationDonorShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID)
{
    if (records.size() != 2)
    {
        return false;
    }

    const SDRdataRecord* recombinedRecord = nullptr;
    const SDRdataRecord* excisedRecord = nullptr;

    for (const SDRdataRecord* record : records)
    {
        if (!record->linear)
        {
            return false;
        }

        int homeFragmentCount = 0;
        int foreignFragmentCount = 0;

        for (const SDRfragment& fragment : record->fragments)
        {
            if (fragment.oldStrandID == homeOldStrandID)
            {
                if (isReversedFragment(fragment))
                {
                    return false; // The home fragment itself should not be reversed, a real deletion/translocation donor's own material stays in orientation.
                }

                ++homeFragmentCount;
            }
            else
            {
                // Foreign fragment MAY be reversed, represents an
                // inverted exchange (e.g. forming a dicentric/acentric
                // pair), not a reason to reject this shape.
                ++foreignFragmentCount;
            }
        }

        if (homeFragmentCount >= 1 && foreignFragmentCount >= 1 && recombinedRecord == nullptr)
        {
            recombinedRecord = record;
        }
        else if (homeFragmentCount == 1 && foreignFragmentCount == 0 && excisedRecord == nullptr)
        {
            excisedRecord = record;
        }
        else
        {
            return false;
        }
    }

    if (recombinedRecord == nullptr || excisedRecord == nullptr)
    {
        return false;
    }

    // A recombined record that ALSO holds a real gap between two adjacent home-strand fragments (an independent deletion, e.g. the deletion beside a deletion-insertion's
    // inserted piece) is not this shape: its gap has to be drawn as white space in place, which the generic mixed-record drawing does.
    for (std::size_t i = 1; i < recombinedRecord->fragments.size(); ++i)
    {
        const SDRfragment& previous = recombinedRecord->fragments[i - 1];
        const SDRfragment& current = recombinedRecord->fragments[i];

        if (previous.oldStrandID != current.oldStrandID)
        {
            continue;
        }

        if (current.oldStrandID != homeOldStrandID && (isReversedFragment(previous) || isReversedFragment(current) || previous.oldEndPosition > current.oldStartPosition))	// Foreign material is only trusted in ascending order.
        {
            continue;
        }

        const double tolerance = 0.000050;
        const double prevLo = std::min(previous.oldStartPosition, previous.oldEndPosition);
        const double prevHi = std::max(previous.oldStartPosition, previous.oldEndPosition);
        const double curLo = std::min(current.oldStartPosition, current.oldEndPosition);
        const double curHi = std::max(current.oldStartPosition, current.oldEndPosition);

        if ((prevHi <= curLo && !approxEqual(prevHi, curLo, tolerance)) || (curHi <= prevLo && !approxEqual(curHi, prevLo, tolerance)))
        {
            return false;
        }
    }

    return true;
}













/*

// Recognizes the DONOR side of a deletion-insertion event when it's
// drawn alone in its slot: a single home-strand-only record with 2+
// fragments and at least one real gap between them. The material
// that left this strand shows up as a foreign fragment inside a
// completely different slot's record (the recipient), so there's no
// matching excised record here to pair with, this shape has to be
// recognized from a single record alone.
bool Karyogram::isLoneGapShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID)
{
    if (records.size() != 1)
    {
        return false;
    }

    const SDRdataRecord* record = records[0];

    if (!record->linear)
    {
        return false;
    }

    if (record->fragments.size() < 2)
    {
        return false;
    }

    for (const SDRfragment& fragment : record->fragments)
    {
        if (fragment.oldStrandID != homeOldStrandID)
        {
            return false;
        }

        if (isReversedFragment(fragment))
        {
            return false; // Not this shape - an inversion, not a deletion-insertion donor.
        }
    }

    std::vector<SDRfragment> sortedFragments = record->fragments;

    std::sort(sortedFragments.begin(), sortedFragments.end(), [](const SDRfragment& a, const SDRfragment& b)
    {
        return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
    });

    const double tolerance = 0.000050;
    bool hasGap = false;

    for (std::size_t i = 0; i + 1 < sortedFragments.size(); ++i)
    {
        const double higherPosStrandA = std::max(sortedFragments[i].oldStartPosition, sortedFragments[i].oldEndPosition);
        const double lowerPosStrandB = std::min(sortedFragments[i + 1].oldStartPosition, sortedFragments[i + 1].oldEndPosition);

        if (lowerPosStrandB < higherPosStrandA && !approxEqual(higherPosStrandA, lowerPosStrandB, tolerance))
        {
            return false; // Genuine overlap - not a valid shape.
        }

        if (lowerPosStrandB > higherPosStrandA && !approxEqual(higherPosStrandA, lowerPosStrandB, tolerance))
        {
            hasGap = true;
        }
    }

    return hasGap;
}
*/











/*
// Same pairing pattern as isDeletionShape()/isECDNAshape(), but for a
// deletion-inversion: one record has THREE fragments (2 normal + 1
// reversed) from the home strand, with a real gap between two of
// them; a second record supplies the single fragment filling that
// gap. All fragments (both records) must be linear.
bool Karyogram::isDeletionInversionShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID)
{
    if (records.size() != 2)
    {
        return false;
    }

    const SDRdataRecord* threeFragmentRecord = nullptr;
    const SDRdataRecord* oneFragmentRecord = nullptr;

    for (const SDRdataRecord* record : records)
    {
        if (!record->linear)
        {
            return false;
        }

        for (const SDRfragment& fragment : record->fragments)
        {
            if (fragment.oldStrandID != homeOldStrandID)
            {
                return false;
            }
        }

        if (record->fragments.size() == 3 && threeFragmentRecord == nullptr)
        {
            threeFragmentRecord = record;
        }
        else if (record->fragments.size() == 1 && oneFragmentRecord == nullptr)
        {
            oneFragmentRecord = record;
        }
        else
        {
            return false;
        }
    }

    if (threeFragmentRecord == nullptr || oneFragmentRecord == nullptr)
    {
        return false;
    }



    int reversedCount = 0;

    for (const SDRfragment& fragment : threeFragmentRecord->fragments)
    {
        if (fragment.oldStartPosition > fragment.oldEndPosition)
        {
            ++reversedCount;
        }
    }

    if (reversedCount != 1)
    {
        return false;
    }


    std::vector<SDRfragment> sortedThreeFragmentRecord = threeFragmentRecord->fragments;

    std::sort(sortedThreeFragmentRecord.begin(), sortedThreeFragmentRecord.end(), [](const SDRfragment& a, const SDRfragment& b)
    {
        return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
    });

    const SDRfragment& excisedFragment = oneFragmentRecord->fragments[0];

    if (excisedFragment.oldStartPosition > excisedFragment.oldEndPosition)
    {
        return false; // The excised piece itself should not be reversed.
    }

    const double delTolerance = 0.000050;
    int gapCount = 0;
    bool gapMatchesExcised = false;

    for (std::size_t i = 0; i + 1 < sortedThreeFragmentRecord.size(); ++i)
    {
        const double higherPosStrandA = std::max(sortedThreeFragmentRecord[i].oldStartPosition, sortedThreeFragmentRecord[i].oldEndPosition);
        const double lowerPosStrandB = std::min(sortedThreeFragmentRecord[i + 1].oldStartPosition, sortedThreeFragmentRecord[i + 1].oldEndPosition);

        if (approxEqual(higherPosStrandA, lowerPosStrandB, delTolerance))
        {
            continue;
        }

        if (lowerPosStrandB <= higherPosStrandA)
        {
            return false; // Overlapping, not a gap - invalid shape.
        }

        ++gapCount;

        if (approxEqual(excisedFragment.oldStartPosition, higherPosStrandA, delTolerance) &&
            approxEqual(excisedFragment.oldEndPosition, lowerPosStrandB, delTolerance))
        {
            gapMatchesExcised = true;
        }
    }

    return gapCount == 1 && gapMatchesExcised;
}
*/











/*


// Identical shape to isDeletionShape(), except the excised (single fragment)
// records must be CIRCULAR rather than linear, and must be acentric -
// an ecDNA fragment carrying a centromere isn't ecDNA. Generalized to N+1
// ecDNA fragments on the same original strand. So instead of requiring exactly
// one fragment per excised record, this checks that the TOTAL
// fragment count across all excised records matches the total number
// of gaps in the remaining record.

bool Karyogram::isECDNAshape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, const SDRmasterHeader& masterHeader)
{
    if (records.size() < 2)
    {
        return false;
    }

    const SDRdataRecord* multipleFragmentRecord = nullptr;
    std::vector<const SDRdataRecord*> excisedRecords;

    for (const SDRdataRecord* record : records)
    {
        for (const SDRfragment& fragment : record->fragments)
        {
            if (fragment.oldStrandID != homeOldStrandID)
            {
                return false;
            }
        }

        if (record->fragments.size() >= 2 && record->linear && multipleFragmentRecord == nullptr)
        {
            multipleFragmentRecord = record;
        }
        else if (!record->linear && !record->fragments.empty())
        {
            bool hasCentromereFlag = false;

            for (const SDRfragment& fragment : record->fragments)
            {
                if (fragment.hasCentromere)
                {
                    hasCentromereFlag = true;
                    break;
                }
            }

            if (hasCentromereFlag)
            {
                return false;
            }

            excisedRecords.push_back(record);
        }
        else
        {
            return false;
        }
    }

    if (multipleFragmentRecord == nullptr || excisedRecords.empty())
    {
        return false;
    }


    std::size_t totalExcisedFragments = 0;

    for (const SDRdataRecord* record : excisedRecords)
    {
        totalExcisedFragments += record->fragments.size();
    }



    std::vector<SDRfragment> sortedFragments = multipleFragmentRecord->fragments;

    std::sort(sortedFragments.begin(), sortedFragments.end(), [](const SDRfragment& a, const SDRfragment& b)
    {
        return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
    });


    const double delTolerance = 0.000050;
    std::vector<std::pair<double, double>> gaps;

    const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);
    const double chromosomeFullSizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size()) ? masterHeader.intactChromosomeSizes[sizeIndex] : 0.0;

    // ------------------------------------------------------------- //
    // Determine where the deleted segment came from on the chromosome //
    // ------------------------------------------------------------- //

    // Check if fragment came from the start of the chromosome
    if (chromosomeFullSizeMbp > 0.0)
    {
        const double firstFragmentStart = std::min(sortedFragments.front().oldStartPosition, sortedFragments.front().oldEndPosition);

        if (!approxEqual(firstFragmentStart, 0.0, delTolerance) && firstFragmentStart > 0.0)
        {
            gaps.push_back({0.0, firstFragmentStart});
        }
    }

    // Check if fragment came from the middle of the chromosome
    for (std::size_t i = 0; i + 1 < sortedFragments.size(); i++)
    {
        const double higherPosStrandA = std::max(sortedFragments[i].oldStartPosition, sortedFragments[i].oldEndPosition);
        const double lowerPosStrandB = std::min(sortedFragments[i + 1].oldStartPosition, sortedFragments[i + 1].oldEndPosition);

        if (approxEqual(higherPosStrandA, lowerPosStrandB, delTolerance))
        {
            continue;
        }

        if (lowerPosStrandB <= higherPosStrandA)
        {
            return false; // Overlapping, not a gap - invalid shape.
        }

        gaps.push_back({higherPosStrandA, lowerPosStrandB});
    }

    // Check if fragment came from the end of the chromosome

    if (chromosomeFullSizeMbp > 0.0)
    {
        const double lastFragmentEnd = std::max(sortedFragments.back().oldStartPosition, sortedFragments.back().oldEndPosition);

        if (chromosomeFullSizeMbp > lastFragmentEnd && !approxEqual(lastFragmentEnd, chromosomeFullSizeMbp, delTolerance))
        {
            gaps.push_back({lastFragmentEnd, chromosomeFullSizeMbp});
        }
    }



    if (gaps.size() != totalExcisedFragments)
    {
        return false;
    }


    std::vector<bool> gapClaimed(gaps.size(), false);

    for (const SDRdataRecord* excisedRecord : excisedRecords)
    {
        for (const SDRfragment& excisedFragment : excisedRecord->fragments)
        {
            if (excisedFragment.oldStartPosition > excisedFragment.oldEndPosition)
            {
                return false;
            }


            bool matched = false;

            for (std::size_t i = 0; i < gaps.size(); ++i)
            {
                if (gapClaimed[i])
                {
                    continue;
                }

                if (approxEqual(excisedFragment.oldStartPosition, gaps[i].first, delTolerance) &&
                    approxEqual(excisedFragment.oldEndPosition, gaps[i].second, delTolerance))
                {
                    gapClaimed[i] = true;
                    matched = true;
                    break;
                }
            }

            if (!matched)
            {
                return false;
            }
        }
    }

    return true;

}


*/









/*

// Generalized deletion shape to account for multiple deletions on the same strand.
// Shape looks like this, one record contains N + 1 fragments, where N is the number of
// deletions on that original strand. All N + 1 fragments have the same old strand ID.
// Then this record is followed by N records representing the N deletions, each with one
// fragment in the record. All deletion records have one fragment with the same old strand ID.
bool Karyogram::isDeletionShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, const SDRmasterHeader& masterHeader)
{
    if (records.size() < 2)
    {
        return false;
    }

    const SDRdataRecord* multipleFragmentRecord = nullptr;
    std::vector<const SDRdataRecord*> singleFragmentRecords;

    for (const SDRdataRecord* record : records)
    {

        if(!record->linear)
        {
            return false;
        }

        for (const SDRfragment& fragment : record->fragments)
        {
            if (fragment.oldStrandID != homeOldStrandID)
            {
                return false;
            }
        }

        if (record->fragments.size() >= 2 && multipleFragmentRecord == nullptr)
        {
            multipleFragmentRecord = record;
        }
        else if (record->fragments.size() == 1)
        {
            singleFragmentRecords.push_back(record);
        }
        else
        {
            return false;
        }
    }


    // ------------------------------------------------------------ //
    // Determine if deletion is in the middle or at the ends of a chromosome
    // ------------------------------------------------------------ //
    std::vector<SDRfragment> sortedFragments = multipleFragmentRecord->fragments;

    std::sort(sortedFragments.begin(), sortedFragments.end(), [](const SDRfragment& a, const SDRfragment& b)
    {
        return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
    });

    const double delTolerance = 0.000050;
    std::vector<std::pair<double, double>> gaps;


    // A TERMINAL deletion (touching position 0 or the chromosome's true
    // end) leaves no extra "remaining" fragment on that side, so the
    // old rigid N+1 formula only held when every deletion was strictly
    // internal. Checking against the chromosome's true declared length
    // (rather than just gaps between listed fragments) generalizes
    // this to terminal deletions too.
    const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);
    const double chromosomeFullSizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size())
        ? masterHeader.intactChromosomeSizes[sizeIndex]
        : 0.0;

    // Check if deletion is at the start of the chromosome
    if (chromosomeFullSizeMbp > 0.0)
    {
        const double firstFragmentStart = std::min(sortedFragments.front().oldStartPosition, sortedFragments.front().oldEndPosition);

        if (!approxEqual(firstFragmentStart, 0.0, delTolerance) && firstFragmentStart > 0.0)
        {
            gaps.push_back({0.0, firstFragmentStart});
        }
    }



    // Check for deletions in central segments
    for (std::size_t i = 0; i + 1 < sortedFragments.size(); ++i)
    {
        const double higherPosStrandA = std::max(sortedFragments[i].oldStartPosition, sortedFragments[i].oldEndPosition);
        const double lowerPosStrandB = std::min(sortedFragments[i + 1].oldStartPosition, sortedFragments[i + 1].oldEndPosition);

        if (approxEqual(higherPosStrandA, lowerPosStrandB, delTolerance))
        {
            continue;
        }

        if (lowerPosStrandB <= higherPosStrandA)
        {
            return false; // Overlapping, not a gap - invalid shape.
        }

        gaps.push_back({higherPosStrandA, lowerPosStrandB});
    }


    // Check if deletion at the end of the chromosome
    if (chromosomeFullSizeMbp > 0.0)
    {
        const double lastFragmentEnd = std::max(sortedFragments.back().oldStartPosition, sortedFragments.back().oldEndPosition);

        if (chromosomeFullSizeMbp > lastFragmentEnd && !approxEqual(lastFragmentEnd, chromosomeFullSizeMbp, delTolerance))
        {
            gaps.push_back({lastFragmentEnd, chromosomeFullSizeMbp});
        }
    }




    if (gaps.size() != singleFragmentRecords.size())
    {
        return false;
    }

    std::vector<bool> gapClaimed(gaps.size(), false);

    for (const SDRdataRecord* excisedRecord : singleFragmentRecords)
    {
        const SDRfragment& excisedFragment = excisedRecord->fragments[0];

        if (excisedFragment.oldStartPosition > excisedFragment.oldEndPosition)
        {
            return false; // Excised pieces should not be reversed.
        }

        bool matched = false;

        for (std::size_t i = 0; i < gaps.size(); ++i)
        {
            if (gapClaimed[i])
            {
                continue;
            }

            if (approxEqual(excisedFragment.oldStartPosition, gaps[i].first, delTolerance) &&
                approxEqual(excisedFragment.oldEndPosition, gaps[i].second, delTolerance))
            {
                gapClaimed[i] = true;
                matched = true;
                break;
            }
        }

        if (!matched)
        {
            return false;
        }
    }

    return true;
}

*/















// ------------------------------------------------------------------------------------------------------------------------------ //
// Unified "home-strand gap" shape. Replaces the separate deletion / ecDNA / deletion-inversion / lone-gap shape checks for drawing.
//
// A slot qualifies when EVERY fragment in it comes from the home strand and it holds:
//   - exactly ONE "remaining" record: linear, 2+ fragments (below the chromothripsis fragment threshold), with at least one gap
//     (internal or terminal) between its fragments,
//   - zero or more "pieces": linear single-fragment excised records (not reversed), or acentric circular (ecDNA) records.
// Pieces are OPTIONAL and need not match any gap: the SDR file draws what it implies. A gap with no piece is just white space.
// Uses the same gap logic (analyzeFlankingFragments) and 50 bp tolerance as the detectors in SDRutilities.h.
// ------------------------------------------------------------------------------------------------------------------------------ //
bool Karyogram::analyzeHomeStrandGapShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, const SDRmasterHeader& masterHeader,
                                          const SDRdataRecord*& remainingRecord, std::vector<const SDRdataRecord*>& pieceRecords)
{
    const double delTolerance = 0.000050;			// Same as the detectors
    const std::size_t chromothripsisFragmentThreshold = 10;	// Same default as detectChromothripsis(); such strands are drawn as concatenated shatter, not as gaps.

    remainingRecord = nullptr;
    pieceRecords.clear();

    auto reject = [&]() -> bool
    {
        remainingRecord = nullptr;
        pieceRecords.clear();
        return false;
    };

    if (records.empty())
    {
        return false;
    }

    int multiFragmentLinearCount = 0;
    std::vector<const SDRdataRecord*> singleCentricLinear;		// Single-fragment, centromere-bearing, non-reversed linear records: terminal-deletion remaining candidates

    for (const SDRdataRecord* record : records)
    {
        for (const SDRfragment& fragment : record->fragments)
        {
            if (fragment.oldStrandID != homeOldStrandID)
            {
                return reject();
            }
        }

        if (record->fragments.empty())
        {
            return reject();
        }

        if (record->linear && record->fragments.size() >= 2)
        {
            ++multiFragmentLinearCount;
            remainingRecord = record;
        }
        else if (record->linear && record->fragments[0].hasCentromere && !isReversedFragment(record->fragments[0]))
        {
            singleCentricLinear.push_back(record);
        }
    }

    if (multiFragmentLinearCount > 1)
    {
        return reject();
    }

    // Terminal deletions leave a single-fragment remaining record, only accepted when the strand has no multi-fragment remaining record.
    if (multiFragmentLinearCount == 0)
    {
        if (singleCentricLinear.size() != 1)
        {
            return reject();
        }

        remainingRecord = singleCentricLinear[0];
    }

    // Everything else is a PIECE: a linear, single-fragment, non-reversed excised record, or a circular record cut from a gap. A circular piece
    // may be acentric (ecDNA) or keep the centromere (a ring next to an acentric remaining record), drawSingleOriginCircularPiece() draws either.
    for (const SDRdataRecord* record : records)
    {
        if (record == remainingRecord)
        {
            continue;
        }

        if (record->linear)
        {
            if (record->fragments.size() != 1 || isReversedFragment(record->fragments[0]))
            {
                return reject();
            }
        }

        pieceRecords.push_back(record);
    }

    if (remainingRecord->fragments.size() >= chromothripsisFragmentThreshold)
    {
        return reject();
    }

    std::vector<SDRnormalizedFragment> sorted;
    std::vector<SDRflankGap> gaps;

    if (!analyzeFlankingFragments(remainingRecord->fragments, fullSizeForStrand(masterHeader, homeOldStrandID), delTolerance, sorted, gaps) || gaps.empty())
    {
        return reject();
    }

    // Order the pieces top-to-bottom by where their material came from, so the side column reads like the gaps on the bar.
    std::stable_sort(pieceRecords.begin(), pieceRecords.end(), [](const SDRdataRecord* a, const SDRdataRecord* b)
    {
        auto lowest = [](const SDRdataRecord* r)
        {
            double low = r->fragments.front().oldStartPosition;
            for (const SDRfragment& f : r->fragments)
            {
                low = std::min(low, std::min(f.oldStartPosition, f.oldEndPosition));
            }
            return low;
        };
        return lowest(a) < lowest(b);
    });

    return true;
}







} // namespace karyoke
