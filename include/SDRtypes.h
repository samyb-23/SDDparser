#ifndef SDR_TYPES_H
#define SDR_TYPES_H

#include <string>
#include <vector>
#include <map>

namespace karyoke
{

// Represents one fragment used to construct a new SDR strand. - Field 3
struct SDRfragment
{
    // Field 3: oldStrandID/oldStrandStart/oldStrandEnd/hasCentromere
    int oldStrandID;
    double oldStartPosition;
    double oldEndPosition;
    bool hasCentromere;
};


// Represents one strand in the SDR data section.
struct SDRdataRecord
{
    int cellID;						// SDR data field 1
    int newStrandID;					// SDR data field 2
    std::vector<SDRfragment> fragments;			// SDR data field 3, can contain multiple subfields separated by ',' and each subfield contains 4 entries separated by '/'.
    bool linear;
};


// Stores the SDR master header entries
struct SDRmasterHeader
{
    std::string sdrVersion;
    std::string author;
    std::string associatedSDDFile;
    std::vector<double> cellCyclePhase;			// Optional (SDRv2.1+): same coding as SDD 'Cell cycle phase'. Summary only, never used for the karyogram drawing.
    std::vector<double> intactChromosomeSizes;		// First value in the SDR file is the number of chromosomes, followed by the chromosome sizes.
};


// Stores the information contained in one SDR cell subheader.
struct SDRsubHeader
{
    int cellID;

    std::map<int, double> mutatedChromosomeSizes;	// Key = mutated/new strand ID, Value = mutated chromosome/strand size in Mbp

    std::vector<int> intactStrandIDs;			// IDs of strands that remain intact.

    int totalDSBcount;					// A blank value in SDR subheader means not measured, not an error
    int totalMisrepairCount;

    std::string medrasMClog;

    std::vector<SDRdataRecord> dataRecords;		// All SDR data records belonging to this cell.

};



// Use for checking if total DSB count, total Misrepair count, or MEDRAS-MC log are blank, if so pass warning, not error.
constexpr int SDR_FIELD_NOT_MEASURED = -1;




// ------------------------------------------------ //
// Mutation types for summarizing in text file
// ------------------------------------------------ //

// Represents a detected deletion event within a single cell
struct SDRdeletionEvent
{
    int oldStrandID;		// Which chromosome/strand ID the fragments originate from
    double deletionStart;	// Location (in Mbp from start of p arm) for the beginning of the deletion
    double deletionEnd;		// Location (in Mbp from start of p arm) for the beginning of the deletion
    int remainingStrandID;   	// new strand retaining the two flanking fragments
    int excisedStrandID;     	// new strand consisting of just the deleted segment
};


// Represents a detected balanced inversion: a segment of
// the original strand was flipped in orientation in place, with no
// net gain or loss of material.
struct SDRinversionEvent
{
    int oldStrandID;		// Which chromosome/strand ID the fragments originate from
    double inversionStart;
    double inversionEnd;
    int newStrandID;
};


// Represents a detected balanced translocation: two
// original strands were each split once, and the resulting pieces
// were swapped between two new strands, with no material gained or
// lost. --> Generalized to handle swapping of central fragments.
struct SDRtranslocationEvent
{
    int oldStrandAid;				// Which chromosome/strand ID the fragments originate from for strand A
    int oldStrandBid;				// Which chromosome/strand ID the fragments originate from for strand B
    std::vector<double> breakpointsA;   	// position on oldStrandAid where the breaks occurred
    std::vector<double> breakpointsB;   	// position on oldStrandBid where the breaks occurred
    int newStrandID1;				// ID of the new strand containing fragment rearrangements for strand 1
    int newStrandID2;				// ID of the new strand containing fragment rearrangenemts for strand 2
};


// Represents a detected extrachromosomal DNA fragment.
// It has the same fragment signature as long deletions, but
// the isLinear data field 4 is 0 (for circular fragments) rather than
// 1 for linear fragments. Can be made up of multiple fragments from the
// same original strand and now supports ecDNA made up of multiple foreign
// fragments
struct SDRecDNAevent
{
    std::vector<int> oldStrandIDs;				// Which chromosome/strand ID the fragments originate from
    std::vector<std::pair<double, double>> ecDNAsegments; 	// One or more (start, end) segments making up this ecDNA formation (one molecule, or several when one deleted segment shattered into several circles)
    std::vector<int> remainingStrandIDs;			// New strand ID of the strand that lost the excised fragment
    int excisedStrandID;					// New strand ID of the excised strand (the first circle when there are several)
};




// Represents a combined detected deletion-inversion event.
// One new strand contains three fragments with an inversion and a
// missing section. Another new strand contains one fragment that is
// the missing section from the previous new strand.
struct SDRdeletionInversionEvent
{
    int oldStrandID;		// Which chromosome/strand ID the fragments originate from
    double deletionStart;	// Location of the start of the deletion in Mbp
    double deletionEnd;		// Location of the end of the deletion in Mbp
    double inversionStart;	// Location of the start of the inversion in Mbp
    double inversionEnd;	// Location of the end of the inversion in Mbp
    int remainingStrandID;	// New strand containing the three fragments, 2 original and one inverted.
    int excisedStrandID;	// New strand containing the deleted section, one fragment.
};





// Represents a combined detected deletion-translocation event.
// Two strands contain two translocated fragments, and another strand
// contains a deleted segment of one of the original strand IDs, for 
// a total of three new data records per delTra mutation
struct SDRdeletionTranslocationEvent
{
    int oldStrandAid;				// Old strand ID for strand A involved in the translocation
    int oldStrandBid;				// Old strand ID for strand B involved in the translocation
    int deletedOldStrandID;   			// whichever of oldStrandA/oldStrandB has the gap, the
    std::vector<double> cleanBreakPositions;	// One or more breakpoints on the non-deleted strand
    double deletionStart;			// Location of where the deletion started in Mbp
    double deletionEnd;				// Location of where the deletion ended in Mbp
    int newStrandID1;         			// First translocation record ID
    int newStrandID2;         			// Second translocation record ID
    int excisedStrandID;      			// The excised strand's ID

};




// One strand contains three fragments, including one foreign fragment,
// and the other strand contains two fragments with a missing segment
// in the middle, so two new records for Del-Ins
struct SDRdeletionInsertionEvent
{
    int donorOldStrandID;            // Old strand that lost material
    int recipientOldStrandID;        // Old strand that received the inserted material
    double segmentStart;             // Start position (on donorOldStrandID) of the deleted/inserted segment
    double segmentEnd;               // End position (on donorOldStrandID) of the deleted/inserted segment
    double insertionPoint;           // Position (on recipientOldStrandID) where the segment was inserted
    int donorRemainingNewStrandID;   // newStrandID of the donor's post-deletion "remaining" record
    int recipientNewStrandID;        // newStrandID of the recipient's record containing the insertion
};






// Represents a detected chromoplexy event: a connected cluster of 3+
// original strands linked by inter-strand rearrangements (balanced
// translocations, deletion-translocations, deletion-insertions).
// Plain intra-strand deletions on component strands are reported as
// accompanying context, not required for classification.
struct SDRchromoplexyEvent
{
    std::vector<int> involvedStrandIDs;               // The 3+ old strand IDs in this connected component, sorted
    std::vector<int> contributingNewStrandIDs;        // newStrandIDs of every record that connected two or more of these strands
    std::vector<int> strandsWithAccompanyingDeletions;// Component strands that also show a plain long deletion (context only)
};








// Represents a karyogram-drawing decision for a specific chromoplexy shape:
// one original strand (consumedStrandID) had BOTH of its pieces fused onto
// two OTHER strands' centromere-bearing fragments (one piece per partner),
// so no data record is left representing it on its own, it should be
// hidden from the karyogram entirely rather than drawn as an untouched
// baseline chromosome. The two partner strands' own remaining acentric
// fragments (the halves that did NOT fuse with the consumed strand) are
// the "loose ends" of the same 3-way rearrangement, and get combined and
// drawn together in their own dedicated karyogram slot after X and Y.
struct SDRchromoplexyAcentricPlacement
{
    int consumedStrandID;                     // Original strand to omit entirely from the karyogram
    std::vector<int> leftoverNewStrandIDs;    // newStrandID(s) whose fragments make up the new combined slot either
                                              // one record that already combines both partners' leftover pieces, or
                                              // two separate single-fragment records, one per partner.
};








// Represents a detected chromothripsis event: a cluster of 1 or 2
// original strands with an unusually high number of DETECTED mutation
// events touching them. Greater than 10 fragments across all involved 
// chromosomes = chromothripsis.
struct SDRchromothripsisEvent
{
    std::vector<int> involvedStrandIDs;         // The 1-2 old strand IDs in this cluster
    std::vector<int> contributingNewStrandIDs;  // newStrandIDs of every rearranged record touching these strand(s)
    int totalFragmentCount;                     // Total data-field-3 fragments summed across all contributing records
};

// A complex genomic rearrangement is found on the same 1-2 strand clusters as chromothripsis (see detectComplexGenomicRearrangements()),
// so it carries the same information.
using SDRcomplexRearrangementEvent = SDRchromothripsisEvent;

// Represents a detected indel: a small insertion and/or deletion of < 100 bp at one site of one record. Indels are
// counted only, never drawn. A small deletion is a gap of < 100 bp between two consecutive pieces of the same strand
// (e.g. 1/0/50.000000/0, 1/50.000005/132.6/1), with no excised record. A small insertion is a run of one-base pieces
// copied from anywhere in the genome (start == end, e.g. 5/73.187925/73.187925/0) sitting between two pieces of the
// receiving strand. A deletion and an insertion at the same site are one indel. Like a long deletion, the deleted bases
// may also be listed as their own excised record (e.g. 1/50.000001/50.000004/0), linear or circular, in either orientation
// (a < 100 bp ecDNA cut from the gap is scored as this indel, not as ecDNA, and is not drawn), or be absent.
struct SDRindelEvent
{
    int oldStrandID;		// Strand that received the indel (the strand of the piece left of the site), never the source of inserted bases
    int newStrandID;		// Record containing the indel
    double position;		// Location (in Mbp) of the site: the end of the piece left of it
    int deletedBases;		// Bases missing between the two flanking pieces (0 for a pure insertion)
    int insertedBases;		// Inserted bases (0 for a pure deletion); a piece with start == end is one base
    int excisedStrandID;	// Record holding the deleted bases when the file lists them (same as a long deletion's excised record), -1 when there is none
};

} // namespace karyoke







// Represents a detected ring chromosome: a circular record that still
// carries a centromere, distinct from ecDNA (circular, but WITHOUT a
// centromere). Can be made up of fragments from more than one original
// strand (a fused, multi-chromosome ring).
struct SDRchromosomeRingEvent
{
    int newStrandID;                        // The new strand ID of this circular, centromere-bearing record
    std::vector<int> involvedOldStrandIDs;  // Every distinct original strand contributing material, in first-appearance order
    int centromereOldStrandID;              // Which original strand's centromere this ring retained
};








// Represents a detected dicentric chromosome: a single linear record
// carrying two centromeres from two different original strands.
struct SDRdicentricEvent
{
    int newStrandID;             // The new strand ID of this dicentric record
    int centromereStrandAid;     // Original strand ID of the FIRST centromere-bearing fragment
    int centromereStrandBid;     // Original strand ID of the SECOND centromere-bearing fragment
    std::vector<int> involvedOldStrandIDs;  // Every distinct original strand contributing material (centromere-bearing or not), in first-appearance order


};


#endif
