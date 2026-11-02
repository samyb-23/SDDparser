#ifndef KARYOGRAM_H
#define KARYOGRAM_H

#include <cairo/cairo.h>

#include <string>
#include <vector>
#include <utility>

#include "SDDtypes.h"
#include "SDRtypes.h"

namespace sddparser
{

struct RGB                                      // Use to define chromosome colors using Red, Green, and Blue codes.
{
    double r;
    double g;
    double b;
};


struct PaintedSegment
{
    double startFraction;                       // Position along new strand as fraction of total length
    double endFraction;
    RGB color;                                  // Color of original chromosome this segment came from
    bool hasCentromere;                         // Does the centromere exist in this segment
    bool isReversed;                            // Check for inversions to draw chevron inversion marker.
    bool isSecondHomolog;                       // Fragment came from the second homolog: its stretch of chromosome outline is dashed. White gaps inherit from a neighbour.
    double centromereStartFraction;
    double centromereEndFraction;
};


enum class ChromosomeLayout                     // Determine how the user passed the chromosome sizes to modify karyogram plotting logic.
{
    NON_HOMOLOGOUS,       			// 1,2,3,...,22,Y,X
    ADJACENT_HOMOLOGS,    			// 1,1,2,2,...,22,22,Y,X
    SPLIT_HOMOLOGS        			// 1,2,...,22,1,2,...,22,Y,X
};


enum class KaryogramFormat                      // Image format for karyogram output (files or in-memory).
{
    PNG,                                        // Raster image, identical to the PNG files written to disk.
    SVG                                         // Vector image: stays sharp at any zoom level, used by the web viewer.
};


struct CentromerePosition                       // Map the centromere start and end positions to the corresponding chromosome for the human genome alone.
{
    std::string chromosome;                     // 1-22, X, Y
    long long start;                            // Start position in base pairs
    long long end;                              // End position in base pairs

};




struct KaryogramLabelRegion                     // Where a chromosome label was drawn, so the web viewer can make it clickable.
{
    std::string label;                          // Text drawn under the chromosome, e.g. "Chr 5", "X" or "Y".
    double x;                                   // Bounding box of the label, in image pixels (same units as the PNG/SVG).
    double y;
    double width;
    double height;
    std::vector<int> oldStrandIDs;              // Chromosome(s) drawn in this slot: SDR original strand IDs or SDD data Field 3 chromosome numbers.
                                                // Homolog 1 first, then homolog 2 for homologous layouts.
};


// Use Karyogram for all illustrations of chromosome damages
class Karyogram
{
public:

    bool generateSDDkaryogram(				// Function that generates the overall Karyogram structure and calls individual functions to draw each chromosome based on sizes.
        const std::vector<double>& chromosomeSizes,	// Use to scale the chromosome heights and convert to number of pixels on the Karyogram.
	const std::vector<double>& cellCyclePhase, 	// Look at SDD header 'Cell cycle phase' to determine if chromosomes should be replicated or not in the final karyogram.
	const std::vector<double>& doseOrFluence,	// Look at SDD header for dose or fluence to be summarized in Karyogram.
	const std::vector<int>& incidentParticles,	// Look at SDD header for 'Incident Particles' to be summarized in karyogram.
	const std::vector<Exposure>& exposures,		// Loop through Exposures object to label chromosome IDs onto the karyogram
	bool humanGenome,				// Check if the user passed -kh (human) or -ko (other) on the command line. if humanGenome = true, adjust the centromere locations. If humanGenome = false, use generic centromere locations.
        const std::string& outputFilename,		// Specify generic output file name combining the input file prefix with the suffix "_karyogram.png"
        KaryogramFormat format = KaryogramFormat::PNG	// File format: PNG by default, or SVG (vector image, e.g. for "_karyogram.svg").
    );

    bool generateSDRkaryogram(				// Function to generate the karyogram of the rearrangements in an SDR file
	const SDRmasterHeader& masterHeader,		// One master header per SDR file containing most importantly the sizes of the original chromosomes before mutations
	const SDRsubHeader& subHeader,			// One subheader per cell in the SDR file, iterate through each subheader for each cell
	bool humanGenome,				// Check if user specified human genome to adjust centromere positions on karyogram
	const std::string& outputFilename,		// Output file names for the summary file and the karyogram drawing.
	KaryogramFormat format = KaryogramFormat::PNG	// File format: PNG by default, or SVG (vector image).
    );

    // In-memory overloads: identical drawing, but the image is returned in imageBytes instead of being written to disk (used by the web viewer).
    // The format defaults to PNG; pass KaryogramFormat::SVG for a vector image that stays sharp when zoomed in.
    bool generateSDDkaryogram(
        const std::vector<double>& chromosomeSizes,
	const std::vector<double>& cellCyclePhase,
	const std::vector<double>& doseOrFluence,
	const std::vector<int>& incidentParticles,
	const std::vector<Exposure>& exposures,
	bool humanGenome,
        std::vector<unsigned char>& imageBytes,		// Receives the encoded PNG or SVG image.
        KaryogramFormat format = KaryogramFormat::PNG
    );

    bool generateSDRkaryogram(
	const SDRmasterHeader& masterHeader,
	const SDRsubHeader& subHeader,
	bool humanGenome,
	std::vector<unsigned char>& imageBytes,		// Receives the encoded PNG or SVG image.
	KaryogramFormat format = KaryogramFormat::PNG
    );


    // Chromosome label positions from the most recent generateSDRkaryogram() call (empty after an SDD karyogram).
    // Used by the CLI to make the labels clickable in the web viewer.
    const std::vector<KaryogramLabelRegion>& getLabelRegions() const { return labelRegions; }


private:

    std::vector<KaryogramLabelRegion> labelRegions;	// Filled while drawing an SDR karyogram, see getLabelRegions().

    // Draws a chromosome label at the current point and records its bounding box and original strands in labelRegions.
    void drawChromosomeLabel(cairo_t* cr, const std::string& label, double x, double y, const std::vector<int>& oldStrandIDs);

    // Shared drawing code behind both overloads. Writes to outputFilename when imageBytes is null, otherwise encodes the image into imageBytes, in the given format.
    bool renderSDDkaryogram(const std::vector<double>& chromosomeSizes, const std::vector<double>& cellCyclePhase, const std::vector<double>& doseOrFluence, const std::vector<int>& incidentParticles, const std::vector<Exposure>& exposures, bool humanGenome, const std::string& outputFilename, std::vector<unsigned char>* imageBytes, KaryogramFormat format);
    bool renderSDRkaryogram(const SDRmasterHeader& masterHeader, const SDRsubHeader& subHeader, bool humanGenome, const std::string& outputFilename, std::vector<unsigned char>* imageBytes, KaryogramFormat format);
    // Creates the drawing surface: an image (pixel) surface for PNG, or an SVG surface that writes vector output to outputFilename or streams it into imageBytes.
    cairo_surface_t* createKaryogramSurface(int width, int height, const std::string& outputFilename, std::vector<unsigned char>* imageBytes, KaryogramFormat format);
    // Finishes the image: encodes the PNG to the file or buffer, or completes the SVG document.
    bool writeKaryogramImage(cairo_surface_t* surface, const std::string& outputFilename, std::vector<unsigned char>* imageBytes, KaryogramFormat format);
    // Overload for a karyogram drawn on a recording surface (everything except the legend): places the legend a fixed gap below the lowest
    // thing drawn (the labels of the lowest chromosome row), sizes the image to fit, then writes it as above. Takes ownership of `content`.
    bool writeKaryogramImage(cairo_surface_t* content, int width, void (Karyogram::*drawLegend)(cairo_t*, double), const std::string& outputFilename, std::vector<unsigned char>* imageBytes, KaryogramFormat format);

    // Function to generate chromosome colors depending on the number of chromosomes, and to choose colors that are different enough between successive chromosomes.
    RGB generateChromosomeColor(int chromosomeNumber, int totalChromosomes);
    // For SDR tracking of chromosome colors during rearrangements
    RGB getColorForOriginalStrand(int oldStrandID, const SDRmasterHeader& masterHeader);
    bool isSecondHomologStrand(int oldStrandID, const SDRmasterHeader& masterHeader);

    // Helper functions to build the mutated segments onto the original chromosome strands, adjusting the color depending on the origin of the mutation, or white for deletions
    std::vector<PaintedSegment> buildPaintedSegments(const SDRdataRecord& record, bool humanGenome, const SDRmasterHeader& masterHeader);
    std::vector<PaintedSegment> buildDeletionRemainingSegments(const SDRdataRecord& record, int homeOldStrandID, bool humanGenome, const SDRmasterHeader& masterHeader);
    std::vector<PaintedSegment> buildMixedStrandSegmentsWithGaps(const SDRdataRecord& record, int homeOldStrandID, bool humanGenome, const SDRmasterHeader& masterHeader, double& outTotalLengthMbp);

    // The main draw chromosome function, accounting for individual chromosome sizes and centromere ranges and locations.
    void drawChromosome(cairo_t* cr, double x, double y, double height, double width, RGB color, double centromereStart, double centromereEnd, bool dashedOutline = false);
    // The main draw mutated chromosome function, accounting additionally for foreign segments with different colors, inversion markers, white deletion segments.
    void drawPaintedFragment(cairo_t* cr, double x, double y, double height, double width, const std::vector<PaintedSegment>& segments, bool roundTopCap = true, bool roundBottomCap = true, bool drawOutline = true);
    std::vector<PaintedSegment> applyMinimumSegmentHeight(const std::vector<PaintedSegment>& segments, double height, double minPixels);	// Keeps tiny segments at least minPixels tall
    // Function to help stack mutated segments onto original intact chromosomes with their original colors, determines positioning of the mutated segments and deleted segments
    void drawStackedMutations(cairo_t* cr, const std::vector<const SDRdataRecord*>& records, double slotCenterX, double posY, double chromosomeWidth, double maxLengthMbp, double maxRenderHeight, bool humanGenome, const SDRmasterHeader& masterHeader, int homeOldStrandID);

    // Use to obtain the stored single- and double-strand break locations in each exposure and determine their coordinates on the Karyogram.
    std::vector<DamageLocation> getDoubleStrandBreaks(const std::vector<Exposure>& exposures);
    std::vector<DamageLocation> getSingleStrandBreaks(const std::vector<Exposure>& exposures);

    // Convert damage locations in base pairs to a fractional length along the chromosome for easy Karyogram pixel conversion.
    double getDamageFraction(const DamageLocation& damage, double chromosomeSize);

    // Functions to draw the damage locations at a given x and y coordinate on the Karyogram, the entire width of the drawn chromosome.
    void drawDoubleStrandBreakMarker(cairo_t* cr, double x, double y, double chromosomeWidth);
    void drawSingleStrandBreakMarker(cairo_t* cr, double x, double y, double markerLength);
    // Draw symbol for inverted segments arising from mutations in the SDR file
    void drawInversionChevron(cairo_t* cr, double centerX, double topY, double bottomY, double width);
    // Draw symbol for ecDNA mutated segments found in SDR file.
    void drawCircularFragment(cairo_t* cr, double centerX, double centerY, double diameter, RGB color);
    // Draw ecDNA mutated segment combined from multiple foreign strands in SDR file
    void drawPieChartFragment(cairo_t* cr, double centerX, double centerY, double diameter, const std::vector<std::pair<RGB, double>>& wedges);
    void drawRingFragment(cairo_t* cr, double centerX, double centerY, double outerDiameter, RGB color, double centromereStartFraction, double centromereEndFraction, bool dashedOutline = false);
    // Draw a circular piece of a single chromosome: a ring (with centromere) if any of its fragments carries a centromere, otherwise ecDNA.
    void drawSingleOriginCircularPiece(cairo_t* cr, const SDRdataRecord& record, double centerX, double centerY, double diameter, int homeOldStrandID, bool humanGenome, const SDRmasterHeader& masterHeader);
    // Draw the inner circle on ecDNA fragments to resemble a ring fragment.
    void drawCircularFragmentInnerRing(cairo_t* cr, double centerX, double centerY, double outerDiameter);

    // Returns a chromosome's corresponding centromere start and end locations to draw the centromere ellipse on the karyogram.
    const CentromerePosition* getHumanCentromere(int chromosomeID);
    bool getCentromereForOriginalStrand(int oldStrandID, bool humanGenome, const SDRmasterHeader& masterHeader, double&centromereStartBP, double& centromereEndBP);			// Check if original strand has centromere in the given range

    // Helper functions to modify the chromosome heights and widths to adjust for mutated content being carried over into this slot.
    double computeSDRbarHeight(double lengthMbp, double maxLengthMbp, double maxRenderHeight);
    double computeMaxBarHeight(const std::vector<const SDRdataRecord*>& records, double maxLengthMbp, double maxRenderHeight, int homeOldStrandID, const SDRmasterHeader& masterHeader);
    double computeSlotWidth(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, double chromosomeWidth, double maxLengthMbp, double maxRenderHeight, const SDRmasterHeader& masterHeader, double& outColumn2Width);
    int computeExcisedColumnLayout(const std::vector<const SDRdataRecord*>& excisedRecords, double maxLengthMbp, double maxRenderHeight, double maxColumnHeight, double verticalGap, std::vector<double>& outHeights);

    // Function to draw the SDD karyogram legend at the bottom.
    void drawSDDlegend(cairo_t* cr, double legendY);
    // Function to draw the karyogram summary box at the top.
    void drawSDDsummary(cairo_t* cr, const std::vector<double>& cellCyclePhase, const std::vector<Exposure>& exposures, const std::vector<double>& doseOrFluence, const std::vector<int>& incidentParticles);
    // Function to draw the top SDR karyogram summary box summarizing which mutations were in the SDR file
    void drawSDRsummary(cairo_t* cr, const SDRmasterHeader& masterHeader, const SDRsubHeader& subHeader);
    // Draw legend at bottom of the SDR karyogram which has different symbols than the SDD file karyogram
    void drawSDRlegend(cairo_t* cr, double legendY);

    // Function that will determine how the user passed the sizes of the chromosomes (non-homologous, adjacent-homologous, or split-homologous order).
    ChromosomeLayout determineChromosomeLayout(const std::vector<double>& chromosomeSizes);

    // Function to check if the SDR data entries are depciting mutations, or are the original intact strand information
    std::vector<const SDRdataRecord*> filterBaselineIfMutated(const std::vector<const SDRdataRecord*>& records, int numOriginalStrands);

    // Functions to alter the drawing of certain mutations depending on the specific mutation involved
    bool isDeletionShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, const SDRmasterHeader& masterHeader);				// For long deletions
    bool isECDNAshape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, const SDRmasterHeader& masterHeader);				// For ecDNA mutations
    bool isDeletionInversionShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID);							// For deletion-inversion mutations
    bool isDeletionTranslocationDonorShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID);						// For deletion-translocation mutations
    bool analyzeHomeStrandGapShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, const SDRmasterHeader& masterHeader, const SDRdataRecord*& remainingRecord, std::vector<const SDRdataRecord*>& pieceRecords);	// Unified gap shape: deletions (with or without excised), del-inversions, ecDNA, lone gaps
    bool isChromosomeRingShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, const SDRmasterHeader& masterHeader);
    bool isLoneGapShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID);									// For deletion-insertion mutations

    // Helper functions to stack mutated chromosome segments of different colors on top of the original segments and remove the correct amount of 
    // chromosome if a deletion occurs.
    SDRdataRecord concatenateRecordsForDrawing(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID);
    std::vector<const SDRdataRecord*> clusterRecordsForDrawing(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, std::vector<SDRdataRecord>& mergedStorage, const SDRmasterHeader& masterHeader);

    // If intact chromosome strands are not listed in SDR data record, synthesize them using the intact chromosome sizes.
    std::vector<const SDRdataRecord*> synthesizeIntactRecordIfMissing(int oldStrandID, int cellID, const SDRmasterHeader& masterHeader, std::vector<SDRdataRecord>& intactRecordStorage);

};

} //namespace sddparser


#endif
