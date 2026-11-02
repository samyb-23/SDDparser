#include <cairo/cairo.h>
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>

#include "Karyogram.h"
#include "SDDutilities.h"





namespace sddparser
{

// Bookkeeping struct, remembers each drawn chromosome's pixel
// position and height so the DSB/SSB-drawing code further down can
// look it back up by chromosome/chromatid number.
struct ChromosomeGeometry
{
    int chromosomeNumber;						// Chromosome Number (i.e. ID).
    int chromatidNumber;						// Chromatid Number (1 for unduplicated, 1 or 2 for duplicated)
    double x;								// X pixel position
    double y;								// Y pixel position
    double height;							// Chromosome height in pixels
};













// Decides whether the karyogram is drawn pre-replication (one chromatid per
// chromosome) or post-replication (two sister chromatids per chromosome),
// combining the header 'Cell cycle phase' with what data Field 3 actually says.
//
//   Phase 3/4/5 (S/G2/M)  -> two chromatids. Damages with chromatid 1 go on the
//                            left sister, chromatid 2 on the right sister.
//   Phase 1/2   (G0/G1)   -> one chromatid, unless Field 3 references
//                            chromatid 2, which is inconsistent with the header;
//                            warn and draw two chromatids so no damage is lost.
//   Phase 0 / other       -> unspecified, so infer it from Field 3: any
//                            chromatid 2 means the genome was replicated.
bool resolveDoubleChromatid(
    const std::vector<double>& cellCyclePhase,				// SDD header 'Cell cycle phase'
    const std::vector<Exposure>& exposures)				// Damages, to read the chromatid numbers in data Field 3
{
    const int phase = cellCyclePhase.empty() ? 0 : static_cast<int>(cellCyclePhase[0]);

    int secondChromatidDamages = 0;					// Damage sites placed on chromatid 2 by Field 3

    for (const auto& exposure : exposures)
    {
        for (const auto& damage : exposure.damages)
        {
            if (damage.chromosomeID.chromatidNumber == 2)
            {
                ++secondChromatidDamages;
            }
        }
    }

    const bool dataIsReplicated = secondChromatidDamages > 0;

    if (phase == 3 || phase == 4 || phase == 5)				// Header says post-replication
    {
        if (!dataIsReplicated)
        {
            std::cerr << "NOTE: Cell cycle phase is post-replication but no damage in data Field 3 "
                      << "references chromatid 2; all damages are drawn on the first sister chromatid.\n";
        }
        return true;
    }

    if (phase == 1 || phase == 2)					// Header says pre-replication
    {
        if (dataIsReplicated)
        {
            std::cerr << "WARNING: Cell cycle phase is G0/G1 (pre-replication) but "
                      << secondChromatidDamages << " damage site(s) reference chromatid 2 in data Field 3. "
                      << "Drawing two sister chromatids so these damages are not dropped.\n";
        }
        return dataIsReplicated;
    }

    // Phase unspecified (0) or unknown: let data Field 3 decide.
    if (dataIsReplicated)
    {
        std::cerr << "NOTE: Cell cycle phase is unspecified; data Field 3 references chromatid 2, "
                  << "so the karyogram is drawn post-replication.\n";
    }
    return dataIsReplicated;
}












bool Karyogram::generateSDDkaryogram(                                         // Writes the karyogram to a PNG (default) or SVG file.
    const std::vector<double>& chromosomeSizes,
    const std::vector<double>& cellCyclePhase,
    const std::vector<double>& doseOrFluence,
    const std::vector<int>& incidentParticles,
    const std::vector<Exposure>& exposures,
    bool humanGenome,
    const std::string& outputFilename,
    KaryogramFormat format)
{
    return renderSDDkaryogram(chromosomeSizes, cellCyclePhase, doseOrFluence, incidentParticles, exposures, humanGenome, outputFilename, nullptr, format);
}


bool Karyogram::generateSDDkaryogram(                                         // Returns the karyogram as PNG or SVG bytes in memory.
    const std::vector<double>& chromosomeSizes,
    const std::vector<double>& cellCyclePhase,
    const std::vector<double>& doseOrFluence,
    const std::vector<int>& incidentParticles,
    const std::vector<Exposure>& exposures,
    bool humanGenome,
    std::vector<unsigned char>& imageBytes,
    KaryogramFormat format)
{
    return renderSDDkaryogram(chromosomeSizes, cellCyclePhase, doseOrFluence, incidentParticles, exposures, humanGenome, std::string(), &imageBytes, format);
}


bool Karyogram::renderSDDkaryogram(                                           // Function that generates the overall Karyogram structure
    const std::vector<double>& chromosomeSizes,                            // Requires chromosome sizes in the SDD file header for scaling on the image
    const std::vector<double>& cellCyclePhase,                             // Requires cell cycle phase in the SDD file header for determining presence of duplicated chromosomes.
    const std::vector<double>& doseOrFluence,                              // Accesses the SDD dose or fluence listed in the header
    const std::vector<int>& incidentParticles,                             // Accesses the SDD incident particles in the header
    const std::vector<Exposure>& exposures,                                // Accesses the damage locations associated with each chromosome in a given exposure.
    bool humanGenome,                                                      // Boolean to determine if karyogram should be drawn with human centromere ranges or generic centromere locations for non-human genomes.
    const std::string& outputFilename,                                     // Outputs the image to the outputFilename (when imageBytes is null).
    std::vector<unsigned char>* imageBytes,                                // If non-null, receives the image in memory instead.
    KaryogramFormat format)                                                // PNG or SVG.
{

    labelRegions.clear();                                                  // Filled again below by drawChromosomeLabel, one region per chromosome label.

    if (chromosomeSizes.empty())                                           // Make sure SDD header 'Chromosome sizes' Is non-empty.
    {
        std::cerr << "ERROR: No chromosome sizes provided.\n";
        return false;
    }

    std::vector<DamageLocation> dsbVec = getDoubleStrandBreaks(exposures); // Use this vector for plotting DSBs onto karyogram

    std::vector<DamageLocation> ssbVec = getSingleStrandBreaks(exposures);

    std::vector<ChromosomeGeometry> chromosomeGeometry;                    // Store the geometry of the chromosomes in x and y pixel coordinates to draw the DSBs.

    const int chromosomeCount = static_cast<int>(chromosomeSizes[0]);      // The first entry specifies the number of chromosomes.

    //const int homologousSetSize = (chromosomeCount - 2) / 2;             // Store number of homologs for DSB drawing

    // Generic metacentric chromosomes similar to MEDRAS-MC
    const double genericCentromereStart = 0.49;                            // In case non-human chromosomes specified, draw generic centromere regions
    const double genericCentromereEnd = 0.51;                              // Will be replaced later if human genome is specified with actual centromere positions.

    if (chromosomeSizes.size() !=                                          // Check that the vector contains the expected number of chromosome sizes.
        static_cast<size_t>(chromosomeCount + 1))
    {
        std::cerr << "ERROR: Chromosome count does not match "
                  << "the number of chromosome sizes provided.\n";
        return false;
    }

    const bool doubleChromatid = resolveDoubleChromatid(cellCyclePhase, exposures);                                          // Determine whether one or two chromatids per chromosome should be drawn

    // ------------------------------------------------------------
    // Determine chromosome layout
    // ------------------------------------------------------------

    const ChromosomeLayout chromosomeLayout = determineChromosomeLayout(chromosomeSizes); // Determine which layout the chromosome sizes are listed in for plotting logic.

    const bool hasHomologs = chromosomeLayout != ChromosomeLayout::NON_HOMOLOGOUS;  // Default layout is non-homologous 1,...,N,Y,X

    const int homologousPairs = hasHomologs ? (chromosomeCount - 2) / 2 : 0;   // Number of homologous pairs if homologs present.

    const int drawableGroups = hasHomologs ? homologousPairs : chromosomeCount - 2; // Number of chromosome groups that will be drawn before the Y and X

    const int totalColorGroups = drawableGroups + 2;                       // Use totalColorGroups to have evenly distrbuted colors across all chromosomes.

    double maxBP = 0.0;                                                    // Find largest chromosome, scale all chromosomes with respect to the largest

    for (size_t i = 1 ; i < chromosomeSizes.size(); i++)
    {
        if (chromosomeSizes[i] > maxBP)
        {
            maxBP = chromosomeSizes[i];
        }
    }

    if (maxBP <= 0.0)
    {
        std::cerr << "ERROR: Chromosome sizes must contain positive values.\n";
        return false;
    }

    // ----------------------------------------------------
    // Karyogram dimensions
    // ----------------------------------------------------

    const int imgWidth = 1000;                                             // Image width in number of pixels
    const int columns = 4;                                                 // Desired number of columns in the karyogram
    const double colWidth = static_cast<double>(imgWidth) / columns;       // Column width also in number of pixels
    const double rowHeight = 250.0;                                        // Row height in number of pixels
    const double startY = 188.0;                                           // Choosing the starting Y position for the first row of chromosomes
    const double maxRenderHeight = 180.0;                                  // Limit the maximum height of the chromosome, also in pixels
    const double chromosomeWidth = 20.0;                                   // Individual chromosomes are 20 pixels wide.
    const double chromatidGap = 1.5;                                       // The gaps between homologs in a replicated chromosome are 1.5 pixels.


    // ----------------------------------------------------
    // Generate the background Karyogram template
    // ----------------------------------------------------
    // Everything except the legend is drawn on a recording surface first: the image height is only known once the lowest label
    // has been drawn, since the legend starts a fixed gap below it (see writeKaryogramImage).
    cairo_surface_t* content = cairo_recording_surface_create(CAIRO_CONTENT_COLOR_ALPHA, nullptr);

    cairo_t* cr = cairo_create(content);

    cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);

    // ----------------------------------------------------
    // Draw SDD summary
    // ----------------------------------------------------

    drawSDDsummary(cr, cellCyclePhase, exposures, doseOrFluence, incidentParticles);


    // -----------------------------------------------------------
    // Draw chromosomes depending on chromosome size layout given
    // -----------------------------------------------------------

    for (int i = 0; i < drawableGroups; ++i)
    {
        size_t firstIndex = 0;                                             // Declare index of first homolog to be determined depending on chromosome sizes layout.
        size_t secondIndex = 0;                                            // Second homolog begins after the first

        // -----------------------------------------------------
        // Determine indexing based off chromosome sizes layout
        // -----------------------------------------------------

        if (chromosomeLayout == ChromosomeLayout::SPLIT_HOMOLOGS)          // Chromosomes listed [count,1,2,3,..., 1,2,3,...,Y,X]
        {
            firstIndex = static_cast<size_t>(i + 1);
            secondIndex = static_cast<size_t>(i + 1 + homologousPairs);
        }
        else if (chromosomeLayout == ChromosomeLayout::ADJACENT_HOMOLOGS)  // Chromosomes listed [count,1,1,2,2,...,Y,X]
        {
            firstIndex = static_cast<size_t>(2 * i + 1);
            secondIndex = static_cast<size_t>(2 * i + 2);
        }
        else                                                               // Chromosomes listed [count,1,2,3,4...,N-2,Y,X]
        {
            firstIndex = static_cast<size_t>(i + 1);
        }


        // Determine position of chromosome group
        int col = i % columns;
        int row = i / columns;

        double groupCenterX = (col * colWidth) + (colWidth / 2.0);
        double posY = startY + (row * rowHeight);

        RGB chromosomeColor = generateChromosomeColor(i, totalColorGroups);// Homologs receive same color, non-homologous chromosomes receive their own color.


        // -----------------------------------------------------------------
        // NON-HOMOLOGOUS CHROMOSOME DRAWING
        // -----------------------------------------------------------------

        if (chromosomeLayout == ChromosomeLayout::NON_HOMOLOGOUS)
        {
            const double chromosomeSize = chromosomeSizes[firstIndex];     // Only care about firstIndex because there are no homologs.

            double chromosomeHeight = (chromosomeSize / maxBP) * maxRenderHeight; // Normalize chromosome heights relative to the largest chromosoem maxBP

            if (chromosomeHeight < 30.0)                                   // Prevent too small chromosome sizes.
            {
                chromosomeHeight = 30.0;
            }

            // If non-human genome passed, use generic centromere start and end locations
            double centromereStart = genericCentromereStart;
            double centromereEnd = genericCentromereEnd;

            // Adjust centromere locations to the ones listed at the start if a humanGenome specified.
            if (humanGenome)
            {
                int humanChromosomeID = static_cast<int>(firstIndex);

                if (chromosomeCount == 24)
                {
                    if (firstIndex == 23)                                  // When chromosome sizes index 23 is reached, associate that with the Y chromosome
                    {
                        humanChromosomeID = 45;   // Y
                    }
                    else if (firstIndex == 24)                             // When chromosome sizes index 24 is reached, associate that with the X chromosome.
                    {
                        humanChromosomeID = 46;   // X
                    }
                }


                const CentromerePosition* centromere = getHumanCentromere(humanChromosomeID);       // Get corresponding centromere positions per chromosome ID

                if (centromere != nullptr)
                {
                    // Set corresponding centromere ranges and convert to fraction of the chromosome length for conversion to pixels later.
                    centromereStart = static_cast<double>(centromere->start) / (chromosomeSize * 1000000.0);        // chromosome sizes listed in Mbp, convert to bp.
                    centromereEnd = static_cast<double>(centromere->end)/ (chromosomeSize*1000000.0);
                }
            }


            // ------------------------------------------------
            // Single-chromatid drawing
            // ------------------------------------------------
            if (!doubleChromatid)
            {
                const double posX = groupCenterX - chromosomeWidth / 2.0;  // Set the horizontal position of the single chromatid drawing

                // Draw single-chromatid chromosome and store the ID, chromatid number, position and size in chromosomeGeometry.
                drawChromosome(cr, posX, posY, chromosomeHeight, chromosomeWidth, chromosomeColor, centromereStart, centromereEnd);
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 1, posX, posY, chromosomeHeight});

            }

            // ------------------------------------------------
            // Two chromatids drawing
            // ------------------------------------------------
            else
            {
                // Determine horizontal X positions in the karyogram of the two chromatids for the corresponding chromosome.
                const double posX = groupCenterX - chromosomeWidth - chromatidGap / 2.0;
                const double posX2 = groupCenterX + chromatidGap / 2.0;

                // Draw both chromatids of the chromosome in their respective X positions.
                drawChromosome(cr, posX, posY, chromosomeHeight, chromosomeWidth, chromosomeColor, centromereStart, centromereEnd);
                drawChromosome(cr, posX2, posY, chromosomeHeight, chromosomeWidth, chromosomeColor, centromereStart, centromereEnd);

                // Both chromatids in the non-homologous case are assigned the same chromosome ID, it is just in the post-replicated cell phase, so there are two.
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 1, posX, posY, chromosomeHeight});
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 2, posX2, posY, chromosomeHeight});
            }

            // ------------------------------------------------
            // Label chromosome
            // ------------------------------------------------
            cairo_set_source_rgb(
                cr,
                0.1,
                0.1,
                0.1
            );

            cairo_select_font_face(
                cr,
                "Sans",
                CAIRO_FONT_SLANT_NORMAL,
                CAIRO_FONT_WEIGHT_BOLD
            );

            cairo_set_font_size(cr, 16.0);

            std::string label =
                "Chr " + std::to_string(i + 1);

            // Draws the label and records where it is, so the web viewer can make it clickable (damage summary for this chromosome).
            drawChromosomeLabel(cr, label, groupCenterX - 20.0, posY + chromosomeHeight + 25.0,
                                {static_cast<int>(firstIndex)});
        }

        // -----------------------------------------------------
        // HOMOLOGOUS CHROMOSOMES DRAWING
        // -----------------------------------------------------

        else
        {
            // Initialize centromere positions per chromosome ID for each homolog
            double firstCentromereStart = genericCentromereStart;

            double firstCentromereEnd = genericCentromereEnd;

            double secondCentromereStart = genericCentromereStart;

            double secondCentromereEnd = genericCentromereEnd;

            // If human centromeres, get the positions corresponding to the correct ID.
            if (humanGenome)
            {
                const CentromerePosition* centromere = getHumanCentromere(i + 1);

                if(centromere != nullptr)
                {
                // Obtain centromre start and end positions for both homologs from the CentromerePosition struct.
                    firstCentromereStart = static_cast<double>(centromere->start) / (chromosomeSizes[firstIndex] * 1000000.0);
                    firstCentromereEnd = static_cast<double>(centromere->end) / (chromosomeSizes[firstIndex] * 1000000.0);
                    secondCentromereStart = static_cast<double>(centromere->start) / (chromosomeSizes[secondIndex] * 1000000.0);
                    secondCentromereEnd = static_cast<double>(centromere->end) / (chromosomeSizes[secondIndex] * 1000000.0);
                }

            }


            // ----------------------------------------------
            // Scale chromosome heights
            // ----------------------------------------------

            double firstHeight = (chromosomeSizes[firstIndex] / maxBP) * maxRenderHeight;   // Normalize height of first homolog to the largest chromosome height

            double secondHeight = (chromosomeSizes[secondIndex] / maxBP) * maxRenderHeight; // Normalize height of second homolog to the largest chromosome height.

	    // Clamp minimum chromosome height to 30.0 pixels
            if (firstHeight < 30.0)
            {
                firstHeight = 30.0;
            }

            if (secondHeight < 30.0)
            {
                secondHeight = 30.0;
            }


            // ----------------------------------------------
            // Single-chromatid homologs
            // ----------------------------------------------

            if (!doubleChromatid)
            {
                const double homologGap = 18.0;                            			// Define spacing between homologs in pixels.

                const double leftChromPosX = groupCenterX - chromosomeWidth - homologGap / 2.0;	// Determine the X positions of the left homolog

                const double rightChromPosX = groupCenterX + homologGap / 2.0;              	// Determine the X positions of the right homolog

                // Draw first homolog at the desired location with the desired dimensions
                drawChromosome(cr, leftChromPosX, posY, firstHeight, chromosomeWidth, chromosomeColor, firstCentromereStart, firstCentromereEnd);

                // Draw second homolog at the desired location with the desired dimensions
                drawChromosome(cr, rightChromPosX, posY, secondHeight, chromosomeWidth, chromosomeColor, secondCentromereStart, secondCentromereEnd, true);

                // Store actual SDD chromosome IDs and their corresponding dimensions and locations --> for damage marker drawings later on.
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 1, leftChromPosX, posY, firstHeight});
                chromosomeGeometry.push_back({static_cast<int>(secondIndex), 1, rightChromPosX, posY, secondHeight});
            }


            // ----------------------------------------------
            // Double-chromatid homologs
            // ----------------------------------------------

            else
            {
                const double homologGap = 40.0;                            			// Increase homolog gap to account for second chromatids

                // Determine homolog 1 left and right chromatid positions
                const double leftChromatid1PosX = groupCenterX - homologGap / 2.0 - 2.0 * chromosomeWidth - chromatidGap;
                const double leftChromatid2PosX = leftChromatid1PosX + chromosomeWidth + chromatidGap;

                // Determine homolog 2 left and right chromatid positions
                const double rightChromatid1PosX = groupCenterX + homologGap / 2.0;
                const double rightChromatid2PosX = rightChromatid1PosX + chromosomeWidth + chromatidGap;


                // Draw homolog 1 chromatids at the desired locations with the desired dimensions
                drawChromosome(cr, leftChromatid1PosX, posY, firstHeight, chromosomeWidth, chromosomeColor, firstCentromereStart, firstCentromereEnd);
                drawChromosome(cr, leftChromatid2PosX, posY, firstHeight, chromosomeWidth, chromosomeColor, firstCentromereStart, firstCentromereEnd);


                // Draw homolog 2 chromatids at the desired locations with the desired dimensions
                drawChromosome(cr, rightChromatid1PosX, posY, secondHeight, chromosomeWidth, chromosomeColor, secondCentromereStart, secondCentromereEnd, true);
                drawChromosome(cr, rightChromatid2PosX, posY, secondHeight, chromosomeWidth, chromosomeColor, secondCentromereStart, secondCentromereEnd, true);

                // Store actual SDD chromosome IDs with their corresponding dimensions and locations --> for damage marker drawings later on
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 1, leftChromatid1PosX, posY, firstHeight});
                chromosomeGeometry.push_back({static_cast<int>(firstIndex), 2, leftChromatid2PosX, posY, firstHeight});
                chromosomeGeometry.push_back({static_cast<int>(secondIndex), 1, rightChromatid1PosX, posY, secondHeight});
                chromosomeGeometry.push_back({static_cast<int>(secondIndex), 2, rightChromatid2PosX, posY, secondHeight});
            }


            // ----------------------------------------------
            // Label homologous pair
            // ----------------------------------------------

            const double labelHeight = std::max(firstHeight, secondHeight);

            cairo_set_source_rgb(
                cr,
                0.1,
                0.1,
                0.1
            );

            cairo_select_font_face(
                cr,
                "Sans",
                CAIRO_FONT_SLANT_NORMAL,
                CAIRO_FONT_WEIGHT_BOLD
            );

            cairo_set_font_size(cr, 16.0);

            std::string label =
                "Chr " + std::to_string(i + 1);

            // Draws the label and records where it is, so the web viewer can make it clickable (damage summary per homolog).
            drawChromosomeLabel(cr, label, groupCenterX - 20.0, posY + labelHeight + 25.0,
                                {static_cast<int>(firstIndex), static_cast<int>(secondIndex)});
        }
    }


    // --------------------------------------------------
    // X and Y chromosome drawing
    // --------------------------------------------------

    // Always assume the user gives the Y chromosome then the X chromosome as the last two listed chromosome sizes.
    const size_t xIndex = static_cast<size_t>(chromosomeCount);            // Determining the index of the X chromosome
    const size_t yIndex = static_cast<size_t>(chromosomeCount - 1);        // Determining the index of the Y chromosome

    // Set the centromere locations to the generic positions, replaced if user passed the '-kh' (human) option.
    double xCentromereStart = genericCentromereStart;
    double xCentromereEnd = genericCentromereEnd;
    double yCentromereStart = genericCentromereStart;
    double yCentromereEnd = genericCentromereEnd;

    if (humanGenome)                                                       // If user specified a human genome '-kh'
    {
        int xHumanID = 46;
        int yHumanID = 45;

        // Get the base pair locations of the centromeres for the X and Y chromosomes
        const CentromerePosition* xCentromere = getHumanCentromere(xHumanID);
        const CentromerePosition* yCentromere = getHumanCentromere(yHumanID);

        if (xCentromere != nullptr)
        {
            // Scale centromere locations from Mbp to bp
            xCentromereStart = static_cast<double>(xCentromere->start) / (chromosomeSizes[xIndex] * 1000000.0);
            xCentromereEnd = static_cast<double>(xCentromere->end) / (chromosomeSizes[xIndex] * 1000000.0);
        }

        if (yCentromere != nullptr)
        {
            // Scale centromere locations from Mbp to bp
            yCentromereStart = static_cast<double>(yCentromere->start) / (chromosomeSizes[yIndex] * 1000000.0);
            yCentromereEnd = static_cast<double>(yCentromere->end) / (chromosomeSizes[yIndex] * 1000000.0);
        }
    }

    // Scale the X and Y chromosome heights relative to the largest chromosome.
    double xHeight = (chromosomeSizes[xIndex] / maxBP) * maxRenderHeight;
    double yHeight = (chromosomeSizes[yIndex] / maxBP) * maxRenderHeight;

    // Prevent X and Y chromosomes from being smaller than 30 pixels
    if (xHeight < 30.0)
    {
        xHeight = 30.0;
    }
    if (yHeight < 30.0)
    {
        yHeight = 30.0;
    }

    // -------------------------------------------------------------
    // X and Y chromosome pixel dimensions
    // -------------------------------------------------------------
    // X and Y occupy the remaining positions in the final row.
    const int sexChromosomeRow = (drawableGroups > 0) ? (drawableGroups - 1) / columns : 0;
    // Plot X first then Y next.
    const int xColumn = 2;
    const int yColumn = 3;
    const double xGroupCenterX = (xColumn * colWidth) + (colWidth / 2.0);
    const double yGroupCenterX = (yColumn * colWidth) + (colWidth / 2.0);
    const double sexChromPosY = startY + (sexChromosomeRow * rowHeight);
    const double xChromPosX = xGroupCenterX - chromosomeWidth / 2.0;
    const double yChromPosX = yGroupCenterX - chromosomeWidth / 2.0;


    // Set X and Y chromosome Colors
    RGB yColor = generateChromosomeColor(drawableGroups, totalColorGroups);
    RGB xColor = generateChromosomeColor(drawableGroups + 1, totalColorGroups);

    if (!doubleChromatid)                                                  // Depending on if SDD header 'Cell cycle phase' specified a post-replicated chromosome.
    {
        // Draw X
        drawChromosome(cr, xChromPosX, sexChromPosY, xHeight, chromosomeWidth, xColor, xCentromereStart, xCentromereEnd);
        // Store geometry of X chromosome for damage drawing
        chromosomeGeometry.push_back({chromosomeCount, 1, xChromPosX, sexChromPosY, xHeight});

        // Draw Y
        drawChromosome(cr, yChromPosX, sexChromPosY, yHeight, chromosomeWidth, yColor, yCentromereStart, yCentromereEnd);
        // Store geometry of Y chromosome for damage drawing
        chromosomeGeometry.push_back({chromosomeCount - 1, 1, yChromPosX, sexChromPosY, yHeight});

    }
    else
    {
	// Sister chromatids centred on the X and Y group centres, separated by chromatidGap.
        const double xChromatid1PosX = xGroupCenterX - chromosomeWidth - chromatidGap / 2.0;
        const double xChromatid2PosX = xGroupCenterX + chromatidGap / 2.0;
        const double yChromatid1PosX = yGroupCenterX - chromosomeWidth - chromatidGap / 2.0;
        const double yChromatid2PosX = yGroupCenterX + chromatidGap / 2.0;

        // Draw X chromatids 1 and 2
        drawChromosome(cr, xChromatid1PosX, sexChromPosY, xHeight, chromosomeWidth, xColor, xCentromereStart, xCentromereEnd);
        drawChromosome(cr, xChromatid2PosX, sexChromPosY, xHeight, chromosomeWidth, xColor, xCentromereStart, xCentromereEnd);
        // Store X chromatids 1 and 2 geometries for damage drawing
        chromosomeGeometry.push_back({chromosomeCount, 1, xChromatid1PosX, sexChromPosY, xHeight});
        chromosomeGeometry.push_back({chromosomeCount, 2, xChromatid2PosX, sexChromPosY, xHeight});

        // Draw Y chromatids 1 and 2
        drawChromosome(cr, yChromatid1PosX, sexChromPosY, yHeight, chromosomeWidth, yColor, yCentromereStart, yCentromereEnd);
        drawChromosome(cr, yChromatid2PosX, sexChromPosY, yHeight, chromosomeWidth, yColor, yCentromereStart, yCentromereEnd);

        // Store Y chromatids 1 and 2 geometries for damage drawing
        chromosomeGeometry.push_back({chromosomeCount - 1, 1, yChromatid1PosX, sexChromPosY, yHeight});
        chromosomeGeometry.push_back({chromosomeCount - 1, 2, yChromatid2PosX, sexChromPosY, yHeight});

    }





    // --------------------------------------------------
    // Damage placement using data Field 3
    // --------------------------------------------------
    // Field 3 gives the chromosome number and the chromatid number. Field 4 gives the position along the genetic length (p-arm -> q-arm). Sister chromatids are copies
    // of the same sequence, so the Field 4 position maps identically onto either one; the chromatid number alone decides which sister the marker is drawn on. A damage
    // is drawn on exactly one chromatid, never copied onto both sisters.

    int unplacedDamageSites = 0;                                           // Damage sites whose Field 3 chromosome/chromatid has no drawn chromatid

    // Look up the drawn chromatid matching a damage's Field 3 chromosome and chromatid numbers.
    const auto findChromatidGeometry = [&](const ChromosomeID& id) -> const ChromosomeGeometry*
    {
        for (const auto& geometry : chromosomeGeometry)
        {
            if (geometry.chromosomeNumber == id.chromosomeNumber &&
                geometry.chromatidNumber == id.chromatidNumber)
            {
                return &geometry;
            }
        }
        return nullptr;
    };

    // Convert a damage to marker pixel coordinates on its chromatid. Returns false if it cannot be placed.
    const auto locateDamage = [&](const DamageLocation& damage, double& markerX, double& markerY) -> bool
    {
        const ChromosomeGeometry* geometry = findChromatidGeometry(damage.chromosomeID);

        const size_t sizeIndex = static_cast<size_t>(damage.chromosomeID.chromosomeNumber);

        if (geometry == nullptr || sizeIndex == 0 || sizeIndex >= chromosomeSizes.size())
        {
            ++unplacedDamageSites;
            return false;
        }

       // Fraction of the chromosome length from the p-arm end, clamped to the drawn chromatid.
       const double damageFraction = std::clamp(getDamageFraction(damage, chromosomeSizes[sizeIndex]), 0.0, 1.0);

       markerX = geometry->x + chromosomeWidth / 2.0;
       markerY = geometry->y + damageFraction * geometry->height;
       return true;
    };






    // --------------------------------------------------
    // Draw double-strand break markers
    // --------------------------------------------------

    for (const auto& dsb : dsbVec)                                      // Loop through the double-strand break vector containing exposure data (damages per chromosome per exposure)
    {
       double markerX = 0.0;
       double markerY = 0.0;

        if (locateDamage(dsb, markerX, markerY))
        {
            drawDoubleStrandBreakMarker(cr, markerX, markerY, chromosomeWidth);   // Line spans slightly more than the chromatid width
        }

    }

    // --------------------------------------------------
    // Draw single-strand break markers
    // --------------------------------------------------

    const double ssbBaseLength = 4.0;                                      // Marker length for a site with one SSB
    const double ssbLengthPerBreak = 4.0;                                  // Extra marker length per additional SSB in the site
    const double maxSSBLength = chromosomeWidth;                           // Never longer than the chromatid is wide


    for (const auto& ssb : ssbVec)
    {
	double markerX = 0.0;
        double markerY = 0.0;

	if (!locateDamage(ssb, markerX, markerY))
        {
            continue;
        }

        double markerLength = ssbBaseLength + ssbLengthPerBreak * (ssb.numSingleStrandBreaks - 1);
        markerLength = std::min(markerLength, maxSSBLength);

        drawSingleStrandBreakMarker(cr, markerX, markerY, markerLength);

    }



    if (unplacedDamageSites > 0)                                           // Never drop damages silently: the summary box still counts them.
    {
        std::cerr << "WARNING: " << unplacedDamageSites << " DSB/SSB marker(s) could not be drawn because their "
                  << "data Field 3 chromosome/chromatid number does not match a drawn chromatid "
                  << "(chromatid must be 1" << (doubleChromatid ? " or 2" : "") << ").\n";
    }



    // X label
    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);

    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);

    cairo_set_font_size(cr, 16.0);

    drawChromosomeLabel(cr, "X", xGroupCenterX - 4.0, sexChromPosY + xHeight + 25.0, {static_cast<int>(xIndex)});

    // Y label
    drawChromosomeLabel(cr, "Y", yGroupCenterX - 4.0, sexChromPosY + yHeight + 25.0, {static_cast<int>(yIndex)});

    // ----------------------------------------------------
    // Draw Karyogram legend at bottom, a fixed gap below the lowest label, and write the
    // karyogram to the desired output filename "'SDDfileName'_karyogram.png"
    // ----------------------------------------------------
    cairo_destroy(cr);

    return writeKaryogramImage(content, imgWidth, &Karyogram::drawSDDlegend, outputFilename, imageBytes, format);
}








void Karyogram::drawChromosome(
    cairo_t* cr,
    double x,                                           // Pixel X coordinate
    double y,                                           // Pixel Y coordinate
    double height,                                      // Scaled chromosome pixel height
    double width,                                       // Chromosome width in pixels
    RGB color,                                          // Associated chromosome color
    double centromereStart,                             // Centromere start location as fraction of the chromosome length from p arm to q arm
    double centromereEnd,                               // Centromere end location as fraction of the chromosome length from p arm to q arm.
    bool dashedOutline)                                 // Second homologs get a dashed outline so they can be told apart from the first homolog (same color).
{

    // --------------------------------------
    // Centromere position and dimensions
    // --------------------------------------
    const double capRadius = width / 2.0;
    const double centromereCenter = (centromereStart + centromereEnd) / 2.0;
    const double centromereY = y + (height * centromereCenter);
    const double constrictionHeight = 7.0;
    const double constrictionAmount = 2.0;
    const double constrictionTop = centromereY - constrictionHeight / 2.0;
    const double constrictionBottom = centromereY + constrictionHeight / 2.0;

    cairo_new_path(cr);

    // ---------------------------------------
    // Top Cap
    // ---------------------------------------
    cairo_arc(
        cr,
        x + capRadius,
        y + capRadius,
        capRadius,
        M_PI,
        2 * M_PI
    );


    // -----------------------------------------
    // Right side approaching centromere
    // -----------------------------------------

    cairo_line_to(
        cr,
        x + width,
        constrictionTop
    );

    // Gradually move inward
    cairo_line_to(
        cr,
        x + width - constrictionAmount,
        centromereY
    );

    // Gradually move outward
    cairo_line_to(
        cr,
        x + width,
        constrictionBottom
    );


    // -----------------------------------------
    // Right side below centromere
    // -----------------------------------------

    cairo_line_to(
        cr,
        x + width,
        y + height - capRadius
    );

    // -----------------------------------------
    // Bottom cap
    // -----------------------------------------

    cairo_arc(
        cr,
        x + capRadius,
        y + height - capRadius,
        capRadius,
        0,
        M_PI
    );

    // -----------------------------------------
    // Left side approaching centromere
    // -----------------------------------------

    cairo_line_to(
        cr,
        x,
        constrictionBottom
    );

    // Gradually move inward
    cairo_line_to(
        cr,
        x + constrictionAmount,
        centromereY
    );

    // Gradually move outward
    cairo_line_to(
        cr,
        x,
        constrictionTop
    );

    // Close chromosome
    cairo_close_path(cr);

    // Fill
    cairo_set_source_rgb(
        cr,
        color.r,
        color.g,
        color.b
    );

    cairo_fill_preserve(cr);

    // Outline
    cairo_set_source_rgb(
        cr,
        0.1,
        0.1,
        0.1
    );

    cairo_set_line_width(cr, 1.5);

    if (dashedOutline)
    {
        const double dashPattern[2] = { 4.0, 3.0 };
        cairo_set_dash(cr, dashPattern, 2, 0.0);
    }

    cairo_stroke(cr);

    if (dashedOutline)
    {
        cairo_set_dash(cr, nullptr, 0, 0.0);                  // Back to solid for everything drawn afterwards
    }


    // -----------------------------------------
    // Drawing centromere ellipse
    // -----------------------------------------

    const double centromereStartY = y + height * centromereStart;

    const double centromereEndY = y + height * centromereEnd;

    const double centromereEllipseCenterY = (centromereStartY + centromereEndY) / 2.0;

    const double ellipseWidth = width + 2.0;

    const double ellipseHeight = std::max(5.0, centromereEndY - centromereStartY);

    cairo_save(cr);

    cairo_translate(cr, x + width / 2.0, centromereEllipseCenterY);

    cairo_scale(cr, ellipseWidth / 2.0, ellipseHeight / 2.0);

    cairo_arc(cr, 0.0, 0.0, 1.0, 0.0, 2.0 * M_PI);

    cairo_restore(cr);

    cairo_set_source_rgb(cr, 0.6, 0.6, 0.6); // Gray centromeres

    cairo_fill(cr);

}














void Karyogram::drawDoubleStrandBreakMarker(					// Function to draw black lines as DSB markers on th chromosomes
    cairo_t* cr,
    double x,									// DSB marker X location
    double y,									// DSB marker Y location
    double chromosomeWidth)							// DSB marker will span slightly more than the width of the chromosome for visual clarity.
{

    const double markerWidth = chromosomeWidth + 8.0; 				// Extend DSB markers past the chromosome width for clarity.
    const double markerLineWidth = 1.0;						// Thickness of DSB line marker

    // Black DSB line marker
    cairo_set_source_rgb(
        cr,
        0.0,
        0.0,
        0.0
    );

    cairo_set_line_width(
        cr,
        markerLineWidth
    );

    cairo_move_to(
        cr,
        x - markerWidth / 2.0,
        y
    );

    cairo_line_to(
        cr,
        x + markerWidth / 2.0,
        y
    );

    cairo_stroke(cr);

}










void Karyogram::drawSingleStrandBreakMarker(						// Function to draw white line SSB markers with marker length depending on the number of SSBs in the damage site.
    cairo_t* cr,
    double x,										// SSB marker x pixel position
    double y,										// SSB marker y pixel position
    double markerLength)								// The scaled SSB white line marker lengths
{
    // White SSB marker
    cairo_set_source_rgb(
        cr,
        1.0,
        1.0,
        1.0
    );

    cairo_set_line_width(cr, 1.5);

    cairo_move_to(
        cr,
        x - markerLength / 2.0,
        y
    );

    cairo_line_to(
        cr,
        x + markerLength / 2.0,
        y
    );

    cairo_stroke(cr);
}





} // namespace sddparser
