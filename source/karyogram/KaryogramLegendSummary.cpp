#include <cairo/cairo.h>
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <map>
#include <set>

#include "Karyogram.h"
#include "SDDutilities.h"
#include "SDRutilities.h"



namespace karyoke
{


// ------------------------------------------------------------ //
//      DRAW SDD SUMMARY BOX AND LEGEND IN KARYOGRAM            //
// ------------------------------------------------------------ //


void Karyogram::drawSDDsummary(                                            // Box at top of karyogram to summarize the following details in the associated SDD file:
    cairo_t* cr,
    const std::vector<double>& cellCyclePhase,                             // Cell cycle phase from SDD file header to be summarized
    const std::vector<Exposure>& exposures,                                // Vector of exposures to summarize the total number of SSBs and DSBs
    const std::vector<double>& doseOrFluence,                              // Dose or fluence from SDD file to be summarized
    const std::vector<int>& incidentParticles)                             // Incident particles from SDD file header to be summarized
{

    // --------------------------------------------------
    // Summary box position and dimensions
    // --------------------------------------------------
    const double summaryX = 50.0;
    const double summaryY = 25.0;

    const double summaryWidth = 900.0;
    const double summaryHeight = 128.0;


    // --------------------------------------------------
    // Count DSBs and SSBs
    // --------------------------------------------------
    int numberOfDSBs = 0;
    int numberOfSSBs = 0;

    for (const auto& exposure : exposures)
    {
        for (const auto& damage : exposure.damages)
        {
            if (damage.damageType.presenceOfDoubleStrandBreaks == 1)
            {
                ++numberOfDSBs;
            }

            numberOfSSBs += damage.damageType.numSingleBackboneBreaks;
        }
    }


    // -------------------------------------
    // Determine cell cycle phase
    // -------------------------------------
    std::string phaseLabel = cellCyclePhaseMeaning(cellCyclePhase);        // Helper function in SDDutilities to interpret cell cycle phase code in SDD header.


    // -------------------------------------
    // Determine Dose or Fluence
    // -------------------------------------
    std::string doseOrFluenceSummary = doseOrFluenceMeaning(doseOrFluence);// Helper function in SDDutilities to interpret dose or fluence in SDD header


    // -------------------------------------
    // Determine incident particle(s)
    // -------------------------------------
    std::string incidentParticlesLabel;
    if (incidentParticles.empty())
    {
        incidentParticlesLabel = "N/A";
    }
    else
    {
        for (std::size_t i = 0; i < incidentParticles.size(); i++)
        {
            incidentParticlesLabel += incidentParticlesMeaning(incidentParticles[i]);
            if (i + 1 < incidentParticles.size())
            {
                incidentParticlesLabel += ", ";
            }
        }
    }


    // --------------------------------------------------
    // Draw summary box
    // --------------------------------------------------

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

    cairo_set_line_width(cr, 1.5);

    cairo_rectangle(cr, summaryX, summaryY, summaryWidth, summaryHeight);

    cairo_stroke(cr);


    // --------------------------------------------------
    // Draw summarized SDD header field summaries in box
    // --------------------------------------------------

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);

    cairo_set_font_size(cr, 20.0);

    // Three rows of text
    const double firstRowY = summaryY + 30.0;
    const double secondRowY = summaryY + 68.0;
    const double thirdRowY = summaryY + 106.0;


    // ------------------------------------------
    // TOP ROW
    // ------------------------------------------

    // Cell ID
    std::string cellIDlabel = exposures.empty() ? "N/A" : std::to_string(exposures[0].exposureID - 1);              // 0-index cell ID by subtracting 1 from the number of exposures.

    cairo_move_to(cr, summaryX + 15.0, firstRowY);

    cairo_show_text(cr, ("Cell ID: " + cellIDlabel).c_str());


    // Incident Particle(s)
    cairo_move_to(cr, summaryX + 470.0, firstRowY);

    cairo_show_text(cr, ("Incident Particle(s): " + incidentParticlesLabel).c_str());


    // -----------------------------------------
    // SECOND ROW
    // -----------------------------------------

    // Dose / fluence
    cairo_move_to(cr, summaryX + 15.0, secondRowY);

    cairo_show_text(cr, ("Dose/Fluence: " + doseOrFluenceSummary).c_str());

    // Cell cycle phase
    cairo_move_to(cr, summaryX + 470.0, secondRowY);

    cairo_show_text(cr, ("Cell Cycle Phase: " + phaseLabel).c_str());


    // --------------------------------------------------
    // THIRD ROW
    // --------------------------------------------------

    // SSBs
    cairo_move_to(cr, summaryX + 15.0, thirdRowY);

    cairo_show_text(cr, ("SSBs: " + std::to_string(numberOfSSBs)).c_str());

    // DSBs
    cairo_move_to(cr, summaryX + 470.0, thirdRowY);

    cairo_show_text(cr, ("DSBs: " + std::to_string(numberOfDSBs)).c_str());
}








void Karyogram::drawSDDlegend(                                     // Function to draw a legend of symbols at the bottom of the karyogram
    cairo_t* cr,
    double legendY)                                             // Height of the legend in pixels, adjusted depending on the height of the image, not hardcoded
{
    // --------------------------------------------------
    // Legend position
    // --------------------------------------------------
    const double legendX = 50.0;                                // Hardcoded for now, may decide to adjust based on number of symbols like I did for legendY.


    // --------------------------------------------------
    // Legend border
    // --------------------------------------------------
    const double legendWidth = 900.0;
    const double legendHeight = 70.0;

    const double borderPaddingX = 15.0;
    const double borderPaddingY = 10.0;

    // Center chromosome vertically inside legend
    const double legendTop = legendY - 25.0 - borderPaddingY;
    const double legendCenterY = legendTop + legendHeight / 2.0;

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_set_line_width(cr, 1.5);
    cairo_rectangle(cr, legendX - borderPaddingX, legendY - 25.0 - borderPaddingY, legendWidth, legendHeight);
    cairo_stroke(cr);

    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 16.0);

    // Centering legend labels
    cairo_font_extents_t fontExtents;
    cairo_font_extents(cr, &fontExtents);

    const double textBaseline = legendCenterY +
    (fontExtents.ascent - fontExtents.descent) / 2.0;


    // --------------------------------------------------
    // Legend title
    // --------------------------------------------------
    cairo_move_to(cr, legendX, textBaseline);
    cairo_show_text(cr, "Symbols:");


    // -------------------------------------------------
    // Draw generic chromosome symbol
    // -------------------------------------------------
    const double chromosomeWidth = 14.0;
    const double chromosomeHeight = 40.0;

    // Position chromosome by its center
    const double chromosomeCenterX = 150.0;
    const double chromosomeX = chromosomeCenterX - chromosomeWidth / 2.0;
    const double chromosomeY = legendCenterY - chromosomeHeight / 2.0;

    // Generic centromere location
    const double centromereCenter = 0.5;
    const double constrictionAmount = 1.5;
    const double capRadius = chromosomeWidth / 2.0;
    const double centerY = chromosomeY + chromosomeHeight * centromereCenter;
    const double constrictionHeight = 5.0;
    const double constrictionTop = centerY - constrictionHeight / 2.0;
    const double constrictionBottom = centerY + constrictionHeight / 2.0;

    cairo_new_path(cr);

    // Top cap
    cairo_arc(cr, chromosomeX + capRadius, chromosomeY + capRadius, capRadius, M_PI, 2.0 * M_PI);

    // Right side
    cairo_line_to(cr, chromosomeX + chromosomeWidth, constrictionTop);
    cairo_line_to(cr, chromosomeX + chromosomeWidth - constrictionAmount, centerY);
    cairo_line_to(cr, chromosomeX + chromosomeWidth, constrictionBottom);
    cairo_line_to(cr, chromosomeX + chromosomeWidth, chromosomeY + chromosomeHeight - capRadius);

    // Bottom cap
    cairo_arc(cr, chromosomeX + capRadius, chromosomeY + chromosomeHeight - capRadius, capRadius, 0.0, M_PI);

    // Left side
    cairo_line_to(cr, chromosomeX, constrictionBottom);
    cairo_line_to(cr, chromosomeX + constrictionAmount, centerY);
    cairo_line_to(cr, chromosomeX, constrictionTop);
    cairo_close_path(cr);

    // Generic chromosome colour
    RGB chromosomeColor = generateChromosomeColor(0, 24);       // total number of chromosomes does not matter for generic chromosome drawing
    cairo_set_source_rgb(cr, chromosomeColor.r, chromosomeColor.g, chromosomeColor.b);
    cairo_fill_preserve(cr);

    // Outline
    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    cairo_set_line_width(cr, 1.0);
    cairo_stroke(cr);

    // Chromosome label
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_move_to(cr, chromosomeX + chromosomeWidth + 10.0, textBaseline);
    cairo_show_text(cr, "Chromosome");


    // --------------------------------------------------
    // DSB marker symbol
    // --------------------------------------------------
    const double dsbX = 330.0;
    const double dsbY = legendCenterY;

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_set_line_width(cr, 2.0);
    cairo_move_to(cr, dsbX - 10.0, dsbY);
    cairo_line_to(cr, dsbX + 10.0, dsbY);
    cairo_stroke(cr);

    // DSB label
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_move_to(cr, dsbX + 14.0, textBaseline);
    cairo_show_text(cr, "Double-Strand Break");

    // --------------------------------------------------
    // Single-Strand Break symbol
    // --------------------------------------------------
    const double ssbX = 550.0;
    const double ssbY = legendCenterY;
    const double ssbMarkerWidth = 14.0;
    const double ssbMarkerHeight = 3.0;

    // Draw white rectangle
    cairo_rectangle(cr, ssbX - ssbMarkerWidth / 2.0, ssbY - ssbMarkerHeight / 2.0, ssbMarkerWidth, ssbMarkerHeight);
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_fill_preserve(cr);

    // Black outline
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_set_line_width(cr, 1.0);
    cairo_stroke(cr);

    // SSB label
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_move_to(cr, ssbX + 10.0, textBaseline);
    cairo_show_text(cr, "Single-Strand Break");


    // --------------------------------------------------
    // Centromere marker symbol
    // --------------------------------------------------
    const double centromereX = 770.0;
    const double centromereY = legendCenterY;
    const double ellipseWidth = 14.0;
    const double ellipseHeight = 8.0;

    cairo_save(cr);

    cairo_translate(cr, centromereX, centromereY);
    cairo_scale(cr, ellipseWidth / 2.0, ellipseHeight / 2.0);
    cairo_arc(cr, 0.0, 0.0, 1.0, 0.0, 2.0 * M_PI);
    cairo_restore(cr);

    // Gray fill
    cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
    cairo_fill_preserve(cr);

    // Centromere label
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_move_to(cr, centromereX + 10.0, textBaseline);
    cairo_show_text(cr, "Centromere");
}












// ------------------------------------------------------------ //
//      DRAW SDR SUMMARY BOX AND LEGEND IN KARYOGRAM            //
// ------------------------------------------------------------ //

void Karyogram::drawSDRsummary(cairo_t* cr, const SDRmasterHeader& masterHeader, const SDRsubHeader& rawSubHeader)
{
    const int numOriginalStrands = masterHeader.intactChromosomeSizes.empty() ? 0 : static_cast<int>(masterHeader.intactChromosomeSizes[0]);

    // Indels (<= 50 bp) are counted on the cell as written; every other mutation type is detected with the indels taken out (see removeIndels).
    std::vector<SDRindelEvent> indels;
    const SDRsubHeader subHeader = removeIndels(rawSubHeader, numOriginalStrands, &indels);

    const std::vector<SDRdeletionEvent> deletions = detectDeletions(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRinversionEvent> inversions = detectInversions(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRtranslocationEvent> translocations = detectTranslocations(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRecDNAevent> ecDNA = detectECDNA(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRdeletionInversionEvent> deletionInversion = detectDeletionInversions(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRdeletionTranslocationEvent> deletionTranslocation = detectDeletionTranslocations(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRdeletionInsertionEvent> deletionInsertion = detectDeletionInsertions(subHeader, numOriginalStrands);
    const std::vector<SDRchromoplexyEvent> chromoplexy = detectChromoplexy(subHeader, numOriginalStrands, masterHeader);
    const std::vector<SDRchromothripsisEvent> chromothripsis = detectChromothripsis(subHeader, numOriginalStrands);
    const std::size_t complexRearrangements = countComplexGenomicRearrangements(subHeader, masterHeader);	// Includes changes no other type accounts for
    const std::vector<SDRchromosomeRingEvent> chromosomeRings = detectChromosomeRings(subHeader, numOriginalStrands);
    const std::vector<SDRdicentricEvent> dicentrics = detectDicentrics(subHeader, numOriginalStrands);

    const double summaryX = 50.0;
    const double summaryY = 25.0;
    const double summaryWidth = 900.0;
    const double summaryHeight = 193.0;						// 5 rows of counts; renderSDRkaryogram's startY sits 70 px below the box.

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_set_line_width(cr, 1.5);
    cairo_rectangle(cr, summaryX, summaryY, summaryWidth, summaryHeight);
    cairo_stroke(cr);

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 20.0);

    const double firstRowY = summaryY + 30.0;
    const double secondRowY = summaryY + 68.0;
    const double thirdRowY = summaryY + 106.0;
    const double fourthRowY = summaryY + 144.0;
    const double fifthRowY = summaryY + 182.0;

    cairo_move_to(cr, summaryX + 15.0, firstRowY);
    cairo_show_text(cr, ("Cell ID: " + std::to_string(subHeader.cellID)).c_str());

    cairo_move_to(cr, summaryX + 300.0, firstRowY);
    cairo_show_text(cr, ("Chromosome Rings: " + std::to_string(chromosomeRings.size())).c_str());

    cairo_move_to(cr, summaryX + 615.0, firstRowY);
    cairo_show_text(cr, ("Dicentrics: " + std::to_string(dicentrics.size())).c_str());

    cairo_move_to(cr, summaryX + 15.0, secondRowY);
    cairo_show_text(cr, ("Long Deletions: " + std::to_string(deletions.size())).c_str());

    cairo_move_to(cr, summaryX + 300.0, secondRowY);
    cairo_show_text(cr, ("Balanced Inversions: " + std::to_string(inversions.size())).c_str());

    cairo_move_to(cr, summaryX + 615.0, secondRowY);
    cairo_show_text(cr, ("Balanced Translocations: " + std::to_string(translocations.size())).c_str());

    cairo_move_to(cr, summaryX + 15.0, thirdRowY);
    cairo_show_text(cr, ("ecDNA: " + std::to_string(ecDNA.size())).c_str());

    cairo_move_to(cr, summaryX + 300.0, thirdRowY);
    cairo_show_text(cr, ("Deletion-Inversions: " + std::to_string(deletionInversion.size())).c_str());

    cairo_move_to(cr, summaryX + 615.0, thirdRowY);
    cairo_show_text(cr, ("Deletion-Translocations: " + std::to_string(deletionTranslocation.size())).c_str());

    cairo_move_to(cr, summaryX + 15.0, fourthRowY);
    cairo_show_text(cr, ("Deletion-Insertions: " + std::to_string(deletionInsertion.size())).c_str());

    cairo_move_to(cr, summaryX + 300.0, fourthRowY);
    cairo_show_text(cr, ("Chromoplexy: " + std::to_string(chromoplexy.size())).c_str());

    cairo_move_to(cr, summaryX + 615.0, fourthRowY);
    cairo_show_text(cr, ("Chromothripsis: " + std::to_string(chromothripsis.size())).c_str());

    cairo_move_to(cr, summaryX + 15.0, fifthRowY);
    cairo_show_text(cr, ("Complex Genomic Rearrangements: " + std::to_string(complexRearrangements)).c_str());

    cairo_move_to(cr, summaryX + 615.0, fifthRowY);						// Third column: the CGR label is too wide to share its row's middle column.
    cairo_show_text(cr, ("Indels (<= 50 bp): " + std::to_string(indels.size())).c_str());

}










void Karyogram::drawSDRlegend(cairo_t* cr, double legendY)
{
    const double legendX = 50.0;
    const double legendWidth = 930.0;
    const double legendHeight = 70.0;
    const double borderPaddingX = 15.0;
    const double borderPaddingY = 10.0;

    const double legendTop = legendY - 25.0 - borderPaddingY;
    const double legendCenterY = legendTop + legendHeight / 2.0;

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_set_line_width(cr, 1.5);
    cairo_rectangle(cr, legendX - borderPaddingX, legendY - 25.0 - borderPaddingY, legendWidth, legendHeight);
    cairo_stroke(cr);

    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 20.0);

    cairo_font_extents_t fontExtents;
    cairo_font_extents(cr, &fontExtents);
    const double textBaseline = legendCenterY + (fontExtents.ascent - fontExtents.descent) / 2.0;

    cairo_move_to(cr, legendX, textBaseline);
    cairo_show_text(cr, "Symbols:");

    // --------------------------------------------------
    // Generic chromosome symbol - identical shape/logic to drawLegend()'s
    // --------------------------------------------------
    const double chromosomeWidth = 14.0;
    const double chromosomeHeight = 40.0;
    const double chromosomeCenterX = 180.0;
    const double chromosomeX = chromosomeCenterX - chromosomeWidth / 2.0;
    const double chromosomeY = legendCenterY - chromosomeHeight / 2.0;

    const double centromereCenter = 0.5;
    const double constrictionAmount = 1.5;
    const double capRadius = chromosomeWidth / 2.0;
    const double centerY = chromosomeY + chromosomeHeight * centromereCenter;
    const double constrictionHeight = 5.0;
    const double constrictionTop = centerY - constrictionHeight / 2.0;
    const double constrictionBottom = centerY + constrictionHeight / 2.0;

    cairo_new_path(cr);
    cairo_arc(cr, chromosomeX + capRadius, chromosomeY + capRadius, capRadius, M_PI, 2.0 * M_PI);
    cairo_line_to(cr, chromosomeX + chromosomeWidth, constrictionTop);
    cairo_line_to(cr, chromosomeX + chromosomeWidth - constrictionAmount, centerY);
    cairo_line_to(cr, chromosomeX + chromosomeWidth, constrictionBottom);
    cairo_line_to(cr, chromosomeX + chromosomeWidth, chromosomeY + chromosomeHeight - capRadius);
    cairo_arc(cr, chromosomeX + capRadius, chromosomeY + chromosomeHeight - capRadius, capRadius, 0.0, M_PI);
    cairo_line_to(cr, chromosomeX, constrictionBottom);
    cairo_line_to(cr, chromosomeX + constrictionAmount, centerY);
    cairo_line_to(cr, chromosomeX, constrictionTop);
    cairo_close_path(cr);

    RGB chromosomeColor = generateChromosomeColor(0, 24);

    cairo_set_source_rgb(cr, chromosomeColor.r, chromosomeColor.g, chromosomeColor.b);
    cairo_fill_preserve(cr);

    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    cairo_set_line_width(cr, 1.0);
    cairo_stroke(cr);

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_move_to(cr, chromosomeX + chromosomeWidth + 10.0, textBaseline);
    cairo_show_text(cr, "Chromosome");

    // --------------------------------------------------
    // Centromere symbol
    // --------------------------------------------------
    const double centromereX = 400.0;
    const double ellipseWidth = 14.0;
    const double ellipseHeight = 8.0;

    cairo_new_path(cr);

    cairo_save(cr);
    cairo_translate(cr, centromereX, legendCenterY);
    cairo_scale(cr, ellipseWidth / 2.0, ellipseHeight / 2.0);
    cairo_arc(cr, 0.0, 0.0, 1.0, 0.0, 2.0 * M_PI);
    cairo_restore(cr);


    cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
    cairo_fill_preserve(cr);

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_move_to(cr, centromereX + 10.0, textBaseline);
    cairo_show_text(cr, "Centromere");

    // --------------------------------------------------
    // Balanced inversion symbol - a small colored square with the
    // same chevron used on the karyogram itself
    // --------------------------------------------------
    const double inversionX = 600.0;
    const double inversionSize = 16.0;

    cairo_new_path(cr);

    cairo_set_source_rgb(cr, chromosomeColor.r, chromosomeColor.g, chromosomeColor.b);
    cairo_rectangle(cr, inversionX - inversionSize / 2.0, legendCenterY - inversionSize / 2.0, inversionSize, inversionSize);
    cairo_fill_preserve(cr);

    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    cairo_set_line_width(cr, 1.0);
    cairo_stroke(cr);

    drawInversionChevron(cr, inversionX, legendCenterY - inversionSize / 2.0, legendCenterY + inversionSize / 2.0, inversionSize);	// Black chevron spanning the whole square, like on the karyogram

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_move_to(cr, inversionX + 15.0, textBaseline);
    cairo_show_text(cr, "Reversed Segment");


    // --------------------------------------------------
    // ecDNA symbol - same ring shape used on the karyogram itself
    // --------------------------------------------------
    const double ecDNAX = 850.0;
    const double ecDNADiameter = 28.0;

    cairo_new_path(cr);

    drawCircularFragment(cr, ecDNAX, legendCenterY, ecDNADiameter, chromosomeColor);

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    cairo_move_to(cr, ecDNAX + 18.0, textBaseline);
    cairo_show_text(cr, "ecDNA");

}




} // namespace karyoke
