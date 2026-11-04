#include "Karyogram.h"

#include <cmath>
#include <algorithm>




namespace karyoke
{

// Function to generate colors for each chromosome by sampling the color hues
// equally around the entire color wheel for each chromosome listed, used to 
// generate SDD and SDR karyograms

RGB Karyogram::generateChromosomeColor(
    int chromosomeNumber, 							// Assign a different color to each chromosome number
    int totalChromosomes)							// Sample 'totalChromosomes' different colors for each chromosome
{

/*    // Golden-ratio color hue spacing gives better visual separation between consecutive chromosome colors.
    const double goldenRatio = 0.618033988749895;

    double hue = std::fmod(chromosomeNumber * goldenRatio, 1.0);
*/

    const double hue = static_cast<double>(chromosomeNumber) / static_cast<double>(totalChromosomes);		// Generate N different equally spaced colors for N chromosomes.

    const double saturation = 0.70;
    const double value = 0.90;

    double c = value * saturation;

    double h = hue * 6.0;

    double x = c * (1.0 - std::abs(std::fmod(h, 2.0) - 1.0));

    double r = 0.0;
    double g = 0.0;
    double b = 0.0;

    if (h < 1.0)
    {
        r = c;
        g = x;
    }
    else if (h < 2.0)
    {
        r = x;
        g = c;
    }
    else if (h < 3.0)
    {
        g = c;
        b = x;
    }
    else if (h < 4.0)
    {
        g = x;
        b = c;
    }
    else if (h < 5.0)
    {
        r = x;
        b = c;
    }
    else
    {
        r = c;
        b = x;
    }

    double m = value - c;

    r += m;
    g += m;
    b += m;

    return {r, g, b};
}













// When mapping rearranged chromosome fragments, need to keep track of the chromosome
// the rearranged fragment originally came from. The fragment can then be drawn on 
// the new chromosome with the original color to easily track the structural variations.
RGB Karyogram::getColorForOriginalStrand(
    int oldStrandID,						// Determine original strand ID color
    const SDRmasterHeader& masterHeader)			// Need SDR master header for the sizes of the chromosomes to generate the colors
{
    const std::vector<double>& sizes = masterHeader.intactChromosomeSizes;

    const int chromosomeCount = sizes.empty() ? 0 : static_cast<int>(sizes[0]);

    if (chromosomeCount < 2)
    {
        return generateChromosomeColor(oldStrandID, std::max(1, chromosomeCount));
    }

    const ChromosomeLayout layout = determineChromosomeLayout(sizes);					// Check if chromosome sizes were listed in split homologous, adjacent-homologous, or non-homologous layouts.
    const bool hasHomologs = layout != ChromosomeLayout::NON_HOMOLOGOUS;				// Default chromosome sizes layout is non-homologous.
    const int homologousPairs = hasHomologs ? (chromosomeCount - 2) / 2 : 0;				// Number of homologous pairs assuming a Y and X chromosome are specified.
    const int drawableGroups = hasHomologs ? homologousPairs : chromosomeCount - 2;			// Drawable groups concerns only the non-sex chromosomes
    const int totalColorGroups = drawableGroups + 2;							// Total color groups across all homologs and sex chromosomes
    const int yStrandID = chromosomeCount - 1; 								// Assumes Y is listed second last (1-indexed)
    const int xStrandID = chromosomeCount;								// Assumes X is listed last (1-indexed)

    int group = oldStrandID;                    							// Default is non-homologous layout, every strand gets their own color

    // Determine the color associated according to a given chromosome/homolog group.
    if(oldStrandID == yStrandID)
    {
        group = drawableGroups;
    }
    else if (oldStrandID == xStrandID)
    {
        group = drawableGroups + 1;
    }
    else if (layout == ChromosomeLayout::SPLIT_HOMOLOGS && homologousPairs > 0)
    {
        group = oldStrandID % homologousPairs;
    }
    else if (layout == ChromosomeLayout::ADJACENT_HOMOLOGS && homologousPairs > 0)
    {
        group = (oldStrandID - 1) / 2;				// Account for 1-indexing with the -1
    }

    if (totalColorGroups <= 0)
    {
        return generateChromosomeColor(oldStrandID, std::max(1, chromosomeCount));
    }

    return generateChromosomeColor(group, totalColorGroups);
}

// Homologs share one color; the SECOND copy of each chromosome is told apart by a dashed outline instead (see drawPaintedFragment). X and Y have no partner.
bool Karyogram::isSecondHomologStrand(int oldStrandID, const SDRmasterHeader& masterHeader)
{
    const std::vector<double>& sizes = masterHeader.intactChromosomeSizes;
    const int chromosomeCount = sizes.empty() ? 0 : static_cast<int>(sizes[0]);

    if (chromosomeCount < 2)
    {
        return false;
    }

    const ChromosomeLayout layout = determineChromosomeLayout(sizes);
    const int homologousPairs = (layout != ChromosomeLayout::NON_HOMOLOGOUS) ? (chromosomeCount - 2) / 2 : 0;

    if (homologousPairs <= 0 || oldStrandID == chromosomeCount - 1 || oldStrandID == chromosomeCount || oldStrandID < 1)	// No pairs, or Y / X
    {
        return false;
    }

    if (layout == ChromosomeLayout::SPLIT_HOMOLOGS)
    {
        return oldStrandID > homologousPairs;				// Strands 1..N first copies, N+1..2N second copies
    }

    return oldStrandID % 2 == 0;						// ADJACENT_HOMOLOGS: 1,3,5... first copy, 2,4,6... second copy
}

} // namespace karyoke
