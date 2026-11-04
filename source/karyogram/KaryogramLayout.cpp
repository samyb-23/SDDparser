#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>

#include "Karyogram.h"



namespace karyoke 
{

ChromosomeLayout Karyogram::determineChromosomeLayout(
    const std::vector<double>& chromosomeSizes)							// Access chromosome sizes to determine how they were passed by the user, adjacent-homologs, split-homologs, non-homologous.
{
    if (chromosomeSizes.empty())
    {
        std::cerr << "Require Chromosome sizes in mega base pairs to plot karyogram.\n";
        return ChromosomeLayout::NON_HOMOLOGOUS;
    }

    const int chromosomeCount = static_cast<int>(chromosomeSizes[0]);						// The first entry is the number of chromosomes.

    if (chromosomeCount < 2 || chromosomeSizes.size() != static_cast<size_t>(chromosomeCount + 1))		// Need at least two chromosomes to determine a mapping.
    {
        return ChromosomeLayout::NON_HOMOLOGOUS;
    }

    // --------------------------------------------------
    // Helper for comparing chromosome sizes.
    // Add a tolerance to the homologous chromosome sizes
    // of 1e-6 * chromosome size <= 250 bases.
    // --------------------------------------------------
    const auto approximatelyEqual = [](double a, double b)
        {
            const double tolerance = 1e-6;
            return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
        };

    // --------------------------------------------------
    // Case 1: Adjacent homologs
    //
    // chromosomeSizes:
    // [count, 1, 1, 2, 2, 3, 3, ...]
    // --------------------------------------------------

    if (chromosomeCount >= 2 && approximatelyEqual(chromosomeSizes[1], chromosomeSizes[2]))			// If neighbouring chromosome sizes are approximately equal --> adjacent homologous layout
    {
        return ChromosomeLayout::ADJACENT_HOMOLOGS;
    }


    // --------------------------------------------------
    // Case 2: Split homologs
    //
    // chromosomeSizes:
    // [count,1,2,3,...,N,1,2,3,...,N]
    //
    // Compare the first chromosome with the chromosome
    // halfway through the chromosome list.
    // --------------------------------------------------

    if (chromosomeCount >= 4 && (chromosomeCount - 2) % 2 == 0)
    {
        const int homologousSetSize = (chromosomeCount - 2) / 2;

        const size_t firstHomologIndex = 1;

        const size_t secondHomologIndex = static_cast<size_t>(homologousSetSize + 1);

        if (approximatelyEqual(chromosomeSizes[firstHomologIndex], chromosomeSizes[secondHomologIndex]))	// If the first chromosome size is approx equal to the chromosome size halfway through the list (minus the two sex chromosomes) --> split homologous set
        {
            return ChromosomeLayout::SPLIT_HOMOLOGS;
        }
    }


    // --------------------------------------------------
    // Case 3: Non-homologous
    //
    // chromosomeSizes:
    // [count,1,2,3,4,...,N]
    // --------------------------------------------------

    return ChromosomeLayout::NON_HOMOLOGOUS;									// Default case.

}
















const std::vector<CentromerePosition> humanCentromeres =                   // Vector storing the associated chromosome ID centromere positions to
{                                                                          // draw the centromere location if the human option of '--karyotype' is passed by the user in the command-line
    {"1", 121535434, 124535434},                                           // Each vector element contains {'Chromosome ID', Centromere Start Position (bp), Centromere End Position (bp)}
    {"2",  92326171,  95326171},
    {"3",  90504854,  93504854},
    {"4",  49660117,  52660117},
    {"5",  46405641,  49405641},
    {"6",  58830166,  61830166},
    {"7",  58054331,  61054331},
    {"8",  43838887,  46838887},
    {"9",  47367679,  50367679},
    {"10", 39254935,  42245935},
    {"11", 51644205,  54644205},
    {"12", 34856694,  37856694},
    {"13", 16000000,  19000000},
    {"14", 16000000,  19000000},
    {"15", 17000000,  20000000},
    {"16", 35335801,  38335801},
    {"17", 22263006,  25263006},
    {"18", 15460898,  18460898},
    {"19", 24681782,  27681782},
    {"20", 26369569,  29369569},
    {"21", 11288129,  14288129},
    {"22", 13000000,  16000000},
    {"Y", 10100001, 10327000},
    {"X", 58632012, 61632012}
};










const CentromerePosition* Karyogram::getHumanCentromere(				// Function to return the centromere start and end positions for each chromosome number. 
    int chromosomeID)									// Obtain centromere ranges for a given chromosome ID

{
    std::string chromosome;

    if (chromosomeID >= 1 && chromosomeID <= 22)
    {
        chromosome = std::to_string(chromosomeID);
    }
    else if (chromosomeID >= 23 && chromosomeID <= 44)						// Chromosome number 23-44 correspond to the homologs of chromosomes 1-22. Y is chromosome 45 and X is chromosome 46 --> MUST MAINTAIN THIS ORDERING
    {
        // Second homologous set.
        chromosome = std::to_string(chromosomeID - 22);
    }
    else if (chromosomeID == 45)                                           			// Chromosome Number 45 = Y chromosome
    {
        chromosome = "Y";
    }
    else if (chromosomeID == 46)                                           			// Chromosome Number 46 = X chromosome
    {
        chromosome = "X";
    }
    else
    {
        return nullptr;                                                    			// If other Chromosome Numbers encountered, return nullptr.
    }

    for (const auto& centromere : humanCentromeres)                        			// Return the centromere start and end locations for each chromosome ID.
    {
        if (centromere.chromosome == chromosome)
        {
            return &centromere;
        }
    }

    return nullptr;
}











// For SDR file, looks up centromere location (in bp) for the original strand ID, return false if no centromere information
bool Karyogram::getCentromereForOriginalStrand(
    int oldStrandID,										// Need old strand ID to associate a given centromere range
    bool humanGenome,										// Did user pass a human genome? If not, use generic centromere regions at 50% of the chromosome length
    const SDRmasterHeader& masterHeader,							// Access SDR master header for chromosome sizes to decide where to draw the centromeres
    double& centromereStartBP,									// Access centromere range start in BP, converted later to pixels
    double& centromereEndBP)									// Access centromere range end in BP, converted later to pixels
{

    if (humanGenome)
    {
	// Modify strand ID based on how chromosome sizes are passed
	int relativeStrandID = oldStrandID;

        const std::vector<double>& sizes = masterHeader.intactChromosomeSizes;
        const int chromosomeCount = sizes.empty() ? 0 : static_cast<int>(sizes[0]);

	if (chromosomeCount >= 2)
        {
	    const ChromosomeLayout layout = determineChromosomeLayout(sizes);
            const bool hasHomologs = layout != ChromosomeLayout::NON_HOMOLOGOUS;
            const int homologousPairs = hasHomologs ? (chromosomeCount - 2) / 2 : 0;
            const int yStrandID = chromosomeCount - 1;
            const int xStrandID = chromosomeCount;

            if (layout == ChromosomeLayout::ADJACENT_HOMOLOGS)
            {
                if (oldStrandID == yStrandID)
                {
                    relativeStrandID = 2 * homologousPairs + 1;
                }
                else if (oldStrandID == xStrandID)
                {
                    relativeStrandID = 2 * homologousPairs + 2;
                }
		else if (homologousPairs > 0)
                {
                    // Adjacent pairs: strand 2i+1 is chromosome i+1's first
                    // copy, strand 2i+2 is its second copy - both map to
                    // the same i via integer division, since (2i+1-1)/2=i
                    // and (2i+2-1)/2=i as well.
                    const int i = (oldStrandID - 1) / 2;
                    const bool isSecondCopy = (oldStrandID % 2 == 0);

                    relativeStrandID = isSecondCopy ? (i + 1 + homologousPairs) : (i + 1);
                }
            }
	    else if (layout == ChromosomeLayout::NON_HOMOLOGOUS)
            {
                // No homolog duplication at all - autosomes (1..drawableGroups)
                // already match getHumanCentromere's convention directly, but
                // Y and X sit wherever the file's total count puts them
                // (e.g. 23/24 for a 24-strand file), NOT at the table's
                // fixed 45/46 - those two need explicit remapping.
                if (oldStrandID == yStrandID)
                {
                    relativeStrandID = 45;
                }
                else if (oldStrandID == xStrandID)
                {
                    relativeStrandID = 46;
                }
            }
	    else if (layout == ChromosomeLayout::SPLIT_HOMOLOGS)
	    {
    		if (oldStrandID == yStrandID)
    		{
        	    relativeStrandID = 45;
    		}
    		else if (oldStrandID == xStrandID)
    		{
        	    relativeStrandID = 46;
    		}
    		else if (homologousPairs > 0 && oldStrandID > homologousPairs)
    		{
        	    // Second-copy homolog: strand (homologousPairs + k) is chromosome k's
        	    // second copy. getHumanCentromere expects second copies at k+22
        	    // specifically, regardless of how many pairs THIS file declares
        	    // (only coincides with oldStrandID directly when homologousPairs==22,
        	    // i.e. a full 46-chromosome file).
        	    const int k = oldStrandID - homologousPairs;
        	    relativeStrandID = k + 22;
    		}
    		// else: first-copy homolog (1 <= oldStrandID <= homologousPairs) already
    		// matches getHumanCentromere's 1-22 convention directly - no conversion needed.
	    }

        }

	const CentromerePosition* centromere = getHumanCentromere(relativeStrandID);

        if (centromere == nullptr)
        {
            return false;
        }

        centromereStartBP = static_cast<double>(centromere->start);
        centromereEndBP = static_cast<double>(centromere->end);

        return true;

    }
    else
    {
        const std::size_t sizeIndex = static_cast<std::size_t>(oldStrandID); // 0th index is number of chromosomes listed

        if (sizeIndex >= masterHeader.intactChromosomeSizes.size())
        {
            return false;
        }

        const double chromosomeSizeBP = masterHeader.intactChromosomeSizes[sizeIndex] * 1000000.0;  // Convert chromosome size in Mbp to bp

	// Generic centromere start and end locations halfway through the chromosome length for non-human chromosomes.
        centromereStartBP = 0.49 * chromosomeSizeBP;
        centromereEndBP = 0.51 * chromosomeSizeBP;

        return true;
    }
}


void Karyogram::drawChromosomeLabel(						// Draw a chromosome label and remember where it is, for the web viewer's clickable labels.
    cairo_t* cr,
    const std::string& label,
    double x,										// Text start position (baseline), as passed to cairo_move_to
    double y,
    const std::vector<int>& oldStrandIDs)
{
    cairo_move_to(cr, x, y);
    cairo_show_text(cr, label.c_str());

    // The box is padded a little so the label is easy to click even when zoomed out.
    cairo_text_extents_t extents;
    cairo_text_extents(cr, label.c_str(), &extents);
    const double padding = 6.0;

    KaryogramLabelRegion region;
    region.label = label;
    region.x = x + extents.x_bearing - padding;
    region.y = y + extents.y_bearing - padding;
    region.width = extents.width + 2.0 * padding;
    region.height = extents.height + 2.0 * padding;
    region.oldStrandIDs = oldStrandIDs;
    labelRegions.push_back(region);
}


// Cairo stream callback: appends each chunk of encoded image data (PNG or SVG) to the std::vector passed as closure.
static cairo_status_t appendImageChunk(void* closure, const unsigned char* data, unsigned int length)
{
    std::vector<unsigned char>* buffer = static_cast<std::vector<unsigned char>*>(closure);
    buffer->insert(buffer->end(), data, data + length);
    return CAIRO_STATUS_SUCCESS;
}


cairo_surface_t* Karyogram::createKaryogramSurface(					// Create the surface every karyogram is drawn on.
    int width,
    int height,
    const std::string& outputFilename,
    std::vector<unsigned char>* imageBytes,
    KaryogramFormat format)
{
    if (format == KaryogramFormat::SVG)
    {
        // SVG keeps every chromosome, marker and letter as a vector shape instead of pixels, so the image can be zoomed
        // without becoming blurry. Cairo writes the SVG text as it is drawn: into imageBytes, or straight to the file.
        cairo_surface_t* surface;

        if (imageBytes)
        {
            imageBytes->clear();
            surface = cairo_svg_surface_create_for_stream(appendImageChunk, imageBytes, width, height);
        }
        else
        {
            surface = cairo_svg_surface_create(outputFilename.c_str(), width, height);
        }

#if CAIRO_VERSION >= CAIRO_VERSION_ENCODE(1, 16, 0)
        // Newer cairo versions default to points (pt) for the SVG size, which browsers would show 1.33x larger than the PNG.
        // Using pixels keeps the SVG exactly the same size as the PNG.
        cairo_svg_surface_set_document_unit(surface, CAIRO_SVG_UNIT_PX);
#endif

        return surface;
    }

    return cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
}


bool Karyogram::writeKaryogramImage(							// Add the legend below a karyogram drawn on a recording surface, then write the image.
    cairo_surface_t* content,								// Recording surface holding everything but the legend; destroyed here.
    int width,
    void (Karyogram::*drawLegend)(cairo_t*, double),					// &Karyogram::drawSDRlegend or &Karyogram::drawSDDlegend
    const std::string& outputFilename,
    std::vector<unsigned char>* imageBytes,
    KaryogramFormat format)
{
    const double legendTopGap = 30.0;							// Space between the lowest label and the legend border, safe to tune.
    const double legendBottomMargin = 40.0;						// Space between the legend border and the bottom of the image, safe to tune.
    const double legendTopOffset = 35.0;						// drawSDRlegend/drawSDDlegend draw their border from legendY - 35 (25 + borderPaddingY) ...
    const double legendBottomOffset = 35.0;						// ... to legendY + 35 (legendHeight 70), keep in sync if those change.

    // Ink extents = the box around everything actually drawn, so the lowest label is found whatever the row heights, derivative lengths or layout
    // (no need to predict where each label lands). Safeguard: an empty recording reports zero size, the legend then simply starts at the top.
    double inkX, inkY, inkWidth, inkHeight;
    cairo_recording_surface_ink_extents(content, &inkX, &inkY, &inkWidth, &inkHeight);

    const double contentBottom = (inkHeight > 0.0) ? inkY + inkHeight : 0.0;
    const double legendY = std::ceil(contentBottom + legendTopGap + legendTopOffset);
    const int height = static_cast<int>(legendY + legendBottomOffset + legendBottomMargin);

    cairo_surface_t* surface = createKaryogramSurface(width, height, outputFilename, imageBytes, format);
    cairo_t* cr = cairo_create(surface);

    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);						// White background (the recording itself is transparent)
    cairo_paint(cr);

    cairo_set_source_surface(cr, content, 0.0, 0.0);					// Replay the karyogram at its original coordinates, so label hotspots still line up.
    cairo_paint(cr);									// Replayed as vector drawing commands, so SVG output stays sharp.

    cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);
    (this->*drawLegend)(cr, legendY);

    const bool written = writeKaryogramImage(surface, outputFilename, imageBytes, format);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    cairo_surface_destroy(content);

    return written;
}


bool Karyogram::writeKaryogramImage(							// Finish the karyogram image: to a file, or into memory if imageBytes is given.
    cairo_surface_t* surface,
    const std::string& outputFilename,
    std::vector<unsigned char>* imageBytes,
    KaryogramFormat format)
{
    cairo_status_t status;

    if (format == KaryogramFormat::SVG)
    {
        // The SVG was already written (to the file or imageBytes) while drawing; finishing the surface writes the closing tags.
        cairo_surface_finish(surface);
        status = cairo_surface_status(surface);
    }
    else if (imageBytes)
    {
        imageBytes->clear();
        status = cairo_surface_write_to_png_stream(surface, appendImageChunk, imageBytes);
    }
    else
    {
        status = cairo_surface_write_to_png(surface, outputFilename.c_str());
    }

    if (status != CAIRO_STATUS_SUCCESS)
    {
        std::cerr << "ERROR: Failed to write karyogram image: " << cairo_status_to_string(status) << "\n";
        return false;
    }

    return true;
}


} // namespace karyoke
