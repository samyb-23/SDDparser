#include "Karyogram.h"




namespace karyoke
{

// Function to return the location in base pairs for a double strand break associated with a given chromosome Number and exposure.
std::vector<DamageLocation> Karyogram::getDoubleStrandBreaks(
    const std::vector<Exposure>& exposures)				// Access exposures to determine the damage locations according to each chromosome ID
{
    std::vector<DamageLocation> dsbVec;                                 // Initialize the double-strand break vector to store base pair locations of DSBs per chromosome per exposure.

    for (const auto& exposure : exposures)                              // Loop through all exposures
    {
        for (const auto& damage : exposure.damages)                     // Loop through the damages in each exposure
        {
            if (damage.damageType.presenceOfDoubleStrandBreaks == 1)    // If double strand break is present, store chromosome ID, base pair location to dsbVec
            {
                DamageLocation dsb;                                     // Instantiate the DamageLocation object dsb for easy access to DSB location and chromosome ID.
                dsb.chromosomeID = damage.chromosomeID;
                dsb.chromosomePosition = damage.chromosomePosition;
                dsbVec.push_back(dsb);                                  // Store in dsbVec the following; {'chromosomeID', chromosomePosition}
            }
        }
    }
    return dsbVec;
}










// Function to return the location in base pairs for a number of single-strand breaks associated with a damage site in a given chromosome and exposure
std::vector<DamageLocation> Karyogram::getSingleStrandBreaks(
    const std::vector<Exposure>& exposures)				// Access exposures to determine the damage locations according to each chromosome ID
{
    std::vector<DamageLocation> ssbVec;					// Initialize the single-strand break vector to store number of SSBs, location in Mbp, and exposure.

    for (const auto& exposure : exposures)				// Loop through each exposure to draw the damages associated per exposure.
    {
        for (const auto& damage : exposure.damages)			// Loop through all the damage entries in a given exposure.
        {
            if (damage.damageType.numSingleBackboneBreaks > 0)		// SSB = single-strand break present, but no double-strand break.
            {
                DamageLocation ssb;					// Instantiate the DamageLocation object ssb for easy access to SSB location and chromosome ID.
                ssb.chromosomeID = damage.chromosomeID;
                ssb.chromosomePosition = damage.chromosomePosition;
                ssb.numSingleStrandBreaks = damage.damageType.numSingleBackboneBreaks;
                ssbVec.push_back(ssb);					// Store in ssbVec the following: {'chromosomeID', 'chromPosition(Mbp)', 'numSSBs'}
            }
        }
    }

    return ssbVec;
}












// Function to convert the damage location stored in the SDD damage block in base pairs to the fractional length of the chromosome.
double Karyogram::getDamageFraction(
    const DamageLocation& damage, 				// Access damage location object to get damage location in Mbp
    double chromosomeSize)					// Obtain chromosomSize to convert damage location to fraction of chromosome length
{
    double position = damage.chromosomePosition.position;

    if (damage.chromosomePosition.isFractional)                 // Check if SDD data field 4 is already in a fractional field format
    {
        return position;                                        // If so, return the fractional position, no conversion needed.
    }

    if (chromosomeSize <= 0.0)
    {
        return 0.0;
    }

    double chromosomeSizeBP = chromosomeSize * 1000000.0;       // p arm = top, q arm = bottom, chromosome sizes stored in Mbp

    return position / chromosomeSizeBP;                         // If SDD data field 4 is not fractional, return scaled position in fractional format.
}

} // namespace karyoke
