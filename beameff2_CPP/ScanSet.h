#ifndef SCANSET_H_
#define SCANSET_H_
/*******************************************************************************
* ALMA - Atacama Large Millimeter Array
* (c) Associated Universities Inc., 2015
*
*This library is free software; you can redistribute it and/or
*modify it under the terms of the GNU Lesser General Public
*License as published by the Free Software Foundation; either
*version 2.1 of the License, or (at your option) any later version.
*
*This library is distributed in the hope that it will be useful,
*but WITHOUT ANY WARRANTY; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*Lesser General Public License for more details.
*
*You should have received a copy of the GNU Lesser General Public
*License along with this library; if not, write to the Free Software
*Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
*
*
*/

#include <string>
#include "ALMAConstants.h"
#include "EfficiencyData.h"

// Forward declare implementation classes:
class ScanData;
struct _dictionary_;
typedef _dictionary_ dictionary;

/// A collection of up to 5 ScanData objects in the following roles:
/// Pol0 Copol, Pol0 Xpol, Pol1 Copol, Pol1 Xpol, 180-degree scan (either pol).
class ScanSet {
public:
    ScanSet(unsigned scanSet = 0);
    ///< default construct
    ///< @param scanSet number from ini file.  Not a database key.

    ~ScanSet();
    ///< destructor

    void clear();
    ///< reset the object to it's just-constructed state.

    void setDatabaseKeys(unsigned scanSetId, unsigned FEConfigId, unsigned TestDataHeaderId);
    ///< store the additional database keys which identify this ScanSet
    ///< @param scanSetId database key to ScanSetDetails table.
    ///< @param FEConfigId database key to FE_Config table.
    ///< @param TestDataHeaderId database key to TestData_Header table.

    bool loadScan(const dictionary *dict, const std::string inputSection,
                  const std::string delim, ALMAConstants::InvertPhaseOptions invertPhaseOption);
    ///< load a scan from the input file into one of the ScanData slots.
    ///< this class determines which slot depending on the properties of the scan
    ///< @param dict: an iniparser dictionary
    ///< @param inputSection: name of the section to load
    ///< @param delim: the delimiter to use when loading listing files.
    ///< @param invertPhaseOption: option for when to invert phase and rotate FF scans on load
    ///< @return true if no errors loading and parsing the section.

    bool calcEfficiencies(ALMAConstants::PointingOptions pointingOption, bool unwrapPhaseOption);
    ///< calculate beam efficiencies, etc. from the loaded data.
    ///< @return true if no errors during calculation.

    bool writeOutputFile(dictionary *dict, const std::string outputFilename);
    ///< save the calcualted efficiency results into the given dictionary and then to the specified output file name
    ///< @return true if no errors writing the results

    bool makePlots(const std::string &outputDirectory, const std::string &gnuplotPath);
    ///< produce the output amplitude, phase, and pointing angle plots
    ///< @param outputDirectory: where to create the plots and temporary data files
    ///< @param gnuplotPath: full path to the system gnuplot command, from the input file.
    ///< @return true if no errors occurred producing the plot.  TODO: doesn't trap gnuplot errors.

    void print(int indent = 0);
    ///< output object state

private:
    unsigned scanSet_m;         ///< 'scanset' number from the input file.  Not related to any database table.
    unsigned scanSetId_m;       ///< keyID from ScanSetDetails table for plot label
    unsigned FEConfigId_m;      ///< keyID from FE_Config table
    unsigned TestDataHeaderId_m;///< keyID from TestData_Header table

    int band_m;                 // ALMA band number (1-10)
    ALMAConstants::PointingOptions pointingOption_m;
                                // pointing option passed to calcEfficiencies()
    bool unwrapPhaseOption_m;   // unwrap phase option passed to analyzeBeams()

    ScanData *CopolPol0_m;      // Data set for pol0 copol
    ScanData *XpolPol0_m;       // Data set for pol0 xpol
    ScanData *CopolPol1_m;      // Data set for pol1 copol
    ScanData *XpolPol1_m;       // Data set for pol1 xpol
    ScanData *Copol180_m;       // Data set for 180-degree copol scan, can be either pol

    EfficiencyData::OnePol Pol0Eff_m;           ///< efficiency data for pol0
    EfficiencyData::OnePol Pol1Eff_m;           ///< efficiency data for pol0
    EfficiencyData::CombinedPols CombinedEff_m; ///< efficiency data for both pols combined

    bool calcEfficiencies_impl(ScanData *copol, ScanData *xpol, EfficiencyData::OnePol &eff);
    ///< Phase 1 calculate efficiencies for a single pol from copol and xpol ScanData.

    bool calcEfficiencies_impl2(const ScanData *copol, const ScanData *xpol, EfficiencyData::OnePol &eff);
    ///< Phase 2 calculate efficiencies for a single pol from copol and xpol ScanData.

    bool analyzeCopol_impl(ScanData *copol, float &azPointing, float &elPointing);
    ///< implementation of copol analysis, reused for copol and copol180 data sets.

    bool calcSquint_impl();
    ///< calculate beam squint for both polarizations plus a third 180-degree scan.

    bool updateCopolSection(dictionary *dict, const std::string &section, const ScanData *copol, const EfficiencyData::OnePol &eff);
    ///< helper to update the dictionary for a copol scan

    bool updateXpolSection(dictionary *dict, const std::string &section, const ScanData *xpol, const EfficiencyData::OnePol &eff);
    ///< helper to update the dictionary for a xpol scan

    static void updateDictionary(dictionary *dict, const std::string &section, const std::string &key, const std::string &value);
    ///< private helper to update the dictionary

    bool writePlotDataFiles(const std::string &outputDirectory, const ScanData *scan,
                            std::string &fileNameFF, std::string &fileNameNF,
                            float copolPeakAmpFF, float copolPeakAmpNF);
    ///< private helper to write output files for plotting
    ///< filenames are returned in fileNameFF and fileNameNF

    bool makeOnePlot(const std::string &outputDirectory, const std::string &gnuplotPath,
                     const std::string &dataFilename, const ScanData *scan,
                     bool nf, bool phase, bool unwrapped,
                     float azPointing, float elPointing, std::string &fileNamePlot);
    ///< private helper to produce an amplitude or phase plot, FF or NF
    ///< prites a Gnuplot command file then executes it
    ///< @param outputDirectory where the plot should be created
    ///< @param gnuplotPath path to the Gnuplot excutable
    ///< @param dataFilename text file containing the raster data to plot
    ///< @param scan ScanData object holding all the metadata to include in the plot
    ///< @param nf true if this is a nearfield plot, false if farfield
    ///< @param phase true is this is a phase plot, false if amplitude
    ///< @param unwrapped true is plot phase unwrapped, false if as-measured
    ///< @param azPointing, elPointing where to center the subreflector circle on the FF plot
    ///< @param fileNamePlot out: the filename of the plot produced
    ///< @return true if no errors setting up the plot.  Won't catch Gnuplot errors.

    bool makePhaseFitPlot(const std::string &outputDirectory, const std::string &gnuplotPath,
                          const ScanData *scan,
                          float azPointing, float elPointing, std::string &fileNamePlot);
    ///< private helper to produce a plot of the FF theoretical phase found by FitPhase
    ///< writes a Gnuplot command file then executes it
    ///< Parameters are a subset of those to makeOnePlot()

    bool makePointingAnglesPlot(const std::string &outputDirectory, const std::string &gnuplotPath,
                                std::string &fileNamePlot);
    ///< private helper to make the pointing angles plot showing both FF copol scans in this ScanSet.

    const std::string &getMeasInfoLabel(std::string &toFill, const ScanData &scan) const;
    ///< private helper to make the measurement info label

    const std::string &getDatabaseKeysLabel(std::string &toFill, unsigned scanId = 0) const;
    ///< private helper to make the database keys label
};


#endif /* SCANSET_H_ */


