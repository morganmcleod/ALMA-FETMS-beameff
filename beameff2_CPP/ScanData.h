#ifndef SCANDATA_H_
#define SCANDATA_H_
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

#include "dictionary.h"
#include "ALMAConstants.h"
#include <math.h>
#include <string>

class ScanDataRaster;

/// \brief contains the data and metadata for a single scan

/// This class contains the data and metadata for a single scan or a pair of
///   scans taken at two Z distances to be combined numerically into a single
///   scan for analysis.
class ScanData {
public:

    /// Mutually exclusive type of a ScanData object:
    enum ScanTypes {
        UNDEFINED,  ///< undefined scan type
        COPOL,      ///< a copol scan
        XPOL,       ///< a crosspol scan
        COPOL180,   ///< a copol scan taken 180 degrees rotated from ether pol0 or pol1
        DUAL        ///< a scan take at 45 degrees between polarizations
    };

    ScanData();
    ///< default constructor

    ~ScanData();
    ///< destructor

    void clear();
    ///< free resources and reset the object to its just-constructed state.

    bool loadFromIni(const dictionary *dict, const std::string inputSection);
    ///< load the specified ini file section
    ///< @param dict: an iniparser dictionary
    ///< @param inputSection: name of the section to load
    ///< @return true if no errors loading and parsing the section.

    bool loadListings(const std::string &delim, ALMAConstants::InvertPhaseOptions invertPhaseOption);
    ///< load the beam scan listing files which were specified in the input file section.
    ///< @param delim: the delimiter to use when parsing the listing files, typically "," or "tab"
    ///< @param invertPhaseOption: option for when to invert phase and rotate FF scans on load
    ///< @return true if no errors parsing the listing files.

    void combineDualZScans();
    ///< If dual Z scans were specified in the input file (NF2 and FF2), numerically combine them.

    void calcCenterOfMass();
    ///< calculate the center of mass of the beam.
    ///< results are stored and can be retrieved with getCenterOfMass

    void getFFCenterOfMass(float &xCenter, float &yCenter) const;
    ///< retrieve the caculated center of mass of the far field beam.
    ///< @param xCenter will contain the x center of mass
    ///< @param yCenter will contain the y center of mass

    void analyzeBeams(ALMAConstants::PointingOptions pointingOption,
                      float azPointing = 0.0, float elPointing = 0.0,
                      float copolPeakAmp = 0.0, bool doUnwrapPhase = false);
    ///< Compute peak amplitude and phase, NF and FF beam centers of mass, etc.
    ///< @param pointingOption: to use for calculating efficiencies
    ///< @param azPointing, elPointing: if noth non-zero use for calculating efficiencies;
    ///<                                else use nominal or actual from this beam.
    ///< @param copolPeakAmp: If non-zero, use to normalize the FF beam xpol to copol peak.
    ///< @param doUnwrapPhase: If true, unwrap the FF phase

    const std::string &getInputSection() const
      { return inputSection_m; }
    ///< @return the input file section from which this scan was loaded

    const std::string &getMeasDateTime() const
      { if (!measDateTime_m.empty())
          return measDateTime_m;
        else
          return NSIDateTime_m; }
    ///< @return measurement date/time string

    unsigned getScanId() const
      { return scanId_m; }
    ///< @return the scan ID

    int getBand() const
      { return band_m; }
    ///< @return the ALMA band number

    ScanTypes getScanType() const
      { return scanType_m; }
    ///< @return the type of scan: COPOL, XPOL, etc.

    std::string getScanTypeString() const;
    ///< @return the type of scan as a string.

    int getPol() const
      { return pol_m; }
    ///< @return the polarization (0 or 1)

    int getTilt() const
      { return tilt_m; }
    ///< @return the tilt angle of the measurement

    float getZDistance() const
      { return zDistance_m; }
    ///< @return the probe z distance

    float getIfAtten() const
      { return ifAtten_m; }
    ///< @return the IF attenuatom

    float getFFPeak(float *phase = NULL) const;
    ///< @return the peak power and optionally phase at peak in the farfield scan
    ///< @param phase: if not NULL, where to put the phase value: radians.

    float getNFPeak(float *phase = NULL) const;
    ///< @return the peak power and optionally phase at peak in the nearfield scan
    ///< @param phase: if not NULL, where to put the phase value: radians.

    float getRFGhz() const
      { return rfGHz_m; }
    ///< @return the RF

    float getKWaveNumber() const
      { return 2 * M_PI * rfGHz_m * 1.0e9 / ALMAConstants::c; }
    ///< @return k, the wavenumber in rad/m computed from rfGHz_m

    const ScanDataRaster *getFFScan() const
      { return FF_m; }
    ///< @return const pointer to the FF scan data

    const ScanDataRaster *getNFScan() const
      { return NF_m; }
    ///< @return const pointer to the NF scan data

    void setPhaseFitResults(float deltaX, float deltaY, float deltaZ, float etaPhase);
    ///< Store the phase center and efficiency computed in FitPhase()

    void setAmpFitResults(float ampFitAmp, float ampFitWidthDeg, float ampFitUOffDeg,
                          float ampFitVOffDeg, float ampFitD_0_90, float ampFitD_45_135);
    ///< Store the amplitude fit results computed in FitAmplitude()

    void subtractForAttenuator(float ifAttenDiff);
    ///< adjust all NF an FF points for the difference in attenuation.
    ///< @param ifAttenDiff: copol - xpol attenuator setting in dB.

    void print(int indent = 0) const;
    ///< Ouput object state.
    ///< @param indent number of spaces

private:

    std::string inputSection_m; ///< section name in input file
    unsigned scanId_m;          ///< keyID from ScanDetails passed in for plot label

    std::string measDateTime_m; ///< date/time string when measured.
    std::string NSIDateTime_m;  ///< date/time string from NSI text file.

    std::string filenameNF_m;   ///< nearfield input file
    std::string filenameFF_m;   ///< farfield input file
    std::string filenameNF2_m;  ///< 2nd Z distance nearfield input file
    std::string filenameFF2_m;  ///< 2nd Z distance farfield input file

    ScanTypes scanType_m;       ///< UNDEFINED, COPOL, XPOL, DUAL
    double rfGHz_m;             ///< RF in GHz
    int band_m;                 ///< ALMA band in 1-10
    int pol_m;                  ///< -1=undefined, 0 or 1
    int sb_m;                   ///< -1=undefined, 1=USB, 2=LSB
    int tilt_m;                 ///< tilt angle in degrees.
    float zDistance_m;          ///< Z distance of scan probe in mm
    float ifAtten_m;            ///< IF processor attenuation in dB for this single scan.

    unsigned startrowNF_m;      ///< number of header rows to skip when reading filenameNF
    unsigned startrowFF_m;      ///< number of header rows to skip when reading filenameFF
    unsigned startrowNF2_m;     ///< rows to skip in filenameNF2.
    unsigned startrowFF2_m;     ///< rows to skip in filenameFF2.

    ScanDataRaster *NF_m;       ///< Input file raster data for NF
    ScanDataRaster *FF_m;       ///< Input file raster data for FF
    ScanDataRaster *NF2_m;      ///< Input file raster data for NF2.  Null if not specified.
    ScanDataRaster *FF2_m;      ///< Input file raster data for FF2.  Null if not specified.

    void freeArrays();
    ///< private helper to free raster array memory.

    void combineDualZScans_impl(ScanDataRaster &z1, ScanDataRaster &z2);
    ///< private helper to numerically combine dual Z scans.
    ///< @param z1 the first scan.  Will be replaced by combined scan.
    ///< @param z2 the second scan.
};

#endif /* SCANDATA_H_ */
