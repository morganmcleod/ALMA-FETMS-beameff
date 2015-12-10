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

// Forward declare implementation classes:
class ScanData;
struct _dictionary_;
typedef _dictionary_ dictionary;

/// All the efficiency data which is collected for a single polarization:
struct EfficiencyData {

    //power differences (copol-xpol):
    float ifatten_diff;             ///< dB difference in IF attenuation in effect during the measurement
    float peak_diff_NF;             ///< dB difference between the nearfield scans peak power levels
    float peak_diff_FF;             ///< dB difference between the farfield scans peak power levels

    //Copol phase fit results:
    float deltaX, deltaY, deltaZ;   ///< phase center coordinates.
    float eta_phase;                ///< phase efficiency

    //Copol efficiency values
    float eta_spillover;            ///< Spillover efficiency value in [0...1]
    float eta_taper;                ///< Taper efficiency in [0...1]
    float eta_illumination;         ///< Illumination efficiency in [0...1]

    //TICRA Polarization efficiencies
    float eta_spill_co_cross;       ///< component of xpol spillover efficiency?
    float eta_pol_on_secondary;     ///< component of xpol spillover efficiency?
    float eta_pol_spill;            ///< xpol spillover efficiency in [0...1]?

    //Additional efficiencies
    float edge_taper_db;            ///< Average power in dB of the pixels falling at the edge of the subreflector.
    float eta_tot_np;               ///< Efficiency other than polarizaion and defocus in [0..1]
                                    ///<  (eta_phase * eta_spillover * eta_taper)
    float eta_pol;                  ///< Polarization efficiency in [0..1]
    float eta_tot_nd;               ///< Efficiency other than defocus in [0..1] (eta_tot_np * eta_pol)
    float eta_defocus;              ///< Defocus efficency in [0..1]
    float total_aperture_eff;       ///< Overall aperture efficiency in [0..1] (eta_tot_nd * eta_defocus)
    float shift_from_focus_mm;      ///< Difference between delta_z and the nominal probe distance of 200 mm.
    float subreflector_shift_mm;    ///< ?
    float defocus_efficiency;       ///< Defocus efficency in [0..1] calculated in GetAdditionalEfficiencies().
                                    ///<  How different from eta_defocus?
    float mean_subreflector_shift;  ///< ?

    void print(int indent = 0);
    ///< output object state
};



/// A collection of up to 5 ScanData objects in the following roles:
/// Pol0 Copol, Pol0 Xpol, Pol1 Copol, Pol1 Xpol, 180-degree scan (either pol).
class ScanSet {
public:
    ScanSet(int id = 0);
    ///< default construct with scanSet ID

    ~ScanSet();
    ///< destructor

    void clear();
    ///< reset the object to it's just-constructed state.

    bool loadScan(const dictionary *dict, const std::string inputSection, const std::string delim);
    ///< load a scan from the input file into one of the ScanData slots.
    ///< this class determines which slot depending on the properties of the scan
    ///< @param dict: an iniparser dictionary
    ///< @param inputSection: name of the section to load
    ///< @param delim: the delimiter to use when loading listing files.
    ///< @return true if no errors loading and parsing the section.

    bool getEfficiencies(ALMAConstants::PointingOptions pointingOption);
    ///< calculate beam efficiencies, etc. from the loaded data.

    void print(int indent = 0);
    ///< output object state

private:
    int id_m;                   // scanSet ID
    int band_m;                 // ALMA band number (1-10)
    ScanData *CopolPol0_m;      // Data set for pol0 copol
    ScanData *XpolPol0_m;       // Data set for pol0 xpol
    ScanData *CopolPol1_m;      // Data set for pol1 copol
    ScanData *XpolPol1_m;       // Data set for pol1 xpol
    ScanData *Copol180_m;       // Data set for 180-degree copol scan, can be either pol

    EfficiencyData Pol0Effs_m;  ///< efficiency data for pol0
    EfficiencyData Pol1Effs_m;  ///< efficiency data for pol0

    float nominal_z_offset_m;   ///< average of delta_z for Pol0 and Pol1 copol scans
    float squint_m;             ///< Beam squint as percentage of FWHM of the Beam.   Not calculated by this program?
    float squint_arcseconds_m;  ///< Beam squint in arcsecods.   Not calculated by this program?

    bool getEfficiencies_impl(ALMAConstants::PointingOptions pointingOption,
                              ScanData *copol, ScanData *xpol, EfficiencyData &eff);
    ///< analyze beams for a single pol from copol and xpol ScanData.

    bool getEfficiencies_impl2(ALMAConstants::PointingOptions pointingOption,
                               const ScanData *copol, const ScanData *xpol, EfficiencyData &eff);
    ///< calculate efficiencies for a single pol from copol and xpol ScanData.
};


#endif /* SCANSET_H_ */


