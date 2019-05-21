#ifndef EFFICIENCYDATA_H_
#define EFFICIENCYDATA_H_
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

/// Data structures for calculated beam efficiency results
namespace EfficiencyData {

    /// These results apply to both polarizations together:
    struct CombinedPols {
        float nominal_z_offset;             ///< average of delta_z for Pol0 and Pol1 copol scans
        float x_corr, y_corr;               ///< Phase center correction to apply to the "second" or "(x90,y90)" polarization.
        float dist_between_centers_mm;      ///< Corrected distance between the phase centers
        int corrected_pol;                  ///< Which polarization was corrected.   -1=NONE
        float squint_arcseconds;            ///< Beam squint in arcsecods.
        float squint_percent;               ///< Beam squint as percentage of FWHM the Beam.

        std::string FFPointingPlot;         ///< Pointing angles plot

        CombinedPols()
          { clear(); }
        ///< default constructor.

        void clear();
        ///< reset to initialized state.

        void print(int indent = 0);
        ///< output object state
    };

    /// All the efficiency data which is calculated for a single polarization:
    struct OnePol {
        //power differences (copol-xpol):
        float ifatten_diff;             ///< dB difference in IF attenuation in effect during the measurement
        float peak_diff_NF;             ///< dB difference between the nearfield scans peak power levels
        float peak_diff_FF;             ///< dB difference between the farfield scans peak power levels

        //Copol phase fit results:
        float deltaX, deltaY, deltaZ;   ///< phase center coordinates.
        float eta_phase;                ///< phase efficiency
        float correctedX, correctedY;   ///< corrected phase center, when there is a 180 degree scan available

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
//        float shift_from_focus_mm;      ///< Difference between delta_z and the nominal probe distance.
//        float subreflector_shift_mm;    ///< Distance subreflector must move in order to get best defocus efficiency.
//        float defocus_efficiency;       ///< Defocus efficency in [0..1] assuming a subreflector shift.

        //Plot filenames
        std::string FFCopolAmpPlot;
        std::string FFCopolPhasePlot;
        std::string NFCopolAmpPlot;
        std::string NFCopolPhasePlot;
        std::string FFCopolPhaseUnwrappedPlot;
        std::string FFCopolPhaseModelPlot;
        std::string FFXpolAmpPlot;
        std::string FFXpolPhasePlot;
        std::string NFXpolAmpPlot;
        std::string NFXpolPhasePlot;

        OnePol()
          { clear(); }
        ///< default constructor.

        void clear();
        ///< reset to initialized state.

        void print(int indent = 0);
        ///< output object state
    };
}; // namespace

#endif /* EFFICIENCYDATA_H_ */
