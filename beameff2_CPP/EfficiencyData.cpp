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

#include "EfficiencyData.h"
#include <iostream>
#include <string>
using namespace std;

namespace EfficiencyData {

    void CombinedPols::clear() {
        nominal_z_offset = 0.0;
        x_corr = y_corr = 0.0;
        dist_between_centers_mm = 0.0;
        corrected_pol = 0;
        squint_arcseconds = 0.0;
        squint_percent = 0.0;
        FFPointingPlot.clear();
    }

    void CombinedPols::print(int _indent) {
        string indent(_indent, ' ');
        cout << indent << "nominal_z_offset = " << nominal_z_offset << endl;
        cout << indent << "phase center correction = " << x_corr << ", " << y_corr << endl;
        cout << indent << "dist_between_centers_mm = " << dist_between_centers_mm << endl;
        cout << indent << "corrected_pol = " << corrected_pol << endl;
        cout << indent << "squint_arcseconds = " << squint_arcseconds << endl;
        cout << indent << "squint_percent = " << squint_percent << endl;
        cout << indent << "FFPointingPlot = '" << FFPointingPlot << "'" << endl;
    }

    void OnePol::clear() {
        ifatten_diff = 0.0;
        peak_diff_NF = 0.0;
        peak_diff_FF = 0.0;
        deltaX = deltaY = deltaZ = 0.0;
        eta_phase = 0.0;
        correctedX = correctedY = 0.0;
        eta_spillover = 0.0;
        eta_taper = 0.0;
        eta_illumination = 0.0;
        eta_spill_co_cross = 0.0;
        eta_pol_on_secondary = 0.0;
        eta_pol_spill = 0.0;
        edge_taper_db = 0.0;
        eta_tot_np = 0.0;
        eta_pol = 0.0;
        eta_tot_nd = 0.0;
        eta_defocus = 0.0;
        total_aperture_eff = 0.0;
        FFCopolAmpPlot.clear();
        FFCopolPhasePlot.clear();
        NFCopolAmpPlot.clear();
        NFCopolPhasePlot.clear();
        FFCopolPhaseUnwrappedPlot.clear();
        FFCopolPhaseModelPlot.clear();
        FFXpolAmpPlot.clear();
        FFXpolPhasePlot.clear();
        NFXpolAmpPlot.clear();
        NFXpolPhasePlot.clear();
    }

    void OnePol::print(int _indent) {
        string indent(_indent, ' ');
        cout << indent << "ifatten_diff = " << ifatten_diff << endl;
        cout << indent << "peak_diff_NF = " << peak_diff_NF << endl;
        cout << indent << "peak_diff_FF = " << peak_diff_FF << endl;
        cout << indent << "phase center = " << deltaX << ", " << deltaY << ", " <<  deltaZ << endl;
        cout << indent << "eta_phase = " << eta_phase << endl;
        cout << indent << "eta_spillover = " << eta_spillover << endl;
        cout << indent << "eta_taper = " << eta_taper << endl;
        cout << indent << "eta_illumination = " << eta_illumination << endl;
        cout << indent << "eta_spill_co_cross = " << eta_spill_co_cross << endl;
        cout << indent << "eta_pol_on_secondary = " << eta_pol_on_secondary << endl;
        cout << indent << "eta_pol_spill = " << eta_pol_spill << endl;
        cout << indent << "edge_taper_db = " << edge_taper_db << endl;
        cout << indent << "eta_tot_np = " << eta_tot_np << endl;
        cout << indent << "eta_pol = " << eta_pol << endl;
        cout << indent << "eta_tot_nd = " << eta_tot_nd << endl;
        cout << indent << "eta_defocus = " << eta_defocus << endl;
        cout << indent << "total_aperture_eff = " << total_aperture_eff << endl;
        cout << indent << "FFCopolAmpPlot = '" << FFCopolAmpPlot << "'" << endl;
        cout << indent << "FFCopolPhasePlot = '" << FFCopolPhasePlot << "'" << endl;
        cout << indent << "NFCopolAmpPlot = '" << NFCopolAmpPlot << "'" << endl;
        cout << indent << "NFCopolPhasePlot = '" << NFCopolPhasePlot << "'" << endl;
        cout << indent << "FFCopolPhaseUnwrappedPlot = '" << FFCopolPhaseUnwrappedPlot << "'" << endl;
        cout << indent << "FFCopolPhaseModelPlot = '" << FFCopolPhaseModelPlot << "'" << endl;
        cout << indent << "FFXpolAmpPlot = '" << FFXpolAmpPlot << "'" << endl;
        cout << indent << "FFXpolPhasePlot = '" << FFXpolPhasePlot << "'" << endl;
        cout << indent << "NFXpolAmpPlot = '" << NFXpolAmpPlot << "'" << endl;
        cout << indent << "NFXpolPhasePlot = '" << NFXpolPhasePlot << "'" << endl;
    }
}; //namespace
