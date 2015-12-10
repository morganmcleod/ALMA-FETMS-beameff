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

#include "ScanSet.h"
#include "ScanData.h"
#include "ScanDataRaster.h"
#include "ALMAConstants.h"
#include "FitPhase.h"
#include "FitAmplitude.h"
#include "iniparser.h"
#include <math.h>
#include <iostream>
using namespace std;

ScanSet::ScanSet(int id)
  : CopolPol0_m(NULL),
    XpolPol0_m(NULL),
    CopolPol1_m(NULL),
    XpolPol1_m(NULL),
    Copol180_m(NULL)
{
    clear();
    id_m = id;
}

ScanSet::~ScanSet() {
    clear();
}

#define WHACK(P) { if (P) delete P; P = NULL; }

void ScanSet::clear() {
    id_m = 0;
    band_m = 0;
    WHACK(CopolPol0_m)
    WHACK(XpolPol0_m)
    WHACK(CopolPol1_m)
    WHACK(XpolPol1_m)
    WHACK(Copol180_m)

    memset(&Pol0Effs_m, 0, sizeof(EfficiencyData));
    memset(&Pol1Effs_m, 0, sizeof(EfficiencyData));

    nominal_z_offset_m = 0.0;
    squint_m = 0.0;
    squint_arcseconds_m = 0.0;
}

bool ScanSet::loadScan(const dictionary *dict, const std::string inputSection, const std::string delim) {
    bool error = false;

    ScanData *scan = new ScanData;
    // load the contents of the section:
    if (!scan -> loadFromIni(dict, inputSection)) {
        cout << "ERROR: ScanSet::loadScan(): loadFromIni failed." << cout;
        error = true;

    // load the listing files:
    } else if (!scan -> loadListings(delim)) {
        cout << "ERROR: ScanSet::loadScan(): loadListings failed." << cout;
        error = true;

    } else {
        // get the band number from the first scan:
        if (band_m == 0)
            band_m = scan -> getBand();

        // combine dual Z scans:
        scan -> combineDualZScans();

        // determine which slot to store the scan object in:

        // check whether this is a 180-degree copol scan for phase center correction (can be either pol):
        if (scan -> getScanType() == ScanData::COPOL180) {
            if (Copol180_m) {
                cout << "ERROR: ScanSet::loadScan(): already a scan loaded for 180 degree copol." << endl;
                error = true;
            } else
                Copol180_m = scan;

        // check for pol 0 scans:
        } else if (scan -> getPol() == 0) {
            if (scan -> getScanType() == ScanData::COPOL) {
                if (CopolPol0_m) {
                    cout << "ERROR: ScanSet::loadScan(): already a scan loaded for copol pol0." << endl;
                    error = true;
                } else
                    CopolPol0_m = scan;
            }
            if (scan -> getScanType() == ScanData::XPOL) {
                if (XpolPol0_m) {
                    cout << "ERROR: ScanSet::loadScan(): already a scan loaded for xpol pol0." << endl;
                    error = true;
                } else
                    XpolPol0_m = scan;
            }

        // check for pol 1 scans:
        } else if (scan -> getPol() == 1) {
            if (scan -> getScanType() == ScanData::COPOL) {
                if (CopolPol1_m) {
                    cout << "ERROR: ScanSet::loadScan(): already a scan loaded for copol pol1." << endl;
                    error = true;
                } else
                    CopolPol1_m = scan;
            }
            if (scan -> getScanType() == ScanData::XPOL) {
                if (XpolPol1_m) {
                    cout << "ERROR: ScanSet::loadScan(): already a scan loaded for xpol pol1." << endl;
                    error = true;
                } else
                    XpolPol1_m = scan;
            }
        }
    }
    if (error) {
        delete scan;
        return false;
    } else
        return true;
}

bool ScanSet::getEfficiencies(ALMAConstants::PointingOptions pointingOption) {
    bool ret = true;

    // calculate for pol0:
    ret = ret && getEfficiencies_impl(pointingOption, CopolPol0_m, XpolPol0_m, Pol0Effs_m);

    // calculate for pol1:
    ret = ret && getEfficiencies_impl(pointingOption, CopolPol1_m, XpolPol1_m, Pol1Effs_m);

    // calculate average Z offset of both scans:
    if (CopolPol0_m && CopolPol1_m) {
        const ScanDataRaster *pol0FF = CopolPol0_m -> getFFScan();
        const ScanDataRaster *pol1FF = CopolPol1_m -> getFFScan();
        if (pol0FF && pol1FF) {
            const AnalyzeResults &pol0 = pol0FF -> getAnalyzeResults();
            const AnalyzeResults &pol1 = pol1FF -> getAnalyzeResults();
            nominal_z_offset_m = (pol0.deltaZ + pol1.deltaZ) / 2.0;
        }
    }

    // calculate phase 2 for pol0:
    ret = ret && getEfficiencies_impl2(pointingOption, CopolPol0_m, XpolPol0_m, Pol0Effs_m);

    // calculate phase 2 for pol1:
    ret = ret && getEfficiencies_impl2(pointingOption, CopolPol1_m, XpolPol1_m, Pol1Effs_m);

    return ret;
}

bool ScanSet::getEfficiencies_impl(ALMAConstants::PointingOptions pointingOption,
                                   ScanData *copolScan, ScanData *xpolScan, EfficiencyData &eff)
{
    // get the secondary reflector radius for the selected pointing option:
    float subreflectorRadius = ALMAConstants::getSubreflectorRadius(pointingOption);

    // get the nominal angles for the selected pointing option:
    float azNominal, elNominal, copolPeakAmp(0);
    if (!ALMAConstants::getNominalAngles(band_m, pointingOption, azNominal, elNominal))
        return false;

    cout << "azNominal = " << azNominal << endl;
    cout << "elNominal = " << elNominal << endl;

    // find peak and atten differences:
    if (copolScan && xpolScan) {
        // Find the difference in copol and xpol IF attenuation:
        eff.ifatten_diff = copolScan -> getIfAtten() - xpolScan -> getIfAtten();
        // Adjust the xpol data to account for the difference:
        xpolScan -> subtractForAttenuator(eff.ifatten_diff);
    }

    // analyze copol:
    if (copolScan) {
        // for ACTUAL pointing option use the center of mass for analysis:
        copolScan -> calcCenterOfMass();
        if (pointingOption == ALMAConstants::ACTUAL) {
            copolScan -> getFFCenterOfMass(azNominal, elNominal);
            cout << "using azActual = " << azNominal << endl;
            cout << "using elActual = " << elNominal << endl;
        }
        copolScan -> analyzeBeams(azNominal, elNominal, subreflectorRadius);
        copolPeakAmp = copolScan -> getFFPeak();
        cout << "copolPeakAmp = " << copolPeakAmp << endl;
        BeamFitting::FitPhase(copolScan, azNominal, elNominal);
        BeamFitting::FitAmplitude(copolScan, azNominal, elNominal);
    }
    // analyze xpol, using azNominal, elNominal from copol if pointing is ACTUAL:
    if (xpolScan)
        xpolScan -> analyzeBeams(azNominal, elNominal, subreflectorRadius, copolPeakAmp);

    return true;
}

bool ScanSet::getEfficiencies_impl2(ALMAConstants::PointingOptions pointingOption,
                                    const ScanData *copolScan, const ScanData *xpolScan, EfficiencyData &eff)
{
    float M, psi_o, psi_m, lambda, delta, beta;
    ALMAConstants::getAntennaParameters(pointingOption, M, psi_o, psi_m);

    if (copolScan && xpolScan) {
        // calculate amplitude peak differences between copol and xpol:
        eff.peak_diff_FF = copolScan -> getFFPeak() - xpolScan -> getFFPeak();
        eff.peak_diff_NF = copolScan -> getNFPeak() - xpolScan -> getNFPeak();

        // get the analysis results to compute the rest of the efficiencies:
        const ScanDataRaster *copolFF = copolScan -> getFFScan();
        const ScanDataRaster *xpolFF = xpolScan -> getFFScan();

        if (copolFF && xpolFF) {
            const AnalyzeResults &copol = copolFF -> getAnalyzeResults();
            const AnalyzeResults &xpol = xpolFF -> getAnalyzeResults();

            // TICRA method for computing spillover:
            eff.eta_spill_co_cross = (copol.sumPowerOnSec + xpol.sumPowerOnSec) / (copol.sumSqE + xpol.sumSqE);

            // TICRA method ratio of copol beam to copol+xpol on secondary:
            eff.eta_pol_on_secondary = (copol.sumPowerOnSec) / (copol.sumPowerOnSec + xpol.sumPowerOnSec);

            // Total polariazaiton*spillover efficiency:
            eff.eta_pol_spill = eff.eta_spill_co_cross * eff.eta_pol_on_secondary;

            // Spillover efficiency, using the 'alternative' definition from R.Hills paper,
            // is the ratio of copol power on the secondary to total copol power:
            eff.eta_spillover = copol.sumPowerOnSec / copol.sumSqE;

            // Polarization efficiency on the secondary using the 'alternative' definition from R.Hills paper,
            // is the ratio of total copol power to total copol+xpol power, NOT masked for the secondary:
            eff.eta_pol = copol.sumSqE / (copol.sumSqE + xpol.sumSqE);

            // Taper efficiency (Amplitude efficiency in R.Hills paper) is the ratio of the illumination
            // of the secondary to a uniform illumination:
            eff.eta_taper = pow(copol.sumEOnSec, 2.0) / (copol.sumMask * copol.sumPowerOnSec);

            // Total efficiency excluding polarization and defocus:
            eff.deltaX = copol.deltaX;
            eff.deltaY = copol.deltaY;
            eff.deltaZ = copol.deltaZ;
            eff.eta_phase = copol.etaPhase;
            eff.eta_tot_np = eff.eta_phase * eff.eta_spillover * eff.eta_taper;

            // Illumination efficiency is the product of taper and spillover efficiency:
            eff.eta_illumination = eff.eta_taper * eff.eta_spillover;

            // Power in dB at edge of the secondary is computed from the average voltage in the edge region:
            eff.edge_taper_db = 20 * log10(copol.sumEOnEdge / copol.edgeNumPts);

            // Total efficiency excluding defocus:
            eff.eta_tot_nd = eff.eta_tot_np * eff.eta_pol;

            // Defocus efficiency:
            eff.eta_defocus = 1 - 0.3 * pow((copolScan -> getKWaveNumber() * (eff.deltaZ - nominal_z_offset_m) / 1000) * pow(32.0, -2.0), 2.0);

            // Total efficiency including defocus:
            eff.total_aperture_eff = eff.eta_tot_nd * eff.eta_defocus;

            // Defocus calculation to find subreflector shift:
            // Note: 0.7197 comes from equation 22 of ALMA MEMO 456 using M=20 and phi_0 = 3.58
            lambda = ALMAConstants::c_mm_per_ns / copolScan -> getRFGhz();
            delta = (eff.deltaZ - copolScan -> getZDistance() + 0.0000000000001) / pow(M, 2.0) / 0.7197;
            beta = (2 * M_PI / lambda) * delta * (1 - cos(psi_o / 57.3));

            eff.defocus_efficiency =
            100*
            (
                pow(ALMAConstants::tau, 2.0)
                * pow(sin(beta / 2.0) / (beta / 2.0), 2.0)
                + (1 - pow(ALMAConstants::tau, 2.0))
                * (
                    pow(sin(beta / 2.0) / (beta / 2.0), 4.0)
                    + 4.0 / pow(beta, 2.0)
                    * pow((sin(beta) / beta) - 1.0, 2.0)
                  )
            );
            eff.shift_from_focus_mm = eff.deltaZ - copolScan -> getZDistance();
            eff.subreflector_shift_mm = fabs(eff.shift_from_focus_mm) / pow(M, 2.0) / 0.7197;
            eff.mean_subreflector_shift = nominal_z_offset_m / pow(M, 2.0) / 0.7197;
            return true;

        }
    }
    return false;
}


void EfficiencyData::print(int _indent) {
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
    cout << indent << "shift_from_focus_mm = " << shift_from_focus_mm << endl;
    cout << indent << "subreflector_shift_mm = " << subreflector_shift_mm << endl;
    cout << indent << "defocus_efficiency = " << defocus_efficiency << endl;
    cout << indent << "mean_subreflector_shift = " << mean_subreflector_shift << endl;
}

void ScanSet::print(int _indent) {
    string indent(_indent, ' ');

    cout << indent << "ScanSet(" << id_m << "):" << endl;
    cout << indent << "band_m = " << band_m << endl;
    cout << indent << "nominal_z_offset_m = " << nominal_z_offset_m << endl;
    cout << indent << "squint_m = " << squint_m << endl;
    cout << indent << "squint_arcseconds_m = " << squint_arcseconds_m << endl;

    cout << indent << "EfficiencyData(pol0):" << endl;
    Pol0Effs_m.print(_indent + 2);

    cout << indent << "EfficiencyData(pol1):" << endl;
    Pol1Effs_m.print(_indent + 2);

    cout << indent << "ScanData(pol0 copol): ";
    if (!CopolPol0_m)
        cout << indent << "NULL" << endl;
    else {
        cout << indent << endl;
        CopolPol0_m -> print(_indent + 2);
    }
    cout << indent << "ScanData(pol0 xpol): ";
    if (!XpolPol0_m)
        cout << indent << "NULL" << endl;
    else {
        cout << indent << endl;
        XpolPol0_m -> print(_indent + 2);
    }
    cout << indent << "ScanData(pol1 copol): ";
    if (!CopolPol1_m)
        cout << indent << "NULL" << endl;
    else {
        cout << indent << endl;
        CopolPol1_m -> print(_indent + 2);
    }
    cout << indent << "ScanData(pol1 xpol): ";
    if (!XpolPol1_m)
        cout << indent << "NULL" << endl;
    else {
        cout << indent << endl;
        XpolPol1_m -> print(_indent + 2);
    }
    cout << indent << "ScanData(180): ";
    if (!Copol180_m)
        cout << indent << "NULL" << endl;
    else {
        cout << indent << endl;
        Copol180_m -> print(_indent + 2);
    }
}



