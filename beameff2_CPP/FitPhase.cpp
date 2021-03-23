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

#include "FitPhase.h"
#include "ScanData.h"
#include "ScanDataRaster.h"
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <math.h>
using namespace std;

extern "C" {
    #include "nr/nr.h"
}

namespace BeamFitting {
    // forward declare fitness functions:
    void dfunction_phase(float p[], float df[]);
    float function_phase(float p[]);
    float function_phase_z(float delta_z);

    // forward declare debug function:
    void MapPhaseEff(float p[]);

    const int nTerms_m(3);              ///< terms of phase center model
    const float linSearchTol(4.40e-4);  ///< linear search tolerance = about the square root of the machine epsilon for float
    const float ftol(1.19e-7);          ///< function evaluation tolerance = about the machine epsilon for float
    ScanData *fitPhaseScan(NULL);       ///< Scan we are currently fitting against
    static float azNominal;             ///< nominal pointing in Az
    static float elNominal;             ///< nominal pointing in El
    static float delta_x_m;             ///< delta_x to hold constant while searching in delta_z
    static float delta_y_m;             ///< delta_y ""
    static bool approxFit(false);       ///< true=use small-angle approximation phase model for fit; false=use exact model
    bool reduceSubreflector(false);     ///< true=use a reduced size subreflector to search

    void FitPhase(ScanData *currentScan, float _azNominal, float _elNominal) {
        BeamFitting::fitPhaseScan = currentScan;
        BeamFitting::azNominal = _azNominal;
        BeamFitting::elNominal = _elNominal;

        int iter_phase;         // how many iterations taken by frprmn()
        float fret_phase;       // fit residual error returned by frprmn() = 1-eta_phase
        float p[nTerms_m + 1];  // terms of the fit search: [0, x, y, z]

        // start from the probe z distance as our guess for delta_z:
        float k = fitPhaseScan -> getKWaveNumber();  // rad/m
        float zRadians = - fitPhaseScan -> getZDistance() * k / 1000.0;
        p[0] = 0.0;
        p[1] = 0.0;
        p[2] = 0.0;
        p[3] = zRadians;

        cout << "StartPos 0: " << 1000.0 * p[1] / k << " mm, "
                               << 1000.0 * p[2] / k << " mm, "
                               << 1000.0 * p[3] / k << " mm, " << endl;

        // first pass multidimensional search with a reduced subreflector mask:
        reduceSubreflector = true;
        frprmn(p, nTerms_m, ftol, &iter_phase, &fret_phase, &function_phase, &dfunction_phase);

        cout << "FitPhase 1: " << 1000.0 * p[1] / k << " mm, "
                               << 1000.0 * p[2] / k << " mm, "
                               << 1000.0 * p[3] / k << " mm, "
                               << "eta=" << 1.0 - fret_phase << ", "
                               << iter_phase << " iterations<br>" << endl;

        // second pass search with the full subreflector maskL
        reduceSubreflector = false;
        frprmn(p, nTerms_m, ftol, &iter_phase, &fret_phase, &function_phase, &dfunction_phase);

        cout << "FitPhase 2: " << 1000.0 * p[1] / k << " mm, "
                               << 1000.0 * p[2] / k << " mm, "
                               << 1000.0 * p[3] / k << " mm, "
                               << "eta=" << 1.0 - fret_phase << ", "
                               << iter_phase << " iterations<br>" << endl;

        // convert back to mm:
        float deltaX = 1000.0 * p[1] / k;
        float deltaY = 1000.0 * p[2] / k;
        float deltaZ = 1000.0 * p[3] / k;

        // save results into the scan object:
        fitPhaseScan -> setPhaseFitResults(deltaX, deltaY, deltaZ, 1.0 - fret_phase);
    }

    void MapPhaseEff(float p[]) {
        // Was used for troubleshooting the phase fit.
        // Outputs a matrix of eta_phase vs delta_x, delta_y at fixed delta_z
        float k = fitPhaseScan -> getKWaveNumber();  // rad/m
        float delta_x = 1000.0 * p[1] / k;
        float delta_y = 1000.0 * p[2] / k;
        float delta_z = 1000.0 * p[3] / k;
        float searchWin = 0.25;
        float stepSize = 0.02;
        float eta_phase;

        cout << "\ndelta_z = " << delta_z << endl;

        for (float x = delta_x - searchWin; x <= delta_x + searchWin; x += stepSize)
            cout << "\t" << x;
        cout << endl;
        for (float y = delta_y - searchWin; y <= delta_y + searchWin; y += stepSize) {
            cout << y;
            for (float x = delta_x - searchWin; x <= delta_x + searchWin; x += stepSize) {
                p[1] = 0.001 * x * k;
                p[2] = 0.001 * y * k;
                eta_phase = 1.0 - function_phase(p);
                cout << '\t' << eta_phase;
            }
            cout << endl;
        }
        cout << endl;
    }

    /// Compute the gradient of the phase fitting function (1.0-eta_phase) at the given coordinates:
    /// p[] = {delta_x, delta_y, delta_z} in radians.
    /// Gradient is returned in df[].
    void dfunction_phase(float p[], float df[]) {
        // Since it is not a analytic function, we must compute the partial derivatives numerically,
        // using:  d(func) / dp[j] = (func(p[j]+del) - func(p[j])) / del

        int i,j;
        float par[nTerms_m + 1];
        float delta = 0.01, del;
        float tol = 4.40e-4;

        const ScanDataRaster *pscan = fitPhaseScan -> getFFScan();
        if (!pscan)
            return;

        for (j = 1; j <= nTerms_m; j++) {
            /* set all the parameters back to the current values */
            for (i = 1; i <= nTerms_m; i++) {
                par[i] = p[i];
            }
            /* apply a small offset to the parameter being adjusted this time through the loop */
            if (fabs(par[j]) >= tol) {
                del = delta * par[j];
            } else {
                /* this takes care of the unique case where the initial guess is zero */
                del = delta;
            }
            par[j] += del;
            df[j] = ((1.0 - pscan -> calcPhaseEfficiency(par, azNominal, elNominal, approxFit, reduceSubreflector))
                   - (1.0 - pscan -> calcPhaseEfficiency(p, azNominal, elNominal, approxFit, reduceSubreflector))) / del;
        }
    }

    /// Return (1.0-eta_phase) if we fit the phase center to the given coordinates:
    /// p[] = {delta_x, delta_y, delta_z} in radians.
    float function_phase(float p[]) {
        const ScanDataRaster *pscan = fitPhaseScan -> getFFScan();
        if (!pscan)
            return 0.0;

        return (1.0 - pscan -> calcPhaseEfficiency(p, azNominal, elNominal, approxFit, reduceSubreflector));
    }

    /// Return (1 - eta_phase) if we fit the phase center delta_z to the given value in radians,
    /// with delta_x and delta_y from the previous search:
    float function_phase_z(float delta_z) {
        float p[nTerms_m + 1];
        p[0] = 0.0;
        p[1] = delta_x_m;
        p[2] = delta_y_m;
        p[3] = delta_z;
        return function_phase(p);
    }

}; //namespace

