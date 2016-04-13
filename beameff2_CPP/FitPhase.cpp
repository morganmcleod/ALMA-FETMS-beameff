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
#include <math.h>

#include <iostream>
#include <iomanip>
#include <limits>
using namespace std;

extern "C" {
    #include "nr/nr.h"
}

namespace BeamFitting {
    // forward declare fitness functions:
    void dfunction_phase(float p[], float df[]);
    float function_phase(float p[]);
    float function_phase_z(float delta_z);

    const int nTerms_m(3);          ///< terms of phase center model
    ScanData *fitPhaseScan(NULL);   ///< Scan we are currently fitting against
    static float azNominal;         ///< nominal pointing in Az
    static float elNominal;         ///< nominal pointing in El
    static float delta_x_m;         ///< delta_x to hold constant while searching in delta_z
    static float delta_y_m;         ///< delta_y ""
    static bool approxFit(false);   ///< true=use approximate phase fitting; false=exact.

    void FitPhase(ScanData *currentScan, float _azNominal, float _elNominal) {
        fitPhaseScan = currentScan;
        azNominal = _azNominal;
        elNominal = _elNominal;

        // set the linear search tolerance at about the square root of the machine epsilon for float:
        float tol = 4.40e-4;

        // set the function evaluation tolerance at about the machine epsilon for float:
        float ftol = 1.93e-7;

        int iter_phase;
        float fret_phase;
        float p[nTerms_m + 1];

        void (*dfunctionphase)(float p[], float df[]);
        float (*functionphase)(float p[]);

        functionphase = &function_phase;
        dfunctionphase = &dfunction_phase;

        // start from the probe z distance as our guess for delta_z:
        p[0] = 0.0;
        p[1] = 0.0;
        p[2] = 0.0;
        p[3] = fitPhaseScan -> getZDistance();

        // convert initial guess to radians:
        float k = fitPhaseScan -> getKWaveNumber();  // rad/m
        p[1] *= 0.001 * k;                           // rad
        p[2] *= 0.001 * k;
        p[3] *= 0.001 * k;

        // first find the phase center minimum closest to the given Z distance:
        frprmn(p, nTerms_m, ftol, &iter_phase, &fret_phase, functionphase, dfunctionphase);

        cout << "FitPhase 1: " << 1000.0 * p[1] / k << ", "
                               << 1000.0 * p[2] / k << ", "
                               << 1000.0 * p[3] / k << ", "
                               << 1.0 - fret_phase << ", "
                               << iter_phase << " iterations<br>" << endl;

        // now do a line search along the z axis:
        delta_x_m = p[1];
        delta_y_m = p[2];
        float ax = p[3];
        float bx = ax + 2.0;
        float cx, fa, fb, fc;

        // bracket the minimum of the phase fit function in delta_z:
        mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, function_phase_z);

        // now find the exact delta_z giving the minimum:
        fb = brent(ax, bx, cx, function_phase_z, tol, &bx);

        // now search the space around the first delta_x, delta_y and the new delta_z:
        p[1] = delta_x_m;
        p[2] = delta_y_m;
        p[3] = bx;
        frprmn(p, nTerms_m, ftol, &iter_phase, &fret_phase, functionphase, dfunctionphase);

        cout << "FitPhase 2: " << 1000.0 * p[1] / k << ", "
                               << 1000.0 * p[2] / k << ", "
                               << 1000.0 * p[3] / k << ", "
                               << 1.0 - fret_phase << ", "
                               << iter_phase << " iterations<br>" << endl;

        // another bracket and linear search in delta_z:
        delta_x_m = p[1];
        delta_y_m = p[2];
        ax = p[3];
        bx = ax + 2.0;
        mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, function_phase_z);
        fb = brent(ax, bx, cx, function_phase_z, tol, &bx);

        // final multidimensional search:
        p[3] = bx;
        frprmn(p, nTerms_m, ftol, &iter_phase, &fret_phase, functionphase, dfunctionphase);

        cout << "FitPhase 3: " << 1000.0 * p[1] / k << ", "
                               << 1000.0 * p[2] / k << ", "
                               << 1000.0 * p[3] / k << ", "
                               << 1.0 - fret_phase << ", "
                               << iter_phase << " iterations<br>" << endl << endl;

        // convert back to mm:
        float deltaX = 1000.0 * p[1] / k;
        float deltaY = 1000.0 * p[2] / k;
        float deltaZ = 1000.0 * p[3] / k;

        // save results into the scan object:
        fitPhaseScan -> setPhaseFitResults(deltaX, deltaY, deltaZ, 1.0 - fret_phase);
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
            df[j] = ((1.0 - pscan -> calcPhaseEfficiency(par, azNominal, elNominal, approxFit))
                   - (1.0 - pscan -> calcPhaseEfficiency(p, azNominal, elNominal, approxFit))) / del;
        }
    }

    /// Return (1.0-eta_phase) if we fit the phase center to the given coordinates:
    /// p[] = {delta_x, delta_y, delta_z} in radians.
    float function_phase(float p[]) {
        const ScanDataRaster *pscan = fitPhaseScan -> getFFScan();
        if (!pscan)
            return 0.0;

        return (1.0 - pscan -> calcPhaseEfficiency(p, azNominal, elNominal, approxFit));
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

