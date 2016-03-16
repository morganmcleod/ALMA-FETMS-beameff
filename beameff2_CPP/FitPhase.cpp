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

extern "C" {
    #include "nr/nr.h"
}

namespace BeamFitting {
    // forward declare fitness functions:
    float function_phase(float p[]);
    void dfunction_phase(float p[], float df[]);

    const int nTerms_m(3);          ///< terms of phase center model
    ScanData *fitPhaseScan(NULL);   ///< Scan we are currently fitting against
    static float azNominal;         ///< nominal pointing in Az
    static float elNominal;         ///< nominal pointing in El

    void FitPhase(ScanData *currentScan, float _azNominal, float _elNominal) {
        fitPhaseScan = currentScan;
        azNominal = _azNominal;
        elNominal = _elNominal;

        float ftol = pow(10, -1.2 * fitPhaseScan -> getBand());

        int iter_phase;
        float fret_phase;
        float p[nTerms_m + 1];

        void (*dfunctionphase)(float p[], float df[]);
        float (*functionphase)(float p[]);

        functionphase = &function_phase;
        dfunctionphase = &dfunction_phase;

        // TODO? better way to FitPhase:
        //  first call with p = {0, 0, 0} and mask radius= 1 deg.
        //  find phase center in rad/deg.
        //  change mask to full subreflector ~3.8 deg.
        //  fimd phase center again using P from first iteration.

        p[0] = 0.0;
        p[1] = 0.0;
        p[2] = 0.0;
        p[3] = fitPhaseScan -> getZDistance();

        // convert initial guess to radians:
        float k = fitPhaseScan -> getKWaveNumber();  // rad/m
        p[1] *= 0.001 * k;                           // rad
        p[2] *= 0.001 * k;
        p[3] *= 0.001 * k;

        frprmn(p, nTerms_m, ftol, &iter_phase, &fret_phase, functionphase, dfunctionphase);

        // convert back to mm:
        float deltaX = p[1] / k * 1000;
        float deltaY = p[2] / k * 1000;
        float deltaZ = p[3] / k * 1000;

        // save results into the scan object:
        fitPhaseScan -> setPhaseFitResults(deltaX, deltaY, deltaZ, 1.0 - fret_phase);
    }

    void dfunction_phase(float p[], float df[]) {
    /* This function should compute the gradients of the chi-squared function,
     * which are stored in the array "df". Since it is not a analytic function,
     * we must compute the partial derivatives numerically, which is done using:
     *    d(chi-square) / dp[j] = *chiSquare(p[j]+del) - chiSquare(p[j])) / del
     *
     */
        int i,j;
        float par[nTerms_m + 1];
        float delta = 0.01, del;

        const ScanDataRaster *pscan = fitPhaseScan -> getFFScan();
        if (!pscan)
            return;

        for (j = 1; j <= nTerms_m; j++) {
            /* set all the parameters back to the current values */
            for (i = 1; i <= nTerms_m; i++) {
                par[i] = p[i];
            }
            /* apply a small offset to the parameter being adjusted this time through the loop */
            if (fabs(par[j]) > (delta / 100000)) {
                del = delta * par[j];
            } else {
                /* this takes care of the unique case where the initial guess is zero */
                del = delta;
            }
            par[j] += del;
            df[j] = ((1.0 - pscan -> calcPhaseEfficiency(par, azNominal, elNominal))
                   - (1.0 - pscan -> calcPhaseEfficiency(p, azNominal, elNominal))) / del;
        }
    }

    float function_phase(float p[]) {
    /* This function should return the chi-squared value for the current model
     * parameters that are passed in as an argument.  In this case, we want
     * to minimize the phase "inefficiency", so we simply calculate that.
     */
        const ScanDataRaster *pscan = fitPhaseScan -> getFFScan();
        if (!pscan)
            return 0.0;

        return (1.0 - pscan -> calcPhaseEfficiency(p, azNominal, elNominal));
    }

}; //namespace

