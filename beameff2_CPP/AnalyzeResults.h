#ifndef ANALYZERESULTS_H_
#define ANALYZERESULTS_H_
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

struct AnalyzeResults {
    // values calculated by calcPeakAndPhase(), calcCenterOfMass(), and analyzeBeam()
    float maxAmp;               ///< maximum amplitude seen
    float phaseAtPeak;          ///< phase at the amplitude peak.
    float xCenterOfMass;        ///< center of mass x coordinate.
    float yCenterOfMass;        ///< center of mass y coordinate.
    float xStepSize;            ///< step size in x coordinate
    float yStepSize;            ///< step size in y coordinate
    double sumMask;             ///< sum of the secondary reflector mask
    double sumE;                ///< sum of E-field voltages, whole scan
    double sumSqE;              ///< sum of squares of E-field voltages, whole scan
    double sumEOnSec;           ///< sum of E-field voltages on the secondary reflector
    double sumSqEOnSec;         ///< sum of squares of E-field voltages on the secondary reflector
    double sumPowerOnSec;       ///< sum of power on the secondary
    double sumSqPowerOnSec;     ///< sum of squares of power on the secondary
    double sumEOnEdge;          ///< sum of voltage in the secondary reflector edge region.
    unsigned long edgeNumPts;   ///< number of points falling in the secondary reflector edge region.

    // FitPhase results:
    float deltaX, deltaY, deltaZ;
                                ///< phase center coordinates.
    float etaPhase;             ///< phase efficiency

    // FitAmplitude results:
    float ampFitAmp;            ///< amplitude
    float ampFitWidthDeg;       ///< width (deg)
    float ampFitUOffDeg;        ///< u_offset (deg)
    float ampFitVOffDeg;        ///< v_offset (deg)
    float ampFitD_0_90;         ///< D_0-90
    float ampFitD_45_135;       ///< D_45-135

    void print(int indent = 0) const;
    ///< output object state
    ///< @param indent number of spaces

};


#endif /* ANALYZERESULTS_H_ */
