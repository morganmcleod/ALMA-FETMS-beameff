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

#include "AnalyzeResults.h"
#include <iostream>
using namespace std;

void AnalyzeResults::print(int _indent) const {
    string indent(_indent, ' ');

    cout << indent << "maxAmp = " << maxAmp << endl;
    cout << indent << "phaseAtPeak = " << phaseAtPeak << endl;
    cout << indent << "xCenterOfMass = " << xCenterOfMass << endl;
    cout << indent << "yCenterOfMass = " << yCenterOfMass << endl;
    cout << indent << "xStepSize = " << xStepSize << endl;
    cout << indent << "yStepSize = " << yStepSize << endl;
    cout << indent << "sumMask = " << sumMask << endl;
    cout << indent << "sumE = " << sumE << endl;
    cout << indent << "sumSqE = " << sumSqE << endl;
    cout << indent << "sumEOnSec = " << sumEOnSec << endl;
    cout << indent << "sumSqEOnSec = " << sumSqEOnSec << endl;
    cout << indent << "sumPowerOnSec = " << sumPowerOnSec << endl;
    cout << indent << "sumSqPowerOnSec = " << sumSqPowerOnSec << endl;
    cout << indent << "sumEOnEdge = " << sumEOnEdge << endl;
    cout << indent << "edgeNumPts = " << edgeNumPts << endl;

    cout << indent << "phaseCenter = " << deltaX << ", " << deltaY << ", " <<  deltaZ << endl;
    cout << indent << "etaPhase = " << etaPhase << endl;
    cout << indent << "ampFitAmp = " << ampFitAmp << endl;
    cout << indent << "ampFitWidthDeg = " << ampFitWidthDeg << endl;
    cout << indent << "ampFitUOffDeg = " << ampFitUOffDeg << endl;
    cout << indent << "ampFitVOffDeg = " << ampFitVOffDeg << endl;
    cout << indent << "ampFitD_0_90 = " << ampFitD_0_90 << endl;
    cout << indent << "ampFitD_45_135 = " << ampFitD_45_135 << endl;

}



