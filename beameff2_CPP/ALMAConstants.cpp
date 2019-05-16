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

#include "ALMAConstants.h"
#include <math.h>
#include <iostream>
using namespace std;

namespace ALMAConstants {

const float c = 2.99792458e8;       // m/s
const float TAU = 0.25;

// Note: FOCAL_DEPTH=0.7197 comes from equation 22 of ALMA MEMO 456 using M=20 and phi_0 = 3.58
// TODO: Is it incorrect for ACA 7m antennas?
const float FOCAL_DEPTH(0.7197);

const float DEG_TO_RAD(M_PI / 180.0);
const float RAD_TO_DEG(180.0 / M_PI);

void getAntennaParameters(PointingOptions pointing,
                          float &M, float &psi_o, float &psi_m,
                          float &plateFactor, float &dishDiameter)
{
    switch (pointing) {
    case ACA7METER:
        M = 21.775537595;
        psi_o = 68.4694425916;
        psi_m = 3.5798212165;
        plateFactor = 3.6833;   // arc-s/mm
        dishDiameter = 7000.0; // mm
        break;

    case ACTUAL:
    case NOMINAL:
    case BAND1TEST:
    default:
        M = 20.0;
        psi_o = 64.0154815383723;
        psi_m = 3.5800849111;
        plateFactor = 2.148;;   // arc-s/mm
        dishDiameter = 12000.0; // mm
        break;
    }
}

float getSubreflectorRadius(PointingOptions pointing) {
    const float subreflector_radius12m = 3.57633437;  // degrees
    const float subreflector_radius7m = 3.5798212165; // degrees
    const float subreflector_reduced = 2.0;

    switch (pointing) {
    case ACA7METER:
        return subreflector_radius7m;

    case REDUCE_SUB:
        return subreflector_reduced;

    case ACTUAL:
    case NOMINAL:
    case BAND1TEST:
    default:
        return subreflector_radius12m;
    }
}

bool getNominalAngles(int band, PointingOptions pointing, float &az, float &el) {

    switch (pointing) {
    case ACTUAL:
        az = 0;
        el = 0;
        break;

    case NOMINAL:
        switch (band) {
        case 1:  az = +1.7553; el = +1.7553; break;
        case 2:  az = +1.7553; el = -1.7553; break;
        case 3:  az = -0.3109; el = -1.7345; break;
        case 4:  az = -0.3109; el = +1.7345; break;
        case 5:  az = -1.6867; el = -1.6867; break;
        case 6:  az = -1.6867; el = +1.6867; break;
        case 7:  az = -0.9740; el =  0.0000; break;
        case 8:  az =  0.0000; el = -0.9740; break;
        case 9:  az =  0.0000; el = +0.9740; break;
        case 10: az = +0.9740; el =  0.0000; break;
        default:
            cout << "ERROR: getNominalAngles() Illegal band number = " << band << endl;
            return false;
        }
        break;

    case ACA7METER:
        switch (band) {
        case 1:  az = +2.943499; el = +2.943499; break;
        case 2:  az = +2.898850; el = -2.898850; break;
        case 3:  az = -0.521949; el = -2.918507; break;
        case 4:  az = -0.549947; el = +3.120381; break;
        case 5:  az = -1.874227; el = -1.874227; break;
        case 6:  az = -1.990572; el = +1.990572; break;
        case 7:  az = -0.764417; el =  0.000000; break;
        case 8:  az =  0.000000; el = -0.757841; break;
        case 9:  az =  0.000000; el = +0.735064; break;
        case 10: az = +0.735064; el =  0.000000; break;
        default:
            cout << "ERROR: getNominalAngles() Illegal band number = " << band << endl;
            return false;
        }
        break;

    case BAND1TEST:
        switch (band) {
        case 1:
            az = 0; el = -2.48;
            break;

        default:
            cout << "ERROR: getNominalAngles() Illegal band number for pointing option BAND1TEST = " << band << endl;
            return false;
        }
        break;

    default:
        cout << "ERROR: getNominalAngles(): Illegal value for pointingOption = " << pointing << endl;
        return false;
    }
    return true;
}


}; // namespace

