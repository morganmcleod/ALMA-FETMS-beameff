#ifndef ALMACONSTANTS_H_
#define ALMACONSTANTS_H_
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

namespace ALMAConstants {

    extern const float c;
    extern const float TAU;
    extern const float FOCAL_DEPTH;

    /// Different pointing options can be used for efficiency calculation:
    enum PointingOptions {
        NOMINAL,    // the direction of the subreflector
        ACTUAL,     // the actual beam pointing direction, even if not centered on the subreflector
        ACA7METER,  // the nominal direction of the ACA 7 meter antenna.  Also controls other constants.
        BAND1TEST   // the band 1 test dewar
    };

    /// Options for when to invert measured phase and rotate FF scan:
    enum InvertPhaseOptions {
        INVERT_LSB, // invert phase and rotate FF scan when in LSB, the FETMS default.
        INVERT_USB, // invert phase and rotate FF scan when in USB.
        INVERT_ALL, // always invert phase and rotate FF scan.
        INVERT_NONE // never invert phase nor rotate FF scan.
    };

    void getAntennaParameters(PointingOptions pointing,
                              float &M, float &psi_o, float &psi_m,
                              float &plateFactor, float &dishDiameter);
    ///< get the parameters for the selected antenna type.
    ///< @param pointing: from enum PointingOptions
    ///< @param M out
    ///< @param psi_o out
    ///< @param psi_m out
    ///< @param plateFactor out
    ///< @param dishDiameter out

    extern float getSubreflectorRadius(PointingOptions pointing);
    ///< @return the subreflector radius in degrees
    ///< @param pointing: from enum PointingOptions

    extern bool getNominalAngles(int band, PointingOptions pointing, float &az, float &el);
    ///< get the nominal ponting angle for the given band and pointing option.
    ///< @param band: the ALMA band number in 1-10
    ///< @param pointing: from enum PointingOptions
    ///< @param az: the nominal azimuth pointing is returned here.
    ///< @param el: the nominal elevation pointing is returned here.
    ///< @return: true if inputs were valid.

}; // namespace

#endif /* ALMACONSTANTS_H_ */
