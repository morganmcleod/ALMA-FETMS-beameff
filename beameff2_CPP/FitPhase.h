#ifndef FITPHASE_H_
#define FITPHASE_H_
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

class ScanData;

namespace BeamFitting {

    void FitPhase(ScanData *currentScan, float azNominal, float elNominal);
    ///< perform phase fit and save results back into currentScan
    ///< @param currentScan: ScanData object containing FF scan for phase fitting
    ///< @param azNominal: nominal pointing angle in Az
    ///< @param elNominal: nominal pointing angle in El

}; // namespace

#endif /* FITPHASE_H_ */
