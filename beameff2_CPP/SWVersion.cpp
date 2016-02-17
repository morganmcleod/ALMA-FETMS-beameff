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

#include "SWVersion.h"

const std::string BEAMEFF_SW_VERSION_STRING("2.0.1");

// Version history
// 2.0.1:  Name pointing angle plots for RF if scanset_id not provided
//         Name [results_xx] section for RF if scanset_id not provided
//         Use el = -2.48 for band1test dewar
//         Calculate squint even when 180-degree scan for correcting probe asymmetry is not provided.
//         For backwards compatibility with LabVIEW wrapper, squint results in output section for pol1, copol.
// 2.0.0:  First release of BeamEff2 C++


