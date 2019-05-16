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

const std::string BEAMEFF_SW_VERSION_STRING("2.0.10.2");

// Version history
// 2.0.10: Testing phase unwrapping:
//         .1 using REDUCED_SUB with wrapped phase.   Folders output2.10.1
//         .2 using default sub with wrapped phase.   Folders output2.10.2
//         .3 using default sub with unwrapped phase. Folders output2.10.3
// 2.0.9:  Removed outputs shift_from_focus_mm, subreflector_shift_mm, and defocus_efficiency_due_to_moving_the_subreflector.
//         Their calculation method was broken and of dubious value.
// 2.0.8:  Display actual pointing angles on the pointing plot.
//
// 2.0.7:  Added input file option [settings] invertphase with options {lsb, usb, all, none}
//         and synonymns {yes, y, t, 1} for all, {no, n, f, 0} for none.  All options case-insensitive.
// 2.0.6:  Using multi-pass search with exact expression for fitted phase:
//         1) CG search from {0, 0, probe-zdistance}
//         2) Line search along delta_z
//         3) CG search, line search, final CG search
//         Fix calculation error in squint_percent
// 2.0.5:  Using exact expression for fitted phase including corrections from Alvaro Gonzalez
//         Invert phase and rotate for LSB now.
//         Invert signs of nominal pointing angles.  Matches scanner coord system now.
//         Added makePhaseFitPlot for visualizing the fitted phase surface.
//         Improved documentation of FitPhase and ScanDataRaster::calcPhaseEfficiency
// 2.0.4:  Don't rotate nearfield scans but do invert phase when appropriate.
//         Back to el = +2.48 for band1test dewar
// 2.0.3:  Adjust labels on pointing angles plot
// 2.0.2:  Added makefile and minor fixes for compilation on Linux
// 2.0.1:  Name pointing angle plots for RF if scanset_id not provided
//         Name [results_xx] section for RF if scanset_id not provided
//         Use el = -2.48 for band1test dewar
//         Calculate squint even when 180-degree scan for correcting probe asymmetry is not provided.
//         For backwards compatibility with LabVIEW wrapper, squint results in output section for pol1, copol.
// 2.0.0:  First release of BeamEff2 C++


