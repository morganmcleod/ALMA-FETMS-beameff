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

const std::string BEAMEFF_SW_VERSION_STRING("2.0.17.1");

//*******                        Be sure to update resource.rc!

// Version history
// 2.0.17.1: Added 5 mrad circle to pointing angles plot.
// 2.0.17:   Added comments to all nr files describing their origin and purpose.
//           Added nr/README-nr.txt
//           64-bit Build using msys64-mingw64 GCC 10.3.0
//           Release build with optimization enabled
//           Fixed parameter mismatch between nr.h and frmprn.c, linmin.c (was causing 32-bit crash)
// 2.0.16b4: From https://github.com/chmunozpardo/ALMA-FETMS-beameff/tree/phase_fit_dev 20acf80,
//           converted to double-precision: nr, ScanDataRaster::calcPhaseEfficiency, calcAmplitudeEfficiency, FitPhase.cpp vars and functions.
//           But other changes not taken from that branch: sign inversions, initial linear Z search.
//           64-bit Build using msys64-mingw64 GCC 10.3.0
// 2.0.16b3: Using positive probe Z distance as nominalFocusZ
//           FitPhase uses linear search between 3D searches
//           Print pointing and invert phase options during load.
// 2.0.16b2: Improve subreflector circle visibility on FF plots.
//           Fix units (m) on NF plots.
//           Clarify code in unwrapPhase()
// 2.0.16b1: Do a linear search in Z between two multidimensional FitPhase.
// 2.0.16b0: Look for GNUPLOT_BIN in environment if not specified in input file.
//         Initialize nominal focus at -zdistance, but invert it to positive if scan phase is inverted.
//         Create output dir if not exists.
// 2.0.15: When using pointingOption::ACTUAL, use the weighted average of the pol0, pol1 beam direction in the phase fit model.
//         Moved call to calcCenterOfMass() from ScanData::AnalyzeBeams() to ScanSet::loadScan() so both pols are precomputed before the first call to FitPhase.
// 2.0.14: FitPhase: ftol=1.19e-7, reduceSubReflector=true for first pass, no z linear searches
// 2.0.13: Rotate xtics in plots.  FitPhase instrumented to make surface plots of eta_phase around the phase center.
//         Print error message if a file can't be loaded.   Fix NF axis units 'mm'
// 2.0.12: Fix bug when using reduced subreflector, but that code is disabled in this version.
//         Prints all steps of the phase center search to stdout.
//         Only makes unwrapped phase plots if "[settings] unwrapphase=1"
//         Fix bug in NF axis tics.
// 2.0.11: Unwrap FF phase defaults to false if not specified with "[settings] unwrapphase=1"
//         Added version resource for Windows.
// 2.0.10: Unwrap FF phase before phase fit; use a reduced subreflector mask for the first pass of phase fit.
//         ScanDataRaster: store phase in radians, not degrees.
//         Added plot_copol_ffphase_unwrap and plot_copol_ffphase_model to output file.
// 2.0.9:  Removed outputs shift_from_focus_mm, subreflector_shift_mm, and defocus_efficiency_due_to_moving_the_subreflector.
//         Their calculation method was broken and of dubious value.
// 2.0.8:  Display actual pointing angles on the pointing plot.
//
// 2.0.7:  Added input file option [settings] invertphase with options {lsb, usb, all, none}
//         and synonyms {yes, y, t, 1} for all, {no, n, f, 0} for none.  All options case-insensitive.
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


