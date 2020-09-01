ALMA FETMS Beam Efficiency Calculator
beam_eff2.exe version 2.0.12

// Version history
// 2.0.12: Fix bug when using reduced subreflector, but that code is disabled in this version.
//         Prints all steps of the phase center search to stdout.
//         Only makes unwrapped phase plots if "[settings] unwrapphase=1"
//         Fix bug in NF axis tics.
// 2.0.11: Unwrap FF phase defaults to false if not specified with "[settings] unwrapphase=1"
//         Added version resource for Windows.
// 2.0.10: Unwrap FF phase before phase fit; use a reduced subreflector mask for the first pass of phase fit.
//         ScanDataRaster: store phase in radians, not degrees.
//         Added plot_copol_ffphase_unwrap and plot_copol_ffphase_model to output file.