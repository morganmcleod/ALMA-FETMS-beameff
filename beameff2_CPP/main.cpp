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
#include "BeamEffInputFile.h"
#include "ScanSet.h"
#include "ScanData.h"

#include <iostream>
using namespace std;

int DEBUGGING_NR = 0;

int main(int argc, char *argv[]) {
    cout << "<br> ********************************************" << endl;
    cout << "<br> Beam Efficiency Calculator Version  " << BEAMEFF_SW_VERSION_STRING << endl;
    cout << "<br> ********************************************\n<br>" << endl;

    // temporary input file name can be put here for debugging:
    string inputFile = "T:\\fecfg866\\ssid1647\\input_file.txt";
                     //"T:\\20151210\\output2\\12_10_2015__9_16 AM_INPUT.txt";
                     //"T:\\ssid1647\\MMoutput\\12_9_2015__5_33 PM_INPUT.txt";
    if (argc > 1)
        inputFile = argv[1];
    if (inputFile.empty()) {
        cout << "Must provide input file as command line parameter.  Stopping." << endl;
        return -1;
    }

    // load the input file:
    BeamEffInputFile inFile(inputFile);
    inFile.print();
    cout << "----" << endl;

    // loop on scanset IDs in the file:
    string section;
    int scanSetId = inFile.nextScanSet(true);
    while (scanSetId) {
        // create a ScanSet object to hold the scans:
        ScanSet SS(scanSetId);

        // loop on sections belonging to the scanset ID.
        section = inFile.nextSection(scanSetId, true);
        while (!section.empty()) {
            // load the section into a ScanData object in the ScanSet:
            SS.loadScan(inFile.getDict(), section, inFile.getDelim());
            // go to next section for the scanSet:
            section = inFile.nextSection(scanSetId);
        }
        SS.getEfficiencies(inFile.getPointingOption());

        SS.print();
        cout << "----" << endl;

        scanSetId = inFile.nextScanSet();
    }

    return 0;
}


