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
    string inputFile = "";
                     //"T:\\validation2.0\\fecfg1076\\ssid2091\\output2\\1_13_2016__4_41 PM_INPUT.txt";
                     //"T:\\validation2.0\\fecfg1076\\ssid1749\\output2\\1_13_2016__2_50 PM_INPUT.txt";
                     //"T:\\validation2.0\\fecfg866\\ssid1647\\output2\\1_13_2016__4_52 PM_INPUT.txt";
                     //"T:\\validation2.0\\RF50\\output2\\1_12_2016__3_59 PM_INPUT.txt";

    if (argc > 1)
        inputFile = argv[1];
    if (inputFile.empty()) {
        cout << "Must provide input file as command line parameter.  Stopping." << endl;
        return -1;
    }

    cout << "input file: " << inputFile << "<br>" << endl;

    // load the input file:
    BeamEffInputFile inFile(inputFile);
//    inFile.print();
//    cout << "----" << endl;

    string outputFile = inFile.getOutputFileName();
    cout << "output file: " << outputFile << "<br>" << endl;

    // loop on scanset IDs in the file:
    string section;

    unsigned scanSet = inFile.nextScanSet(true);
    while (scanSet) {
        cout << "scanset: " << scanSet << "<br>" << endl;

        // create a ScanSet object to hold the scans:
        ScanSet SS(scanSet);

        // initialize database keys to zero for each scanset:
        unsigned scanSetId(0);
        unsigned FEConfigId(0);
        unsigned TestDataHeaderId(0);

        // loop on sections belonging to the scanset ID.
        section = inFile.nextSection(scanSet, true);
        while (!section.empty()) {
            cout << "section: " << section << "<br>" << endl;

            // If we haven't assigned the database keys yet,
            if (scanSetId == 0) {
                // Load them from the current section and store into the ScanSet:
                if (inFile.getDatabaseKeys(section, scanSetId, FEConfigId, TestDataHeaderId))
                    SS.setDatabaseKeys(scanSetId, FEConfigId, TestDataHeaderId);
            }

            // load the section into a ScanData object in the ScanSet:
            SS.loadScan(inFile.getDict(), section, inFile.getDelim(), inFile.getInvertPhaseOption());
            // go to next section for the scanSet:
            section = inFile.nextSection(scanSet);
        }
        // calculate beam efficiencies:
        SS.calcEfficiencies(inFile.getPointingOption());

        // Generate all the plots:
        SS.makePlots(inFile.getOuputDirectory(), inFile.getGnuplotPath());

//        SS.print();
//        cout << "----" << endl;

        // Write the summary output file containing efficiencies and computed results:
        SS.writeOutputFile(inFile.useDict(), outputFile);

        // Go to the next scanset in the input file:
        scanSet = inFile.nextScanSet();
    }

    cout << "done.<br>" << endl;
    return 0;
}


