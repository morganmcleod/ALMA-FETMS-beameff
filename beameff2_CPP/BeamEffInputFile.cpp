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

#include "BeamEffInputFile.h"
#include "iniparser.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

#ifdef LINUX
#define stricmp strcasecmp
#define strnicmp strncasecmp
#endif

//----------------------------------------------------------------------------
// struct scanSetItem
//

/// Contains the scanSet ID and sections metadata.
struct scanSetItem {
    int scanSet_m;                          ///< the unique scanSet ID
    std::vector<string> sections_m;         ///< list of [scan_n] section strings in the scanSet
    std::vector<string>::const_iterator sectionIt_m;
                                            ///< iterator for the section strings

    scanSetItem(int scanset)
      : scanSet_m(scanset)
        {}
    ///< construct with a numeric ID.

    bool operator == (const scanSetItem &other) const {
        return (scanSet_m == other.scanSet_m);
    }
    ///< equality operator required for std::find()

    bool operator < (const scanSetItem &other) const {
        return (scanSet_m < other.scanSet_m);
    }
    ///< LT operator required for std::sort()
};

//----------------------------------------------------------------------------
// class scanSetList
//

/// A list of scanSetItems with some helper functions
class scanSetList : public std::vector<scanSetItem> {
public:
    scanSetList()
      {}
    ///< default construct

    int nextId(bool reset = false) {
        if (reset)
            it = begin();
        else
            ++it;
        return (it != end()) ? (it -> scanSet_m) : (0);
    }
    ///< iterate over scanSet IDs in the list
    ///< @param reset: pass true to start over at the lowest ID.
    ///< @return integer ID or 0 if past last in list.

    void appendId(int scanSet) {
        scanSetItem item(scanSet);
        push_back(item);
    }
    ///< add a new scanset ID to the list.
    ///< @param scanSet: ID to add.

    bool findId(int scanSet) const {
        const_iterator it = std::find(begin(), end(), scanSet);
        return (it != end());
    }
    ///< see if a scanset ID is already in the list.
    ///< @param scanSet: ID to find.
    ///< @return true if found

    void appendSection(int scanSet, const std::string &section) {
        iterator it = std::find(begin(), end(), scanSet);
        if (it != end())
            it -> sections_m.push_back(section);
    }
    ///< Add a [scan_n] section string for a given ID.
    ///< @param scanSet: ID to associate the section with.
    ///< @param section: string to add.

    const std::string &nextSection(int scanSet, bool reset = false) {
        iterator it = std::find(begin(), end(), scanSet);
        if (it == end())
            return emptySection;

        std::vector<string>::const_iterator &it2 = (it -> sectionIt_m);

        if (reset)
            it2 = it -> sections_m.begin();
        else
            it2++;

        if (it2 != it -> sections_m.end())
            return *it2;

        return emptySection;
    }
    ///< iterate over section strings for a given scanSet ID.
    ///< @param scanSet: ID for which to return section strings.
    ///< @param reset: pass true to start at the first section associated with the ID.
    ///< @return section string or "" if past the last section.

private:
    const_iterator it;              ///< iterator for walking the list of scanset IDs.
    const std::string emptySection; ///< empty string to return when past the last section.
};

//----------------------------------------------------------------------------
// class BeamEffInputFile
//

BeamEffInputFile::BeamEffInputFile(const std::string &inputfile)
  : inputFile_m(inputfile),
    delim_m("\t"),
    dict_m(NULL),
    scanSets_m(new scanSetList),
    pointingOption_m(ALMAConstants::NOMINAL)
{
    dict_m = iniparser_load(inputFile_m.c_str());
    if (dict_m) {
        const char *pval;

        // Read the output directory and make the output file name:
        pval = iniparser_getstring(dict_m, "settings:outputdirectory", "");
        outputDir_m = pval;
        outputFile_m = pval;
        outputFile_m += "output.txt";

        // Read the file delimiter to use when loading listing files:
        pval = iniparser_getstring(dict_m, "settings:delimiter", "tab");
        if (!stricmp(pval, "tab"))
            delim_m = "\t";
        else
            delim_m = pval;

        // Read the pointing option to use:
        pval = iniparser_getstring(dict_m, "settings:centers", "nominal");
        if (!stricmp(pval, "actual")) {
            pointingOption_m = ALMAConstants::ACTUAL;
            cout << "Using pointing: ACTUAL" << endl;
        } else if (!stricmp(pval, "7meter")) {
            pointingOption_m = ALMAConstants::ACA7METER;
            cout << "Using pointing: ACA7METER" << endl;
        } else if(!strcmp(pval, "band1test")) {
            pointingOption_m = ALMAConstants::BAND1TEST;
            cout << "Using pointing: BAND1TEST" << endl;
        } else {
            cout << "Using pointing: NOMINAL for 12-meter antennas" << endl;
        }

        // Read the invert phase and rotate FF option to use:
        pval = iniparser_getstring(dict_m, "settings:invertphase", "lsb");
        if (!stricmp(pval, "lsb")) {
            invertPhaseOption_m = ALMAConstants::INVERT_LSB;
            cout << "Inverting phase: LSB" << endl;
        } else if (!stricmp(pval, "usb")) {
            invertPhaseOption_m = ALMAConstants::INVERT_USB;
            cout << "Inverting phase: USB" << endl;
        } else if (!stricmp(pval, "all") || !stricmp(pval, "yes") || !stricmp(pval, "y") || !stricmp(pval, "t") || !stricmp(pval, "1")) {
            invertPhaseOption_m = ALMAConstants::INVERT_ALL;
            cout << "Inverting phase: ALL" << endl;
        } else if (!stricmp(pval, "none") || !stricmp(pval, "no") || !stricmp(pval, "n") || !stricmp(pval, "f") || !stricmp(pval, "0")) {
            invertPhaseOption_m = ALMAConstants::INVERT_NONE;
            cout << "Inverting phase: NONE" << endl;
        } else {
            cout << "Inverting phase: LSB" << endl;
        }

        // Read the unwrap phase option:
        int val = iniparser_getint(dict_m, "settings:unwrapphase", 0);
        unwrapPhaseOption_m = (val != 0);

        // Read the gnuplot path:
        pval = iniparser_getstring(dict_m, "settings:gnuplot", "");
        gnuplotPath_m = pval;

        // Read the gnuplot version, default 4.9:
        pval = iniparser_getstring(dict_m, "settings:gnuplot", "4.9");
        gnuplotVersion_m = pval;

        // Read it from environment variable if it was not specified:
        if (gnuplotPath_m.empty()) {
            cout << "Gnuplot not specified in input file.";
            pval = std::getenv("GNUPLOT_BIN");
            if (pval) {
                gnuplotPath_m = pval;
                cout << " Found from env: GNUPLOT_BIN=" << gnuplotPath_m;
            }
            cout << " <br>" << endl;
        }

        // Read the switch to force +90 or -90 for squint calculation
        squintOption_m = ALMAConstants::SQUINT_DEFAULT;
        pval = iniparser_getstring(dict_m, "settings:squintoption", "default");
        if (!stricmp(pval, "plus90")) {
            squintOption_m = ALMAConstants::SQUINT_PLUS90;
            cout << "Squint option: Plus90 <br>" << endl;
        } else if (!stricmp(pval, "minus90")) {
            squintOption_m = ALMAConstants::SQUINT_MINUS90;
            cout << "Squint option: Minus90 <br>" << endl;
        }
        cout << "Squint option: Default <br>" << endl;
        // Load all the scanset and section metadata:
        loadScanSetIds();
    }
}

BeamEffInputFile::~BeamEffInputFile() {
    iniparser_freedict(dict_m);
    dict_m = NULL;
    delete scanSets_m;
    scanSets_m = NULL;
}

void BeamEffInputFile::loadScanSetIds() {
    // Drop any scanset IDs we already know about:
    scanSetList &sets = *scanSets_m;
    sets.clear();

    // Loop on all sections:
    string sectionName, sectionKey;
    int numSections = iniparser_getnsec(dict_m);
    for (int i = 0; i < numSections; i++) {

        // Check that the section name starts with 'scan', e.g. not the 'settings' section:
        sectionName = iniparser_getsecname(dict_m, i);
        if (sectionName.find("scan") == 0) {

            // Look up the scanset key:
            sectionKey = sectionName;
            sectionKey += ":scanset";
            int scanSet = iniparser_getint(dict_m, sectionKey.c_str(), -1);

            // If found, see if it is already in scanSets_m:
            if (scanSet != -1) {
                if (!sets.findId(scanSet)) {
                    // If not found, append it:
                    sets.appendId(scanSet);
                }
                sets.appendSection(scanSet, sectionName);
            }
        }
        // Sort the list of scanSets:
        std::sort(sets.begin(), sets.end());
    }
}

int BeamEffInputFile::nextScanSet(bool reset) {
    return scanSets_m -> nextId(reset);
}

const std::string &BeamEffInputFile::nextSection(int scanSet, bool reset) {
    return scanSets_m -> nextSection(scanSet, reset);
}

bool BeamEffInputFile::getDatabaseKeys(std::string &section,
                                       unsigned &scanSetId,
                                       unsigned &FEConfigId,
                                       unsigned &TestDataHeaderId)
{
    scanSetId = FEConfigId = TestDataHeaderId = 0;

    if (!dict_m)
        return false;

    string sectionKey = section + ":scanset_id";
    int val = iniparser_getint(dict_m, sectionKey.c_str(), -1);
    if (val > 0)
        scanSetId = static_cast<unsigned>(val);
    else
        return false;

    sectionKey = section + ":fecfg";
    val = iniparser_getint(dict_m, sectionKey.c_str(), -1);
    if (val > 0)
        FEConfigId = static_cast<unsigned>(val);

    sectionKey = section + ":tdh_id";
    val = iniparser_getint(dict_m, sectionKey.c_str(), -1);
    if (val > 0)
        TestDataHeaderId = static_cast<unsigned>(val);

    return true;
}


void BeamEffInputFile::print() const {
    cout << "Input File: '" << inputFile_m << "'" << endl;
    cout << "Delimiter: '" << delim_m << "'" << endl;
    cout << "Pointing: '" << pointingOption_m << "'" << endl;
    cout << "Scanset Ids: '";
    scanSetList &sets = *scanSets_m;
    bool first = true;
    int id = sets.nextId(first);
    while (id > 0) {
        if (first)
            first = false;
        else
            cout << ", ";
        cout << id;
        id = sets.nextId(false);
    }
    cout << "'" << endl;
}



