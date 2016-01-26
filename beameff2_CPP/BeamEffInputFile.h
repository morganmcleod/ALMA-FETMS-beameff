#ifndef BEAMEFFINPUTFILE_H_
#define BEAMEFFINPUTFILE_H_
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

#include <string>

// Forward declare implementation classes:
class scanSetList;
struct _dictionary_;
typedef _dictionary_ dictionary;

/// Handles access to the input file in INI format.
class BeamEffInputFile {
public:
    BeamEffInputFile(const std::string &inputFile);
    ///< constructor
    ///< @param inputFile filename of the input file to parse

    ~BeamEffInputFile();
    ///< destructor.

    ALMAConstants::PointingOptions getPointingOption() const
      { return pointingOption_m; }
    ///< @return the pointing option specified in the input file.

    const std::string &getOuputDirectory() const
      { return outputDir_m; }
    ///< @return the path for output files.

    const std::string &getOutputFileName() const
      { return outputFile_m; }
    ///< @return the file name specified for the output text file.

    const std::string &getDelim() const
      { return delim_m; }
    ///< @return the global delimiter from the ini file.

    const std::string &getGnuplotPath() const
      { return gnuplotPath_m; }
    ///< @return the path to gnuplot from the ini file

    const dictionary* getDict() const
      { return dict_m; }
    ///< @return const access to the ini file dictionary.

    dictionary* useDict() const
      { return dict_m; }
    ///< @return non-const access to the ini file dictionary.

    int nextScanSet(bool reset = false);
    ///< iterate over scansets defined in the file, lowest to highest.
    ///< Each scanset is identified by a unique integer.
    ///< @param reset: pass true to start at the lowest scanset.
    ///< @return integer scanset ID or 0 after the last scanset.

    const std::string &nextSection(int scanSet, bool reset = false);
    ///< iterate over the section names pertaining to a scanset.
    ///< @param scanSet: integer ID of the scanset.
    ///< @param reset: pass true to go to the first section for the given scanset.
    ///< @return section string or "" after the last section.

    bool getDatabaseKeys(std::string &section,
                         unsigned &scanSetId,
                         unsigned &FEConfigId,
                         unsigned &TestDataHeaderId);
    ///< get the database keys from a specified section
    ///< @param section: the section to read from
    ///< @param scanSetId: the value of 'scanset_id' goes here
    ///< @param FEConfigId: the value of 'fecfg' goes here
    ///< @param TestDataHeaderId: the value of 'tdh_id' goes here
    ///< @return true if the section is found and read successfully

    void print() const;
    ///< output object state.

private:
    std::string inputFile_m;            ///< name of the input file to parse.
    std::string outputDir_m;            ///< path to the output directory.
    std::string outputFile_m;           ///< name of the output text file to create.
    std::string gnuplotPath_m;          ///< path to gnuplot executable.
    std::string delim_m;                ///< delimiter specified in the input file [settings] section.
    dictionary *dict_m;                 ///< dictionary holds contents loaded from the input file.
    scanSetList *scanSets_m;            ///< implementation class for scanset and section metadata.

    ALMAConstants::PointingOptions pointingOption_m;
                                        ///< pointing option seepcified in the [settings] section.

    void loadScanSetIds();
    ///< private helper to load the scanset and section metadata
};

#endif /* BEAMEFFINPUTFILE_H_ */
