#ifndef STRINGCONVERT_H_
#define STRINGCONVERT_H_
/********************************************************************************
* ALMA - Atacama Large Millimeter Array
* (c) Associated Universities Inc., 2011
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
*/

/************************************************************************
 * Simple string to/from arbitrary type conversions.
 * adapted from from John Torjo's article:
 * http://articles.techrepublic.com.com/5100-10878_11-1050226.html
 * 
 *----------------------------------------------------------------------
 */
 
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>

template <typename T>
inline T from_string(const std::string &src, std::ios_base &(*f)(std::ios_base&) = std::dec) {
    std::stringstream streamIn(src);
    T ret;
    streamIn >> f >>  ret;
    return ret;
}

// specialization for string:
inline std::string from_string(const std::string &val) {
    return val;
}

template <typename T>
inline std::string to_string(const T &val, std::ios_base &(*f)(std::ios_base&) = std::dec, std::streamsize prec = -1, std::streamsize width = -1) {
    std::stringstream streamOut;
    streamOut << f;
    if (width >= 0)
        streamOut << std::setw(width);
    if (prec >= 0)
        streamOut << std::setprecision(prec);
    streamOut << val;
    return streamOut.str();
}

// specialization for string:
inline std::string to_string(const std::string &val) {
    return val;
}

// extract and return only the numeric portion of the string
extern bool stringConvertIsNotDigit(int c);

template <typename T>
inline T numericPortion(const std::string &src) {
    std::string target(src);
    std::string::iterator it = std::remove_copy_if(src.begin(), src.end(), target.begin(), stringConvertIsNotDigit);
    target.erase(it, target.end());
    return from_string<T>(target);
}

inline std::string escape_string(const std::string &src) {
    std::stringstream streamOut;
    for (std::string::const_iterator it = src.begin(); it != src.end(); ++it) {
        switch (*it) {
            case '\0':      // null
            case '\010':    // backspace
            case '\t':      // tab
            case '\n':      // newline
            case '\r':      // carriage return
            case '\"':      // double quote
            case '\'':      // single quote
            case '\\':      // backslash
            case '\140':    // grave accent
                streamOut << '\\';
                // and drop through

            default:
                streamOut << *it;
        }
    }
    return streamOut.str();
}

inline std::string unescape_string(const std::string &src) {
    std::stringstream streamOut;
    std::string::const_iterator it = src.begin();
    while (it != src.end()) {
        if (*it == '\\')
            ++it;
        streamOut << *it;
    }
    return streamOut.str();
}

#endif /*STRINGCONVERT_H_*/

