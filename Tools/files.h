/*
 * files.h
 *
 */

#ifndef TOOLS_FILES_H_
#define TOOLS_FILES_H_

#include <iostream>

void open_with_check(std::ifstream& stream, const std::string& filename)
{
    stream.open(filename);
    if (not stream.good())
        throw runtime_error("cannot open " + filename);
}

#endif /* TOOLS_FILES_H_ */
