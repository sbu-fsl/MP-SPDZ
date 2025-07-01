/*
 * CodeLocations.h
 *
 */

#ifndef TOOLS_CODELOCATIONS_H_
#define TOOLS_CODELOCATIONS_H_

#include "Lock.h"

#include <set>
#include <tuple>
#include <string>
using namespace std;

class CodeLocations
{
    typedef tuple<string, int, string> location_type;

    static CodeLocations singleton;

    Lock lock;
    set<location_type> done;

public:
    static void maybe_output(const char* file, int line, const char* function);

    void output(const char* file, int line, const char* function);
};

#define CODE_LOCATION CodeLocations::maybe_output(__FILE__, __LINE__, __PRETTY_FUNCTION__);

#endif /* TOOLS_CODELOCATIONS_H_ */
