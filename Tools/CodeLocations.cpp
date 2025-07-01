/*
 * CodeLocations.cpp
 *
 */

#include "CodeLocations.h"
#include "Processor/OnlineOptions.h"

CodeLocations CodeLocations::singleton;

void CodeLocations::maybe_output(const char* file, int line,
        const char* function)
{
    if (OnlineOptions::singleton.code_locations)
        singleton.output(file, line, function);
}

void CodeLocations::output(const char* file, int line,
        const char* function)
{
    location_type location({file, line, function});
    lock.lock();
    if (done.find(location) == done.end())
        cerr << "first call to " << file << ":" << line << ", " << function
                << endl;
    done.insert(location);
    lock.unlock();
}
