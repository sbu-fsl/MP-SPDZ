/*
 * trio-prep-party.cpp
 *
 */

#include "Protocols/TrioShare.h"
#include "GC/AstraSecret.h"

#include "Processor/Machine.hpp"
#include "Processor/RingMachine.hpp"
#include "Protocols/AstraInput.hpp"
#include "Protocols/AstraPrep.hpp"
#include "Protocols/AstraMC.hpp"
#include "Protocols/Trio.hpp"
#include "GC/ShareSecret.hpp"

int main(int argc, const char** argv)
{
    HonestMajorityRingMachine<TrioPrepShare2, TrioPrepShare>(argc, argv);
}
