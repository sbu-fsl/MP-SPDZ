/*
 * trio-party.cpp
 *
 */

#include "Protocols/TrioShare.h"
#include "Protocols/TrioMC.h"
#include "GC/AstraSecret.h"

#include "Processor/Machine.hpp"
#include "Processor/RingMachine.hpp"
#include "Protocols/AstraInput.hpp"
#include "Protocols/AstraPrep.hpp"
#include "Protocols/AstraMC.hpp"
#include "Protocols/Trio.hpp"
#include "Protocols/TrioInput.hpp"
#include "GC/ShareSecret.hpp"

int main(int argc, const char** argv)
{
    HonestMajorityRingMachine<TrioShare2, TrioShare>(argc, argv, 2);
}
