/*
 * astra-party.cpp
 *
 */

#include "Protocols/AstraShare.h"
#include "GC/AstraSecret.h"

#include "Processor/Machine.hpp"
#include "Processor/RingMachine.hpp"
#include "Protocols/AstraInput.hpp"
#include "Protocols/AstraPrep.hpp"
#include "Protocols/AstraMC.hpp"
#include "Protocols/Astra.hpp"
#include "GC/ShareSecret.hpp"

int main(int argc, const char** argv)
{
    HonestMajorityRingMachine<AstraShare2, AstraShare>(argc, argv, 2);
}
