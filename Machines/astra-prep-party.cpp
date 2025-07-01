/*
 * astra-prep-party.cpp
 *
 */

#include "Protocols/AstraShare.h"
#include "GC/AstraSecret.h"

#include "Processor/Machine.hpp"
#include "Processor/RingMachine.hpp"
#include "Protocols/ReplicatedMC.hpp"
#include "Protocols/AstraInput.hpp"
#include "Protocols/AstraPrep.hpp"
#include "Protocols/Astra.hpp"
#include "Protocols/Trio.hpp"
#include "GC/ShareSecret.hpp"

int main(int argc, const char** argv)
{
    HonestMajorityRingMachine<AstraPrepShare2, AstraPrepShare>(argc, argv);
}
