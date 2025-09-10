/*
 * lowgear-offline.cpp
 *
 */

#include "SPDZ.hpp"
#include "Math/gfp.hpp"
#include "Protocols/LowGearShare.h"
#include "Protocols/CowGearOptions.h"
#include "Protocols/CowGearPrep.hpp"
#include "Processor/FieldMachine.hpp"
#include "Processor/OfflineMachine.hpp"

int main(int argc, const char** argv)
{
    ez::ezOptionParser opt;
    CowGearOptions::singleton = CowGearOptions(opt, argc, argv, false);
    DishonestMajorityFieldMachine<LowGearShare, LowGearShare, gf2n_short,
            OfflineMachine<DishonestMajorityMachine>>(argc, argv, opt);
}
