/*
 * AstraMC.h
 *
 */

#ifndef PROTOCOLS_ASTRAMC_H_
#define PROTOCOLS_ASTRAMC_H_

#include "MAC_Check_Base.h"

template<class T>
class AstraMC : public MAC_Check_Base<T>
{
public:
    AstraMC(typename T::mac_key_type = {}, int = 0, int = 0)
    {
    }

    virtual typename T::open_type prepare_summand(const T& secret, int my_num);

    void exchange(const Player& P);

    AstraMC& get_part_MC()
    {
        return *this;
    }
};

#endif /* PROTOCOLS_ASTRAMC_H_ */
