/*
 * TrioMC.h
 *
 */

#ifndef PROTOCOLS_TRIOMC_H_
#define PROTOCOLS_TRIOMC_H_

#include "AstraMC.h"

template<class T>
class TrioMC : public AstraMC<T>
{
public:
    TrioMC(typename T::mac_key_type = {}, int = 0, int = 0)
    {
    }

    typename T::open_type prepare_summand(const T& secret, int my_num)
    {
        return secret[my_num - 1];
    }

    TrioMC& get_part_MC()
    {
        return *this;
    }
};

#endif /* PROTOCOLS_TRIOMC_H_ */
