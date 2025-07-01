/*
 * AstraPrep.cpp
 *
 */

#include "AstraPrep.h"
#include "Rep3Share.h"
#include "ProtocolSet.h"
#include "GC/SemiHonestRepPrep.h"

#include "Rep3Shuffler.hpp"
#include "ReplicatedMC.hpp"
#include "GC/RepPrep.hpp"

template<class T>
AstraPrep<T>::AstraPrep(SubProcessor<T>* proc, DataPositions& usage) :
        BufferPrep<T>(usage), BitPrep<T>(proc, usage),
        ReplicatedRingPrep<T>(proc, usage), RingPrep<T>(proc, usage),
        SemiHonestRingPrep<T>(proc, usage), ReplicatedPrep<T>(proc, usage),
        proc(proc)
{
    use_rep3_prep = OnlineOptions::singleton.has_option("rep3_prep");
}

template<class T>
AstraPrepPrep<T>::AstraPrepPrep(SubProcessor<T>* proc, DataPositions& usage) :
        BufferPrep<T>(usage), BitPrep<T>(proc, usage),
        ReplicatedRingPrep<T>(proc, usage), RingPrep<T>(proc, usage),
        SemiHonestRingPrep<T>(proc, usage), ReplicatedPrep<T>(proc, usage),
        proc(proc), rep3_set(0)
{
}

template<class T>
AstraPrep<T>::~AstraPrep()
{
    this->protocol = 0;
}

template<class T>
AstraPrepPrep<T>::~AstraPrepPrep()
{
    this->protocol = 0;
}

template<class T>
void AstraPrep<T>::set_protocol(typename T::Protocol& protocol)
{
    this->protocol = &protocol;
}

template<class T>
void AstraPrepPrep<T>::set_protocol(typename T::Protocol& protocol)
{
    this->protocol = &protocol;
    if (not T::clear::characteristic_two
            and OnlineOptions::singleton.has_option("rep3_prep"))
    {
        assert(not rep3_set);
        rep3_set = new MixedProtocolSet<rep3_type>(protocol.P, {}, {});
    }
}

template<class T>
void AstraPrep<T>::buffer_bits()
{
    this->base_player = 0;
    this->buffer_bits_without_check();
}

template<class T>
void AstraPrepPrep<T>::buffer_bits()
{
    this->base_player = 0;
    this->buffer_bits_without_check();
}

template<class T>
void AstraPrep<T>::buffer_dabits(ThreadQueues* queues)
{
    if (use_rep3_prep)
    {
        assert(this->protocol);
        octetStream os;
        this->protocol->read(os);
        os.get(this->dabits);
    }
    else
    {
        this->base_player = 0;
        SemiHonestRingPrep<T>::buffer_dabits(queues);
    }
}

template<class T>
void AstraPrepPrep<T>::buffer_dabits(ThreadQueues* queues)
{
    if (rep3_set)
    {
        for (int i = 0; i < T::bit_type::default_length; i++)
        {
            rep3_type a;
            typename rep3_type::bit_type b;
            rep3_set->preprocessing.get_dabit(a, b);
            this->dabits.push_back({this->protocol->from_rep3(a),
                this->protocol->from_rep3(b)});
        }
        assert(this->protocol);
        octetStream os;
        os.store(this->dabits);
        this->protocol->store(os);
    }
    else
    {
        this->base_player = 0;
        SemiHonestRingPrep<T>::buffer_dabits(queues);
    }
}

template<class T>
void AstraPrep<T>::buffer_edabits_with_queues(bool strict, int n_bits)
{
    if (use_rep3_prep)
    {
        assert(this->protocol);
        octetStream os;
        this->protocol->read(os);
        os.get(this->edabits[{strict, n_bits}]);
    }
    else
        SemiHonestRingPrep<T>::buffer_edabits_with_queues(strict, n_bits);
}

template<class T>
void AstraPrepPrep<T>::buffer_edabits_with_queues(bool strict,
        int n_bits)
{
    if (rep3_set)
    {
        assert(this->protocol);
        auto rep3 = rep3_set->preprocessing.get_edabitvec(strict, n_bits);
        this->edabits[{strict, n_bits}].push_back({});
        auto& tmp = this->edabits[{strict, n_bits}].back();
        for (auto& x : rep3.a)
            tmp.a.push_back(this->protocol->from_rep3(x));
        for (auto& x : rep3.b)
            tmp.b.push_back(this->protocol->from_rep3(x));
        octetStream os;
        os.store(this->edabits[{strict, n_bits}]);
        this->protocol->store(os);
    }
    else
        SemiHonestRingPrep<T>::buffer_edabits_with_queues(strict, n_bits);
}
