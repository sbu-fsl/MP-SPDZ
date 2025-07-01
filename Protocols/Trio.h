/*
 * Trio.h
 *
 */

#ifndef PROTOCOLS_TRIO_H_
#define PROTOCOLS_TRIO_H_

#include "Astra.h"
#include "TrioShare.h"

template<class T>
class Trio : public AstraOnlineBase<T>
{
    friend class TrioShare<typename T::clear>;

    typedef AstraOnlineBase<T> super;

    typedef typename T::open_type open_type;
    typedef typename T::clear clear;
    typedef pair<T, typename T::open_type> pre_tuple;

    IteratorVector<pre_tuple> results;
    octetStream prep, os[2];

    size_t n_mults();

    template<int>
    pre_tuple pre_mul(const T& x, const T& y);
    template<int>
    pre_tuple pre_dot(const open_type& input);
    pre_tuple pre_common(const open_type& input);

    void init_reduced_mul(size_t n_mul);
    void exchange_reduced_mul(size_t n_mul);

public:
    Trio(Player& P) :
            AstraOnlineBase<T>(P)
    {
    }

    void exchange();

    T finalize_mul(int = -1) final;
    T finalize_mul_fast()
    {
        return finalize_mul();
    }

    void unsplit1(StackedVector<T>& dest,
            StackedVector<typename T::bit_type>& source,
            const Instruction& instruction);
};

template<class T>
class TrioPrepProtocol : public AstraPrepProtocol<T>
{
    typedef AstraPrepProtocol<T> super;

    typedef typename T::open_type open_type;

    octetStream os, prep_os;

    void pre_P0(open_type& r01, const open_type& input);

public:
    TrioPrepProtocol(Player& P) :
            AstraPrepProtocol<T>(P)
    {
    }

    bool local_mul_for(int player)
    {
            return player == this->my_astra_num() and player == 0;
    }

    void exchange();

    template<class U>
    TrioPrepShare<U> from_rep3(const FixedVec<U, 2>& x);

    T get_random();
    void randoms_inst(StackedVector<T>&, const Instruction&);
};

#endif /* PROTOCOLS_TRIO_H_ */
