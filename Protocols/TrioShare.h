/*
 * TrioShare.h
 *
 */

#ifndef PROTOCOLS_TRIOSHARE_H_
#define PROTOCOLS_TRIOSHARE_H_

#include "AstraShare.h"

template<class T> class Trio;
template<class T> class TrioPrepProtocol;
template<class T> class TrioMC;
template<class T> class TrioInput;

template<class T>
class TrioShare : public AstraShare<T>
{
    typedef TrioShare This;
    typedef AstraShare<T> super;

public:
    // type for clear values in relevant domain
    typedef T clear;
    typedef clear open_type;

    typedef GC::AstraSecret<TrioShare<BitVec>> bit_type;

    // opening facility
    typedef TrioMC<This> MAC_Check;
    typedef MAC_Check Direct_MC;

    // multiplication protocol
    typedef Trio<This> Protocol;

    // preprocessing facility
    typedef AstraPrep<This> LivePrep;

    // private input facility
    typedef TrioInput<This> Input;

    // default private output facility (using input tuples)
    typedef ::PrivateOutput<This> PrivateOutput;

    template<class U>
    using mc_template = TrioMC<U>;
    template<class U>
    using protocol_template = Trio<U>;
    template<class U>
    using input_template = TrioInput<U>;

    static string type_string()
    {
        return "Trio share " + T::type_string();
    }

    // used for preprocessing storage location
    static string type_short()
    {
        return "trio-" + T::type_short();
    }

    TrioShare()
    {
    }

    template<class U>
    TrioShare(const U& other) :
            super(other)
    {
    }

    T local_mul_P1(const This& other) const;
    T local_mul_P2(const This& other) const;

    T for_split(int i) const
    {
        if (i == 0)
            return common_m();
        else
            return (*this)[i];
    }

    T common_m() const
    {
        return this->m(-1) - this->neg_lambda(-1);
    }

    void set_common_m(const T& x)
    {
        (*this)[0] = x + (*this)[1];
    }

    FixedVec<T, 2> to_rep3() const
    {
        This res;
        for (int i = 0; i < 2; i++)
            res[i] = for_split(i);
        return res;
    }

    static This from_rep3(const FixedVec<T, 2>& x)
    {
        This res(x);
        res[0] += res[1];
        return res;
    }

    static This from_astra(FixedVec<T, 2>& x)
    {
        return from_rep3(x);
    }

    template<class U, int MY_NUM>
    static void pre_reduced_mul(U& a, U& b, U& c, Trio<U>& protocol,
            const T& aa, const T& bb);

    template<class U, int MY_NUM>
    static pair<U, clear> post_reduced_mul(typename U::Protocol & protocol);

    template<class U, int MY_NUM>
    static This post_input0(typename U::Protocol& protocol);

    template<class U>
    static void shrsi(SubProcessor<U>& proc, const Instruction& inst)
    {
        for (int i = 0; i < inst.get_size(); i++)
        {
            auto& dest = proc.get_S_ref(inst.get_r(0) + i);
            auto& source = proc.get_S_ref(inst.get_r(1) + i);
            dest = from_rep3(source.to_rep3() >> inst.get_n());
        }
    }
};

template<class T>
class TrioPrepShare : public AstraPrepShare<T>
{
    typedef TrioPrepShare This;
    typedef AstraPrepShare<T> super;

public:
    typedef GC::AstraSecret<TrioPrepShare<BitVec>> bit_type;

    // opening facility
    typedef ReplicatedMC<This> MAC_Check;
    typedef MAC_Check Direct_MC;

    // multiplication protocol
    typedef TrioPrepProtocol<This> Protocol;

    // preprocessing facility
    typedef AstraPrepPrep<This> LivePrep;

    // private input facility
    typedef AstraPrepInput<This> Input;

    // default private output facility (using input tuples)
    typedef ::PrivateOutput<This> PrivateOutput;

    template<class U>
    using protocol_template = TrioPrepProtocol<U>;

    // used for preprocessing storage location
    static string type_short()
    {
        return "trio-" + T::type_short();
    }

    TrioPrepShare()
    {
    }

    template<class U>
    TrioPrepShare(const U& other) :
            super(other)
    {
    }

    T local_mul_P0(const This& other) const;
    T local_mul_P1(const This& other) const;
    T local_mul_P2(const This& other) const;
};

template<int K>
using TrioShare2 = TrioShare<SignedZ2<K>>;
template<int K>
using TrioPrepShare2 = TrioPrepShare<SignedZ2<K>>;

#endif /* PROTOCOLS_TRIOSHARE_H_ */
