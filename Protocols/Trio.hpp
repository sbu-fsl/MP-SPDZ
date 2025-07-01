/*
 * Trio.hpp
 *
 */

#ifndef PROTOCOLS_TRIO_HPP_
#define PROTOCOLS_TRIO_HPP_

#include "Trio.h"
#include "Tools/ranges.h"

#include "Astra.hpp"

template<class T>
T TrioShare<T>::local_mul_P1(const This& other) const
{
    return this->m(1) * other.lambda(1) + other.m(1) * this->lambda(1);
}

template<class T>
T TrioShare<T>::local_mul_P2(const This& other) const
{
    return this->m(2) * other.m(2);
}

template<class T>
T TrioPrepShare<T>::local_mul_P0(const This& other) const
{
    auto& x = *this;
    auto& y = other;
    return x[1] * y[1] - (x[0] - x[1]) * (y[0] - y[1]);
}

template<class T>
T TrioPrepShare<T>::local_mul_P1(const This&) const
{
    return {};
}

template<class T>
T TrioPrepShare<T>::local_mul_P2(const This&) const
{
    return {};
}

template<class T>
template<int Pi>
typename Trio<T>::pre_tuple Trio<T>::pre_mul(const T& x, const T& y)
{
    typename T::clear tmp;
#define X(I) if (Pi == I) tmp = x.local_mul_P##I(y);
    X(1) X(2)
#undef X
    return pre_dot<Pi>(tmp);
}

template<class T>
template<int Pi>
typename Trio<T>::pre_tuple Trio<T>::pre_dot(const open_type& input)
{
    auto tmp = pre_common(input);
    auto& V = tmp.second;
    auto& z = tmp.first;
    if (Pi == 1)
        os[0].store_no_resize(V + z.neg_lambda(this->astra_num));
    else
        os[0].store_no_resize(V - z.neg_lambda(this->astra_num));
    return tmp;
}

template<class T>
typename Trio<T>::pre_tuple Trio<T>::pre_common(const open_type& input)
{
    auto tmp = input + prep.get_no_check<open_type>();
    T z;
    z.neg_lambda(this->astra_num) = prep.get_no_check<open_type>();
    return {z, tmp};
}

template<class T>
size_t Trio<T>::n_mults()
{
    return this->inputs.size() + this->input_pairs.size();
}

template<class T>
void Trio<T>::exchange()
{
    CODE_LOCATION
    auto& P = this->P;
    auto& inputs = this->inputs;
    int my_num = this->astra_num;

    this->read(prep);
    if (prep.left() < open_type::size() * this->n_mults())
        throw runtime_error("not enough preprocessing");

    os[0].clear();
    os[0].reserve(open_type::size() * this->n_mults());

    results.clear();
    results.reserve(this->n_mults());

#define X(Pi) \
    if (my_num == Pi) \
    { \
        for (auto& x : this->input_pairs) \
            results.push_back(pre_mul<Pi>(x[0], x[1])); \
        for (auto& x : inputs) \
            results.push_back(pre_dot<Pi>(x)); \
    }
    X(1) X(2)

    P.pass_around(os[0], os[1], 1);

    if (os[1].left() < open_type::size() * this->n_mults())
        throw runtime_error("insufficient data in Trio");

    if (my_num == 1)
    {
        for (auto& x : results)
            x.first.m(my_num) = os[1].template get_no_check<open_type>() - x.second;
    }
    else
    {
        for (auto& x : results)
            x.first.m(my_num) = x.second - os[1].template get_no_check<open_type>();
    }

    results.reset();
}

template<class T>
T Trio<T>::finalize_mul(int)
{
    return results.next().first;
}

template<class T>
void TrioPrepProtocol<T>::pre_P0(typename T::open_type& r01,
        const typename T::open_type& input)
{
    r01.randomize(this->prng_protocol.shared_prngs[0]);
    this->os.store_no_resize(input + r01);
}

template<class T>
void TrioPrepProtocol<T>::exchange()
{
    CODE_LOCATION
    auto& P = this->P;
    auto& inputs = this->inputs;
    auto& results = this->results;
    int my_num = P.my_num();
    octetStream& os = this->os;
    os.reset_write_head();
    assert(results.size() == 0);
    auto& shared_prngs = this->prng_protocol.shared_prngs;
    typedef typename T::open_type open_type;

    this->prepare_exchange();

    results.reserve(this->n_mults);
    os.reserve(open_type::size() * 2 * this->n_mults);

    T z;
    open_type r01;

    if (my_num == 0)
    {
        for (auto& x : inputs)
        {
            pre_P0(r01, x);
            this->prng_protocol.randomize(z);
            results.push_back(z);
        }
        for (auto& x : this->input_pairs)
        {
            pre_P0(r01, x[0].local_mul_P0(x[1]));
            this->prng_protocol.randomize(z);
            results.push_back(z);
        }
        P.send_to(2, os);
    }
    else if (my_num == 1)
    {
        for (size_t i = 0; i < this->n_mults; i++)
        {
            r01.randomize(shared_prngs[1]);
            z.neg_lambda(my_num).randomize(shared_prngs[1]);
            results.push_back(z);
            os.store_no_resize(r01);
            os.store_no_resize(z.neg_lambda(my_num));
        }
        this->store(os);
    }
    else
    {
        assert(my_num == 2);
        P.receive_player(0, os);
        octetStream prep_os;
        prep_os.reserve(this->n_mults * 2 * open_type::size());
        if (os.left() < this->n_mults * open_type::size())
            throw runtime_error("insufficient data in multiplication");
        for (size_t i = 0; i < this->n_mults; i++)
        {
            prep_os.store_no_resize(os.get_no_check<open_type>());
            z.neg_lambda(my_num).randomize(shared_prngs[0]);
            prep_os.store_no_resize(z.neg_lambda(my_num));
            results.push_back(z);
        }
        this->store(prep_os);
    }

    results.reset();
}

template<class T>
template<class U>
TrioPrepShare<U> TrioPrepProtocol<T>::from_rep3(const FixedVec<U, 2>& x)
{
    TrioPrepShare<U> res;
    for (int i = 0; i < 2; i++)
        res[i] = x[this->rep_index(i)];
    if (this->my_num > 0)
        res[0] += res[1];
    return res;
}

template<class T>
T TrioPrepProtocol<T>::get_random()
{
    auto res = from_rep3(
            this->prng_protocol.template get_random<typename T::open_type>());
    this->store(res);
    return res;
}

template<class T>
void TrioPrepProtocol<T>::randoms_inst(StackedVector<T>& S,
        const Instruction& instruction)
{
    octetStream os;

    for (int j = 0; j < instruction.get_size(); j++)
    {
        auto& res = S[instruction.get_r(0) + j];
        for (int i = 0; i < 2; i++)
            res[i].randomize_part(this->prng_protocol.shared_prngs[i],
                    instruction.get_n());
        res = from_rep3(res);
        res.pack(os);
    }

    this->store(os);
}

template<class T>
template<class U, int MY_NUM>
void AstraPrepShare<T>::pre_reduced_mul(U& a, U&, U& c,
        AstraPrepProtocol<U>& protocol, const T& aa, const T&)
{
    auto& shared_prngs = protocol.prng_protocol.shared_prngs;
    auto& os = protocol.os;
    auto& os_prep = protocol.os_prep;

    switch(MY_NUM)
    {
    case 0:
        a[0].randomize(shared_prngs[0]);
        a[1] = aa - a[0];
        os.store_no_resize(a[1]);
        c = protocol.prng_protocol.template get_random<T>();
        return;
    case 1:
    {
        auto& G = shared_prngs[1];
        a[1].randomize(G);
        c[1].randomize(G);
        os_prep.store_no_resize(a[1]);
        os_prep.store_no_resize(c[1]);
        return;
    }
    case 2:
    {
        auto& G = shared_prngs[0];
        os.get_no_check(a[1]);
        c[1] = G.template get<T>();
        os_prep.store_no_resize(a[1]);
        os_prep.store_no_resize(c[1]);
        return;
    }
    }
}

template<class T>
void AstraPrepProtocol<T>::init_reduced_mul(size_t n_mul)
{
    auto& P = this->P;

    os_prep.reset_write_head();
    os.reset_write_head();
    os.reserve<typename T::clear>(n_mul);
    os_prep.reserve<typename T::clear>(2 * n_mul);

    if (P.my_num() == 2)
    {
        P.receive_player(0, os);
        os.require<typename T::clear>(n_mul);
    }
}

template<class T>
void AstraPrepProtocol<T>::exchange_reduced_mul(size_t)
{
    auto& P = this->P;

    if (P.my_num() == 0)
        P.send_to(2, os);

    if (P.my_num() == 2)
        assert(not os.left());

    this->store(os_prep);
}

template<class T>
void AstraPrepProtocol<T>::unsplit1(StackedVector<T>& dest,
        StackedVector<typename T::bit_type>& source,
        const Instruction& instruction)
{
    CODE_LOCATION
    int n_bits = instruction.get_size();
    assert(instruction.get_start().size() == 1);
    auto& P = this->P;
    init_reduced_mul(n_bits);

    if (P.my_num() == 0)
    {
        auto range = BlockRange(source, instruction.get_r(0), n_bits);
        auto dest_it = dest.iterator_for_size(instruction.get_start()[0],
                n_bits);

        for (auto& x : range)
        {
            for (auto bit : BitLeftRange(x.sum(), x, range))
            {
                T a, b, c;
                T::template pre_reduced_mul<T, 0>(a, b, c, *this, bit, {});
                *dest_it++ = a - 2 * c;
            }
        }
    }
    else if (P.my_num() == 1)
    {
        for (auto& x : Range(dest, instruction.get_start()[0], n_bits))
        {
            T a, b, c;
            T::template pre_reduced_mul<T, 1>(a, b, c, *this, {}, {});
            x[1] = a[1] - c[1] * 2;
        }
    }
    else if (P.my_num() == 2)
    {
        for (auto& x : Range(dest, instruction.get_start()[0], n_bits))
        {
            T a, b, c;
            T::template pre_reduced_mul<T, 2>(a, b, c, *this, {}, {});
            x[1] = a[1] - c[1] * 2;
        }
    }

    exchange_reduced_mul(n_bits);
}

template<class T>
template<class U, int MY_NUM>
void TrioShare<T>::pre_reduced_mul(U& a, U& b, U& c,
        Trio<U>& protocol, const T&, const T& bb)
{
    auto& os = protocol.os;
    auto& prep = protocol.prep;
    auto& results = protocol.results;

    switch(MY_NUM)
    {
    case 0:
        throw not_implemented();
    case 1:
    {
        prep.get_no_check(a[1]);
        prep.get_no_check(c[1]);
        a[0] = a[1];
        b[0] = bb;
        auto V1 = b.m(-1) * a.lambda(-1);
        os[0].store_no_resize(V1 + c.neg_lambda(-1));
        results.push_back({{}, V1});
        return;
    }
    case 2:
    {
        prep.get_no_check(a[1]);
        prep.get_no_check(c[1]);
        a[0] = a[1];
        b[0] = bb;
        auto V2 = a.m(-1) * b.m(-1);
        os[0].store_no_resize(V2 - c.neg_lambda(-1));
        results.push_back({{}, V2});
        return;
    }
    }
}

template<class T>
template<class U, int MY_NUM>
pair<U, T> TrioShare<T>::post_reduced_mul(typename U::Protocol& protocol)
{
    auto& results = protocol.results;
    auto& os = protocol.os;

    switch(MY_NUM)
    {
    case 0:
        throw not_implemented();
    case 1:
    {
        auto res = results.next();
        res.second = os[1].template get_no_check<clear>() - res.second;
        return res;
    }
    case 2:
        auto& res = results.next();
        res.second = res.second - os[1].template get_no_check<clear>();
        return res;
    }
}

template<class T>
void Trio<T>::init_reduced_mul(size_t n_mul)
{
    this->read(prep);
    prep.require<clear>(2 * n_mul);

    os[0].reset_write_head();
    os[0].template reserve<clear>(n_mul);
    results.clear();
    results.reserve(n_mul);
}

template<class T>
void Trio<T>::exchange_reduced_mul(size_t n_mul)
{
    this->P.pass_around(os[0], os[1], 1);
    os[1].template require<clear>(n_mul);
    results.reset();
    assert(results.left() == n_mul);
}

template<class T>
template<class U, int MY_NUM>
TrioShare<T> TrioShare<T>::post_input0(typename U::Protocol& protocol)
{
    This res;
    protocol.cs_prep.get_no_check(res[1]);
    return from_astra(res);
}

template<class T>
void Trio<T>::unsplit1(StackedVector<T>& dest,
        StackedVector<typename T::bit_type>& source,
        const Instruction& instruction)
{
    CODE_LOCATION
    size_t n_bits = instruction.get_size();
    assert(instruction.get_start().size() == 1);

    init_reduced_mul(n_bits);
    auto range = BlockRange(source, instruction.get_r(0), n_bits);

    if (this->my_astra_num() == 1)
    {
        for (auto& x : range)
        {
            for (typename T::clear bit : BitLeftRange(x.sum(), x, range))
            {
                T a, b, c;
                T::template pre_reduced_mul<T, 1>(a, b, c, *this, {}, bit);
                results.back().first = a + b - 2 * c;
            }
        }
    }
    else
    {
        for (auto& x : range)
        {
            for (typename T::clear bit : BitLeftRange(x.sum(), x, range))
            {
                T a, b, c;
                T::template pre_reduced_mul<T, 2>(a, b, c, *this, {}, bit);
                results.back().first = a + b - 2 * c;
            }
        }
    }

    exchange_reduced_mul(n_bits);

    if (this->my_astra_num() == 1)
    {
        for (auto& x : Range(dest, instruction.get_start()[0], n_bits))
        {
            auto res = T::template post_reduced_mul<T, 1>(*this);
            auto& mc = res.second;
            x = res.first;
            x[0] -= mc * 2;
        }
    }
    else
    {
        for (auto& x : Range(dest, instruction.get_start()[0], n_bits))
        {
            auto res = T::template post_reduced_mul<T, 2>(*this);
            auto& mc = res.second;
            x = res.first;
            x[0] -= mc * 2;
        }
    }

    assert(not os[1].left());
    assert(not prep.left());
}

#endif /* PROTOCOLS_TRIO_HPP_ */
