/*
 * AstraInput.cpp
 *
 */

#include "AstraInput.h"

template<class T>
AstraInput<T>::AstraInput(SubProcessor<T>& proc, MAC_Check_Base<T>&) :
        InputBase<T>(&proc), protocol(proc.protocol)
{
}

template<class T>
AstraInput<T>::AstraInput(SubProcessor<T>* proc, Player&) :
        AstraInput<T>(*proc, proc->MC)
{
}

template<class T>
AstraInput<T>::AstraInput(MAC_Check_Base<T>&, Preprocessing<T>&,
        Player&, AstraOnlineBase<T>* protocol) :
        protocol(*protocol)
{
}

template<class T>
AstraPrepInput<T>::AstraPrepInput(SubProcessor<T>& proc, MAC_Check_Base<T>&) :
        InputBase<T>(&proc), protocol(proc.protocol)
{
    n_inputs.resize(3);
}

template<class T>
AstraPrepInput<T>::AstraPrepInput(SubProcessor<T>* proc, Player&) :
        AstraPrepInput<T>(*proc, proc->MC)
{
}

template<class T>
AstraPrepInput<T>::AstraPrepInput(MAC_Check_Base<T>&, Preprocessing<T>&,
        Player&, AstraPrepProtocol<T>* protocol) :
        protocol(*protocol)
{
    n_inputs.resize(3);
}

template<class T>
bool AstraPrepInput<T>::is_me(int player, int)
{
    assert(this->P);
    return player + 1 == this->P->my_num();
}

template<class T>
void AstraInput<T>::reset(int player)
{
    if (protocol.P.my_num() == player)
    {
        send_os.reset_write_head();
        inputs.clear();
        my_results.clear();
    }
    results.clear();
}

template<class T>
void AstraPrepInput<T>::reset(int player)
{
    if (is_me(player))
        prep_os.reset_write_head();
    n_inputs.at(player) = 0;
}

template<class T>
void AstraPrepInput<T>::add_mine(const typename T::open_type&,
        int)
{
    add_other(this->P->my_num() - 1);
}

template<class T>
void AstraPrepInput<T>::add_other(int player, int)
{
    assert(player != 2);
    n_inputs.at(player)++;
}

template<class T>
void AstraInput<T>::add_mine(const typename T::open_type& input,
        int)
{
    inputs.push_back(input);
}

template<class T>
void AstraInput<T>::add_other(int, int)
{
    results.push_back({});
}

template<class T>
void AstraPrepInput<T>::exchange()
{
    CODE_LOCATION
    if (OnlineOptions::singleton.has_option("verbose_astra"))
    {
        for (int i = 0; i < 3; i++)
            fprintf(stderr, "astra input from %d exchange %d\n", i, n_inputs.at(i));
    }

    auto& P = *this->P;

    for (auto& x : results)
        x.clear();

    if (P.my_num() == 0)
        for (int i = 0; i < 2; i++)
        {
            T res;
            int n = n_inputs.at(i);
            results[1 + i].reserve(n);
            for (int j = 0; j < n; j++)
            {
                res[i].randomize(protocol.prng_protocol.shared_prngs[i]);
                results[1 + i].push_back(res);
            }
        }
    else
    {
        // other online player
        int other = P.my_num() == 1;
        int offset = 1 + (P.my_num() == 2);
        results[offset].resize(n_inputs.at(other));

        // my inputs
        T res;
        int my_astra_num = P.my_num() - 1;
        int n = n_inputs.at(my_astra_num);
        prep_os.reserve(n * open_type::size());
        results[0].reserve(n);
        for (int j = 0; j < n; j++)
        {
            auto gamma =
                    protocol.prng_protocol.shared_prngs[other].template get<
                            typename T::open_type>();
            prep_os.store_no_resize(gamma);
            res.neg_lambda(-1) = gamma;
            results[0].push_back(res);
        }
    }

    protocol.store(prep_os);

    for(auto& x : results)
        x.reset();
}

template<class T>
void AstraInput<T>::exchange()
{
    CODE_LOCATION
    if (OnlineOptions::singleton.has_option("verbose_astra"))
        fprintf(stderr, "astra input exchange %zu\n",
                this->inputs.size() * T::open_type::size());

    octetStream prep_os;
    protocol.read(prep_os);
    my_results.reserve(inputs.size());
    send_os.reserve(inputs.size() * T::open_type::size());

    if (prep_os.left() < inputs.size() * open_type::size())
        throw runtime_error("insufficient data in input");

    for (auto& input : inputs)
    {
        auto gamma = prep_os.get_no_check<typename T::open_type>();
        send_os.store_no_resize(input - gamma);
        T res;
        res.neg_lambda(-1) = gamma;
        my_results.push_back(res);
    }

    assert(send_os.left() == my_results.size() * T::open_type::size());

    protocol.P.pass_around(send_os, recv_os, 1);

    if (recv_os.left() < results.size() * T::open_type::size())
        throw runtime_error("insufficient data in Astra input");

    assert(not prep_os.left());

    results.reset();
    my_results.reset();
}

template<class T>
T AstraPrepInput<T>::finalize_offset(int offset)
{
    return results[offset].next();
}

template<class T>
T AstraPrepInput<T>::finalize(int player, int)
{
    return finalize_offset((player - (this->P->my_num() - 1) + 3) % 3);
}

template<class T>
T AstraInput<T>::finalize(int player, int)
{
    return finalize_offset(player - this->P->my_num());
}

template<class T>
T AstraInput<T>::finalize_offset(int offset)
{
    T res;
    octetStream* o;
    if (offset == 0)
    {
        res = my_results.next();
        o = &send_os;
    }
    else
    {
        res = results.next();
        o = &recv_os;
    }
    res.m(-1) = o->get_no_check<typename T::open_type>();
    return res;
}
