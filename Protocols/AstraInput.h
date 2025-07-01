/*
 * AstraInput.h
 *
 */

#ifndef PROTOCOLS_ASTRAINPUT_H_
#define PROTOCOLS_ASTRAINPUT_H_

#include "Processor/Input.h"

template<class T> class AstraOnlineBase;

template<class T>
class AstraInput : public InputBase<T>
{
protected:
    typedef typename T::open_type open_type;

    octetStream send_os, recv_os;
    AstraOnlineBase<T>& protocol;
    IteratorVector<T> results, my_results;
    CheckVector<typename T::open_type> inputs;

public:
    AstraInput(SubProcessor<T>& proc, MAC_Check_Base<T>& MC);
    AstraInput(SubProcessor<T>* proc, Player& P);
    AstraInput(MAC_Check_Base<T>& MC, Preprocessing<T>&, Player& P,
            AstraOnlineBase<T>* protocol);

    void reset(int player);
    void add_mine(const typename T::open_type& input, int n_bits = -1) final;
    void add_other(int player, int n_bits = -1) final;

    void exchange();
    T finalize(int player, int n_bits = -1) final;
    virtual T finalize_offset(int offset);
};

template<class T>
class AstraPrepInput : public InputBase<T>
{
    typedef typename T::open_type open_type;

    AstraPrepProtocol<T>& protocol;
    array<IteratorVector<T>, 3> results;
    octetStream prep_os;
    vector<int> n_inputs;

public:
    AstraPrepInput(SubProcessor<T>& proc, MAC_Check_Base<T>& MC);
    AstraPrepInput(SubProcessor<T>* proc, Player& P);
    AstraPrepInput(MAC_Check_Base<T>& MC, Preprocessing<T>&, Player& P,
            AstraPrepProtocol<T>* protocol);

    bool is_me(int player, int = -1);

    void reset(int player);
    void add_mine(const typename T::open_type& input, int n_bits = -1) final;
    void add_other(int player, int n_bits = -1) final;

    void exchange();
    T finalize(int player, int n_bits = -1) final;
    T finalize_offset(int player);
};

#endif /* PROTOCOLS_ASTRAINPUT_H_ */
