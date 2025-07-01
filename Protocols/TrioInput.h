/*
 * TrioInput.h
 *
 */

#ifndef PROTOCOLS_TRIOINPUT_H_
#define PROTOCOLS_TRIOINPUT_H_

#include "AstraInput.h"

template<class T>
class TrioInput : public AstraInput<T>
{
    typedef AstraInput<T> super;

public:
    TrioInput(SubProcessor<T>& proc, MAC_Check_Base<T>& MC) :
            super(proc, MC)
    {
    }
    TrioInput(SubProcessor<T>* proc, Player& P) :
            super(proc, P)
    {
    }
    TrioInput(MAC_Check_Base<T>& MC, Preprocessing<T>& prep, Player& P,
            Trio<T>* protocol) :
            super(MC, prep, P, protocol)
    {
    }

    T finalize_offset(int offset) final;
};

#endif /* PROTOCOLS_TRIOINPUT_H_ */
