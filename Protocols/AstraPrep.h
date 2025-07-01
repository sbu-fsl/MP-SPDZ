/*
 * AstraPrep.h
 *
 */

#ifndef PROTOCOLS_ASTRAPREP_H_
#define PROTOCOLS_ASTRAPREP_H_

#include "Protocols/ReplicatedPrep.h"

template<class T> class MixedProtocolSet;

template<class T>
class AstraPrep : public ReplicatedPrep<T>
{
    SubProcessor<T>* proc;

    bool use_rep3_prep;

public:
    AstraPrep(SubProcessor<T>* proc, DataPositions& usage);
    ~AstraPrep();

    void set_protocol(typename T::Protocol& protocol);

    void buffer_bits();

    void buffer_dabits(ThreadQueues* queues);

    void buffer_edabits_with_queues(bool strict, int n_bits);
};

template<class T>
class AstraPrepPrep : public ReplicatedPrep<T>
{
    SubProcessor<T>* proc;

    typedef Rep3Share<typename T::clear> rep3_type;
    MixedProtocolSet<rep3_type>* rep3_set;

public:
    AstraPrepPrep(SubProcessor<T>* proc, DataPositions& usage);
    ~AstraPrepPrep();

    void set_protocol(typename T::Protocol& protocol);

    void buffer_bits();

    void buffer_dabits(ThreadQueues* queues);

    void buffer_edabits_with_queues(bool strict, int n_bits);
};

#endif /* PROTOCOLS_ASTRAPREP_H_ */
