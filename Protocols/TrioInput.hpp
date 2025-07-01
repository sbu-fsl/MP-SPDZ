/*
 * TrioInput.hpp
 *
 */

#ifndef PROTOCOLS_TRIOINPUT_HPP_
#define PROTOCOLS_TRIOINPUT_HPP_

#include "TrioInput.h"

template<class T>
T TrioInput<T>::finalize_offset(int offset)
{
    T res;
    if (offset == 0)
    {
        res = this->my_results.next();
        res.m(-1) = this->send_os.template get_no_check<typename T::open_type>()
                + res.neg_lambda(-1);
    }
    else
    {
        res = this->results.next();
        res.m(-1) = this->recv_os.template get_no_check<
                typename T::open_type>();
    }
    return res;
}

#endif /* PROTOCOLS_TRIOINPUT_HPP_ */
