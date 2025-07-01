/*
 * DoubleRange.h
 *
 */

#ifndef TOOLS_DOUBLERANGE_H_
#define TOOLS_DOUBLERANGE_H_

template<class T>
class DoubleIterator
{
    typedef typename vector<T>::iterator iterator;

public:
    iterator left, right;

    DoubleIterator(iterator left, iterator right) :
            left(left), right(right)
    {
    }

    bool operator!=(const DoubleIterator& other) const
    {
        return left != other.left;
    }

    DoubleIterator& operator++()
    {
        left++;
        right++;
        return *this;
    }

    pair<T&, T&> operator*()
    {
        return {*left, *right};
    }
};

template<class T>
class DoubleRange
{
    StackedVector<T>& S;
    size_t left, right, size;

public:
    DoubleRange(StackedVector<T>& S, size_t left, size_t right, size_t size) :
            S(S), left(left), right(right), size(size)
    {
    }

    DoubleIterator<T> begin()
    {
        return {S.iterator_for_size(left, size),
            S.iterator_for_size(right, size)};
    }

    DoubleIterator<T> end()
    {
        return {S.iterator_for_size(left + size, 0),
            S.iterator_for_size(right + size, 0)};
    }
};

#endif /* TOOLS_DOUBLERANGE_H_ */
