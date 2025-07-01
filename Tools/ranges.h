/*
 * ranges.h
 *
 */

#ifndef TOOLS_RANGES_H_
#define TOOLS_RANGES_H_

template<class T>
class BlockRange
{
    typename vector<T>::iterator begin_, end_;

public:
    const int n_bits;

    BlockRange(StackedVector<T>& whole, size_t start, size_t n_bits) : n_bits(n_bits)
    {
        begin_ = whole.begin() + start;
        end_ = begin_ + DIV_CEIL(n_bits, T::default_length);
        assert(end_ <= whole.end());
    }

    auto begin()
    {
        return begin_;
    }

    auto end()
    {
        return end_;
    }
};

template<class T>
class BitLeftIterator
{
    int i;
    const T& element;

public:
    BitLeftIterator(int i, const T& element) :
            i(i), element(element)
    {
    }

    bool operator!=(BitLeftIterator& other) const
    {
        return i != other.i;
    }

    void operator++()
    {
        i++;
    }

    T operator*() const
    {
        return element.get_bit(i);
    }
};

template<class T, class U>
class BitLeftRange
{
    int end_;
    T element;

public:
    BitLeftRange(const T& element, const U& current, BlockRange<U>& range) :
            element(element)
    {
        if (&current < &*range.end() - 1)
            end_ = U::default_length;
        else
            end_ = range.n_bits % U::default_length;

        if (not end_)
            end_ = U::default_length;
    }

    BitLeftIterator<T> begin()
    {
        return {0, element};
    }

    BitLeftIterator<T> end()
    {
        return {end_, element};
    }
};

#endif /* TOOLS_RANGES_H_ */
