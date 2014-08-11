/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file user_bitvector.hpp
    \brief user_bitvector.hpp provides a template class user_bitvector<T>
           which wraps a class T providing access to a sequence of N bits
           via operator [] and size () operators. The wrapper provides 
           additional functionality "data", "capacity", "empty", and "bit_size"
           as well as defines an iterator type t_const_uint64_iter for accessing 
           the bit-vector in 64 bit words.
    \author Shaun Harker
*/
#ifndef INCLUDED_SDSL_USER_BITVECTOR
#define INCLUDED_SDSL_USER_BITVECTOR

#include <algorithm>
#include <functional>
#include <iterator>

//! Namespace for the succinct data structure library.
namespace sdsl
{


//! An iterator class to access 64-bit words of user-defined bitvectors
/*! 
 *  With boost this could be shorted to:
 *  typedef boost::transform_iterator<std::function<uint64_t const(int64_t)>, 
 *          boost::counting_iterator<int64_t>> user_bitvector_const_iterator;
 */
class user_bitvector_const_iterator : public std::iterator<std::random_access_iterator_tag, uint64_t const, int64_t>
{
private:
    typedef typename user_bitvector<t_bitvector>::Functor Functor;
    difference_type m_pos;
    Functor m_fun;
public:
    user_bitvector_const_iterator() {}
    user_bitvector_const_iterator(difference_type pos, Functor const& fun) : m_pos(pos), m_fun(fun) {}
    user_bitvector_const_iterator(difference_type pos, Functor && fun) : m_pos(pos), m_fun(fun) {}
    value_type operator*() {return m_fun(m_pos);}
    value_type operator[](difference_type diff) {return m_fun(m_pos+diff);}
    bool operator==(const user_bitvector_const_iterator& rhs) {return m_pos==rhs.m_pos;}
    bool operator!=(const user_bitvector_const_iterator& rhs) {return m_pos!=rhs.m_pos;}
    bool operator>=(const user_bitvector_const_iterator& rhs) {return m_pos>=rhs.m_pos;}
    bool operator<=(const user_bitvector_const_iterator& rhs) {return m_pos<=rhs.m_pos;}
    bool operator>(const user_bitvector_const_iterator& rhs) {return m_pos>rhs.m_pos;}
    bool operator<(const user_bitvector_const_iterator& rhs) {return m_pos<rhs.m_pos;}
    user_bitvector_const_iterator & operator+= ( difference_type diff ) { m_pos += diff; return *this;}
    user_bitvector_const_iterator & operator-= ( difference_type diff ) { m_pos -= diff; return *this;}
    user_bitvector_const_iterator & operator++ () {++m_pos; return *this;}
    user_bitvector_const_iterator & operator-- () {--m_pos; return *this;}
    user_bitvector_const_iterator operator+ ( difference_type diff ) {return user_bitvector_const_iterator ( m_pos + diff, m_fun );}
    user_bitvector_const_iterator operator- ( difference_type diff ) {return user_bitvector_const_iterator ( m_pos - diff, m_fun );}
    user_bitvector_const_iterator operator++ (int) { user_bitvector_const_iterator tmp(*this); operator++(); return tmp; }
    user_bitvector_const_iterator operator-- (int) { user_bitvector_const_iterator tmp(*this); operator--(); return tmp; }
};

//! A class to support user-defined bitvectors
/*! 
 *  \tparam t_bitvector Bitvector class (must support operator[] and size)
 */
template<class t_bitvector>
class user_bitvector : public t_bitvector 
{
public:
    typedef uint64_t size_type;
    typedef user_bitvector_const_iterator t_const_uint64_iter;
    typedef std::function < uint64_t const ( int64_t ) > Functor;
    //typedef boost::transform_iterator<Functor, boost::counting_iterator<int64_t> > t_const_uint64_iter;

    using t_bitvector::operator[];
    using t_bitvector::size;

    user_bitvector ( t_bitvector const& bv ) : t_bitvector(bv) {}
    user_bitvector ( t_bitvector && bv ) : t_bitvector(bv) {}

    //! Returns the size of the occupied bits of the int_vector.
        /*! The capacity of a int_vector is greater or equal to the
            bit_size of the vector: capacity() >= bit_size().
            \sa size, bit_size, max_size, capacity
         */
    size_type capacity() const {
        return ((size()+63)>>6)<<6;
    }

    //! Iterator to the raw data of the user_bitvector
    /*! \returns Const iterator to the raw data of the user_bitvector
     */
    t_const_uint64_iter data() const {
        return t_const_uint64_iter ( 0, std::bind ( &user_bitvector<t_bitvector>::read_word, this, std::placeholders::_1 ));
    }

    //! Report if user_bitvector is empty or not
    bool empty () const {
        return (size()==0)? true : false;
    }

    //! Return number of bits (alternative syntax to size())
    size_type bit_size () const {
        return size ();
    }
private:
    typedef t_const_uint64_iter::difference_type difference_type;
    //! Read the 64-bit word indexed by the argument
    /*! \returns 64-bit word indexed by argument word_pos
     */
    uint64_t read_word(difference_type word_pos) const {
        return readword ( static_cast<t_bitvector const*>(this), word_pos );
    }

    // SFINAE idiom to dispatch to read_word if it is available
    template<class T>
    static auto readword_imp(T const* obj, difference_type word_pos, int) 
        -> decltype( obj -> read_word ( word_pos ) )
    {
        return obj -> read_word ( word_pos );
    }

    template<class T>
    static auto readword_imp(T const* obj, difference_type word_pos, long) -> uint64_t 
    {
        difference_type bit_pos = word_pos<<6;
        difference_type bit_pos_end = std::min(bit_pos+64,(difference_type)(obj->size()));
        uint64_t result = 0, digit = 0;
        while (bit_pos != bit_pos_end) result |= ( ((uint64_t)(*obj)[bit_pos++]) << digit++ );
        return result;
    }

    template<class T>
    static auto readword(T const* obj, difference_type word_pos) -> uint64_t
    {
        return readword_imp(obj, word_pos, 0);
    }
};

} // end namespace

#endif
