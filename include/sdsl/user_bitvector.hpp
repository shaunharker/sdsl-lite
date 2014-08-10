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
           additional functionality "data" and "capacity" and defines an
           iterator type t_const_uint64_iter for accessing the bit-vector
           in 64 bit words.
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

template<class t_bitvector>
class user_bitvector_const_iterator;

//! A class to support user-defined bitvectors
/*! 
 *  \tparam t_bitvector Bitvector class (must support operator[] and size)
 */
template<class t_bitvector>
class user_bitvector : public t_bitvector 
{
public:
    typedef uint64_t size_type;
    typedef user_bitvector_const_iterator<t_bitvector> t_const_uint64_iter;
    typedef std::function < uint64_t const ( typename t_const_uint64_iter::difference_type ) > Functor;
    //typedef boost::transform_iterator<Functor, boost::counting_iterator<int64_t> > t_const_uint64_iter;

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
        return t_const_uint64_iter ( 0, std::bind ( &user_bitvector<t_bitvector>::read_word, this, _1 ));
    }
private:
    //! Read the 64-bit word indexed by the argument
    /*! \returns 64-bit word indexed by argument word_pos
     */
    uint64_t read_word(typename t_const_uint64_iter::difference_type word_pos) const {
        size_type bit_pos = word_pos>>6;
        size_type bit_pos_end = std::min(bit_pos+64,size());
        uint64_t result = 0;
        while (bit_pos != bit_pos_end) result |= (*this)[bit_pos++];
    }
};

//! An iterator class to access 64-bit words of user-defined bitvectors
/*! 
 *  \tparam t_bitvector Bitvector class (must support operator[] and size)
 */
template<class t_bitvector>
class user_bitvector_const_iterator : public std::iterator<std::random_acess_iterator_tag, uint64_t const, int64_t>
{
private:
    typedef typename user_bitvector<t_bitvector>::Functor Functor;
    difference_type m_pos;
    Functor m_fun;
public:
    using std::iterator<>::value_type;
    using std::iterator<>::difference_type;
    using std::iterator<>::pointer;
    using std::iterator<>::reference;
    using std::iterator<>::iterator_category;
    user_bitvector_const_iterator() {}
    user_bitvector_const_iterator(difference_type pos, Functor const& fun) : m_pos(pos), m_fun(fun) {}
    user_bitvector_const_iterator(difference_type pos, Functor && fun) : m_pos(pos), m_fun(fun) {}
    value_type & operator*() {return m_fun(m_pos);}
    value_type & operator[](difference_type diff) {return m_fun(m_pos+diff);}
    bool operator==(const user_bitvector_const_iterator& rhs) {return m_pos==rhs.m_pos;}
    bool operator!=(const user_bitvector_const_iterator& rhs) {return m_pos!=rhs.m_pos;}
    bool operator>=(const user_bitvector_const_iterator& rhs) {return m_pos>=rhs.m_pos;}
    bool operator<=(const user_bitvector_const_iterator& rhs) {return m_pos<=rhs.m_pos;}
    user_bitvector_const_iterator & operator+= ( difference_type diff ) { m_pos += diff; }
    user_bitvector_const_iterator & operator-= ( difference_type diff ) { m_pos -= diff; }
    user_bitvector_const_iterator & operator++ () {++m_pos; return *this;}
    user_bitvector_const_iterator & operator-- () {--m_pos; return *this;}
    user_bitvector_const_iterator operator+ ( difference_type diff ) { return user_bitvector_const_iterator(*this) += diff; }
    user_bitvector_const_iterator operator- ( difference_type diff ) { return user_bitvector_const_iterator(*this) -= diff; }
    user_bitvector_const_iterator operator++ (int) { user_bitvector_const_iterator tmp(*this); operator++(); return tmp; }
    user_bitvector_const_iterator operator-- (int) { user_bitvector_const_iterator tmp(*this); operator--(); return tmp; }
};

} // end namespace
#endif
