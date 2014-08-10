/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

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
/*! \file select_support.hpp
    \brief select_support.hpp contains classes that support a sdsl::bit_vector with constant time select information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT
#define INCLUDED_SDSL_SELECT_SUPPORT

/** \defgroup select_support_group Select Support (SCS)
 * This group contains data structures which support an sdsl::bit_vector with the select method.
 */

#include "int_vector.hpp"
#include "rank_support.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{
//! The base class of classes supporting select queries for a sdsl::bit_vector in constant time.
/*! Abstract base class for classes supporting select queries.
 *  \tparam t_bitvector Bitvector class (must support operator[], size, capacity, and data)
 */
template<class t_bitvector=int_vector<1>>
class select_support
{
    protected:
        const t_bitvector* m_v; //!< Pointer to the select supported sdsl::bit_vector.
    public:
        typedef typename t_bitvector::size_type size_type;
        const t_bitvector* v;

        //! Constructor of select_support.
        /*! \param v The bit_vector to support rank queries.
         */
        select_support(const t_bitvector* f_v=nullptr):v(f_v) {
            m_v = f_v;
        }
        //! Copy constructor
        /*! Copy the whole select_support including the  pointer
         *  to the supported bit_vector.
         */
        select_support(const select_support& f_v);
        //! Destructor of select_support.
        virtual ~select_support() {};

        //! Select returns the index of the i-th 1-bit in the supported bit_vector.
        /*!	\param i Argument to calculate the index of the i-th 1-bit in the supported bit_vector.
        	\return The index \f$\in [0..v.size()-1]\f$ of the i-th 1-bit in the supported bit_vector.
        	Call init or load to initialize the data structure before the first call of this method.
         	\sa init, load.
         */
        virtual size_type select(size_type i) const = 0;

        //! Alias for select
        virtual size_type operator()(size_type i) const = 0;
        //! Serialize the select_support to an out file stream.
        virtual size_type serialize(std::ostream& out, structure_tree_node* v, std::string name)const = 0;
        //! Load the select_support from an in file stream.
        /*!	Load an previously serialized select_support from a std::istream.
        	\param in The std::istream to load the select_support.
        	\param v The bitvector to be supported.
        	\sa init, select.
         */
        virtual void load(std::istream& in, const t_bitvector* v=nullptr) = 0;

        //! This method sets the supported bit_vector
        virtual void set_vector(const t_bitvector* v=nullptr) = 0;
};


template<uint8_t bit_pattern, uint8_t pattern_len, class t_bitvector>
struct select_support_trait {
    typedef typename select_support<t_bitvector>::size_type size_type;
    typedef decltype(std::declval<t_bitvector const>().data()) t_const_uint64_iter; 
    /* Count the number of arguments for the specific select support */
    static size_type arg_cnt(const t_bitvector&) {
        return 0;
    }

    static size_type args_in_the_first_word(uint64_t, uint8_t, uint64_t) {
        return 0;
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t, size_type, uint8_t, uint64_t) {
        return 0;
    }

    static size_type args_in_the_word(uint64_t, uint64_t&) {
        return 0;
    }

    static size_type ith_arg_pos_in_the_word(uint64_t, size_type, uint64_t) {
        return 0;
    }

    static bool found_arg(size_type, const t_bitvector&) {
        return 0;
    }

    static uint64_t init_carry(t_const_uint64_iter, size_type) {
        return 0;
    }

    static uint64_t get_carry(uint64_t) {
        return 0;
    }
};

template<class t_bitvector>
struct select_support_trait<0,1,t_bitvector> {
    typedef typename select_support<t_bitvector>::size_type size_type;
    typedef decltype(std::declval<t_bitvector const>().data()) t_const_uint64_iter; 

    static size_type arg_cnt(const t_bitvector& v) {
        return v.bit_size()-util::cnt_one_bits(v);
    }
    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t) {
        return bits::cnt((~w)& bits::lo_unset[offset]);
    }
    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t) {
        return bits::sel(~w & bits::lo_unset[offset], i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t&) {
        return bits::cnt(~w);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t) {
        return bits::sel(~w, i);
    }
    static bool found_arg(size_type i, const t_bitvector& v) {
        return !v[i];
    }
    static uint64_t init_carry(t_const_uint64_iter, size_type) {
        return 0;
    }
    static uint64_t get_carry(uint64_t) {
        return 0;
    }
};

template<class t_bitvector>
struct select_support_trait<1,1,t_bitvector> {
    typedef typename select_support<t_bitvector>::size_type	size_type;
    typedef decltype(std::declval<t_bitvector const>().data()) t_const_uint64_iter; 

    static size_type arg_cnt(const t_bitvector& v) {
        return util::cnt_one_bits(v);
    }
    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t) {
        return bits::cnt(w & bits::lo_unset[offset]);
    }
    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t) {
        return bits::sel(w & bits::lo_unset[offset], i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t&) {
        return bits::cnt(w);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t) {
        return bits::sel(w, i);
    }
    static bool found_arg(size_type i, const t_bitvector& v) {
        return v[i];
    }
    static uint64_t init_carry(t_const_uint64_iter, size_type) {
        return 0;
    }
    static uint64_t get_carry(uint64_t) {
        return 0;
    }
};

template<class t_bitvector>
struct select_support_trait<10,2,t_bitvector> {
    typedef typename select_support<t_bitvector>::size_type size_type;
    typedef decltype(std::declval<t_bitvector const>().data()) t_const_uint64_iter; 

    static size_type arg_cnt(const t_bitvector& v) {
        return util::cnt_onezero_bits(v);
    }
    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
        return bits::cnt(bits::map10(w, carry) & bits::lo_unset[offset]);
    }
    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
        return bits::sel(bits::map10(w, carry) & bits::lo_unset[offset], i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return bits::cnt10(w, carry);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
        return bits::sel(bits::map10(w, carry), i);
    }
    static bool found_arg(size_type i, const t_bitvector& v) {
        if (i > 0 and v[i-1] and !v[i])
            return true;
        return false;
    }
    static uint64_t init_carry(t_const_uint64_iter data, size_type word_pos) {
        return word_pos ? (*(data-1)>>63) : 0;
    }
    static uint64_t get_carry(uint64_t w) {
        return w>>63;
    }
};

template<class t_bitvector>
struct select_support_trait<01,2,t_bitvector> {
    typedef typename select_support<t_bitvector>::size_type size_type;
    typedef decltype(std::declval<t_bitvector const>().data()) t_const_uint64_iter; 

    static size_type arg_cnt(const t_bitvector& v) {
        return util::cnt_zeroone_bits(v);
    }
    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
        return bits::cnt(bits::map01(w, carry) & bits::lo_unset[offset]);
    }
    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
        return bits::sel(bits::map01(w, carry) & bits::lo_unset[offset], i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return bits::cnt01(w, carry);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
        return bits::sel(bits::map01(w, carry), i);
    }
    static bool found_arg(size_type i, const t_bitvector& v) {
        if (i > 0 and !v[i-1] and v[i])
            return true;
        return false;
    }
    static uint64_t init_carry(t_const_uint64_iter data, size_type word_pos) {
        return word_pos ? (*(data-1)>>63) : 1;
    }
    static uint64_t get_carry(uint64_t w) {
        return w>>63;
    }
};

} // end namespace sdsl

#include "select_support_mcl.hpp"
#include "select_support_scan.hpp"

#endif
