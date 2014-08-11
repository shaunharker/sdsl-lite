#include <iostream>
#include <utility>
#include <algorithm>
#include <chrono>
#include <stdexcept>

#include "sdsl/util.hpp"
#include "sdsl/user_bitvector.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/rank_support_v5.hpp"

/// slow_bitvector   (basic example -- "slow" since it doesn't have optimized access)
///   We want to make a bitvector z out of bitvectors x and y
///   such that: z[i] := x[i] || y[rank0(x,i)];
///   We provide methods
///      * size
///      * operator []
///   but not the "read_word" method (see fast_bitvector below)
class slow_bitvector {
public:
  slow_bitvector ( sdsl::bit_vector const* x, sdsl::bit_vector const* y ) : m_x(x),m_y(y) {
    m_rankx0 = sdsl::rank_support_v5<0>(m_x); 
    m_size = m_x -> size ();   
  }
  uint64_t size ( void ) const {
    return m_size;
  }
  bool operator [] ( int64_t i ) const {
    return (*m_x)[i] || (*m_y)[m_rankx0(i)];
  }

private:
  sdsl::bit_vector const* m_x;
  sdsl::bit_vector const* m_y;
  sdsl::rank_support_v5<0> m_rankx0;
  uint64_t m_size;
};

/// class TwoCache   (helper class for fast_bitvector)
///   simple caching mechanism (stores last two distinct results)
template<class K, class T>
class TwoCache {
public:
  TwoCache () {}
  TwoCache (K key, T data) : key1(key), data1(data), key2(key), data2(data){}
  void assign (K key, T data) {*this = TwoCache(key,data);}
  void put(K key, T data) {
    std::swap ( key1, key2 );
    std::swap ( data1, data2 );
    key2 = key;
    data2 = data;
  }
  int count ( K key ) const {
    if ( key1 == key || key2 == key ) return 1;
    return 0;
  }
  T const& get(K key) {
    if ( key == key1 ) return data1;
    return data2;    
  }
private:
  K key1; T data1;
  K key2; T data2;
};

/// fast_bitvector
///   Same as slow_bitvector, except we provide also the method
///   *  read_word
///   which provides a faster method for reading a 64-bit word
class fast_bitvector {
public:
  fast_bitvector ( sdsl::bit_vector const* x, sdsl::bit_vector const* y ) : m_x(x),m_y(y) {
    m_rankx0 = sdsl::rank_support_v5<0>(m_x); 
    m_size = m_x -> size ();   
    m_block_cache . assign ( 0, read_block ( 0 ) );
    m_word_cache . assign ( 0, read_block ( 0 ) );
  }
  uint64_t size ( void ) const {
    return m_size;
  }
  bool operator [] ( int64_t i ) const {
    uint64_t bit = i & 0x3F;
    uint64_t block = (i & 0xFFFFFFFFFFFFFFC0);
    if ( m_block_cache . count ( block ) == 0 ) m_block_cache . put ( block, read_block ( block ) );
    return ( m_block_cache . get ( block ) & ( ((uint64_t)1) << bit )) != 0;
  }
  uint64_t read_word(int64_t word_pos) const {
    if ( m_word_cache . count ( word_pos ) == 0 ) { 
      int64_t bit_pos = word_pos<<6;
      uint64_t word = read_block ( bit_pos );
      m_word_cache . put ( word_pos, word );
      return word;
    } else {
      return m_word_cache . get ( word_pos ) ;
    }
  }
private:
  uint64_t read_block(int64_t bit_pos) const {
    uint64_t r = m_rankx0(bit_pos);
    uint64_t bit = 1;
    uint64_t end = std::min(bit_pos + 64, (int64_t)m_size);
    uint64_t result = 0;
    while(bit_pos < end) {
      if ((*m_x)[bit_pos++] || (*m_y)[r++]) result |= bit;
      bit <<= 1;
    }
    return result;
  }
  sdsl::bit_vector const* m_x;
  sdsl::bit_vector const* m_y;
  sdsl::rank_support_v5<0> m_rankx0;
  uint64_t m_size;
  mutable TwoCache<int64_t, uint64_t> m_block_cache;
  mutable TwoCache<int64_t, uint64_t> m_word_cache;

};

/// timetrial -- executes function f with each argument between
///   start and end and reports the average time of execution
template <class Functor>
void timetrial ( Functor const& f, uint64_t start, uint64_t end ) 
{
  auto starttime = std::chrono::high_resolution_clock::now();
  for ( uint64_t i = start; i <= end; ++ i ) f(i);
  auto endtime = std::chrono::high_resolution_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>(endtime - starttime);
  std::cout << "Measured time: " << dur.count()/(double)(end-start+1) << "ns per operation\n";
}

/// main function
///   * creates my_bitvector example
///   * creates fast_bitvector example
///   * creates sdsl::bit_vector example
///   * runs checks and time trials
int main ( void ) {
  /// GENERATE RANDOM BITVECTORS
  uint64_t N = 1000000;
  sdsl::bit_vector x ( N );
  uint64_t M = N - sdsl::util::cnt_one_bits(x);
  sdsl::bit_vector y ( M );
  for ( uint64_t i = 0; i < N; ++ i ) x [ i ] = rand () % 2;
  for ( uint64_t i = 0; i < M; ++ i ) y [ i ] = rand () % 2;

  /// USER BITVECTOR EXAMPLE
  typedef sdsl::user_bitvector<slow_bitvector> slow_user_bitvector;
  typedef sdsl::user_bitvector<fast_bitvector> fast_user_bitvector;
  slow_user_bitvector slow ( slow_bitvector ( &x, &y ) );
  fast_user_bitvector fast ( fast_bitvector ( &x, &y ) );
  sdsl::select_support_mcl<1,1,slow_user_bitvector> select_slow ( &slow );
  sdsl::select_support_mcl<1,1,fast_user_bitvector> select_fast ( &fast );

  /// USUAL SDSL::BITVECTOR EXAMPLE
  sdsl::bit_vector z ( N );
  for ( uint64_t i = 0; i < N; ++ i ) z[i] = fast[i];
  sdsl::select_support_mcl<> select_sdsl ( &z );

  /// TEST sdsl::util::cnt_one_bits 
  int64_t L = sdsl::util::cnt_one_bits(z);
  if ( L != sdsl::util::cnt_one_bits(slow) || 
       L != sdsl::util::cnt_one_bits(fast) ) {
    throw std::logic_error("sdsl::util_cnt_one_bits reported the wrong answer\n");
  }

  /// TEST TO CHECK BITVECTOR CONSTRUCTED PROPERLY
  std::cout << "Testing access of user bit-vector.\n";
  timetrial ( [&] (uint64_t i) { 
    if ( z[i] != slow[i] || z[i] != fast[i] ) { 
      throw std::logic_error("Implementation error. vectors not same!\n");
    } }, 0, N-1 );

  /// TEST TO MAKE SURE RESULTS MATCH
  std::cout << "Ensuring results of user bit-vector and sdsl::bit_vector are same.\n";
  timetrial ( [&] (uint64_t i) { 
      if ( select_sdsl(i) != select_slow(i) || select_sdsl(i) != select_fast(i) ) { 
        throw std::logic_error("Select not implemented properly.\n");
      } }, 1, L );
  std::cout << "Passed test.\n";

  /// FORWARD SEQUENTIAL SPEED TEST
  std::cout << "Sequential speed of select_support on slow user bit-vector.\n";
  timetrial ( [&] (uint64_t i) { select_slow(i); }, 1, L );

  std::cout << "Sequential speed of select_support on fast user bit-vector.\n";
  timetrial ( [&] (uint64_t i) { select_fast(i); }, 1, L );

  std::cout << "Sequential speed of select_support on sdsl::bit_vector\n";
  timetrial ( [&] (uint64_t i) { select_sdsl(i); }, 1, L );

  /// BACKWARD SEQUENTIAL SPEED TEST
  std::cout << "Backwards sequential speed of select_support on slow user bit-vector\n";
  timetrial ( [&] (uint64_t i) { select_slow(L - i + 1); }, 1, L );

  std::cout << "Backwards sequential speed of select_support on fast user bit-vector.\n";
  timetrial ( [&] (uint64_t i) { select_fast(L - i + 1); }, 1, L );

  std::cout << "Backwards sequential speed of select_support on sdsl::bit_vector\n";
  timetrial ( [&] (uint64_t i) { select_sdsl(L - i + 1); }, 1, L );

  /// RANDOM ACCESS SPEED TEST
  std::cout << "Random access speed of select_support on slow user bit-vector\n";
  timetrial ( [&] (uint64_t i) { select_slow( 1 + rand()%L ); }, 1, L );

  std::cout << "Random access speed of select_support on fast user bit-vector.\n";
  timetrial ( [&] (uint64_t i) { select_fast( 1 + rand()%L ); }, 1, L );

  std::cout << "Random access speed of select_support on sdsl::bit_vector\n";
  timetrial ( [&] (uint64_t i) { select_sdsl( 1 + rand()%L ); }, 1, L );

  return 0;
}


