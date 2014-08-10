#include <iostream>
#include <utility>
#include <algorithm>
#include <chrono>
#include <stdexcept>

#include "sdsl/util.hpp"
#include "sdsl/user_bitvector.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/rank_support_v5.hpp"

/// class TwoCache
///   simple caching mechanism (stores last two distinct results)
template<class K, class T>
class TwoCache {
public:
  TwoCache ( void ) {}
  TwoCache ( K key1, T data1,
             K key2, T data2 ) 
  : key1(key1), data1(data1), key2(key2), data2(data2){}
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
/// We want to make a bitvector z out of bitvectors x and y
/// such that:
/// length(y) == rank0(x) (i.e. total number of 0's in x)
/// z[i] := 1 whenever x[i] == 1
/// z[i] := y[rank0[i]] whenever x[i] == 0
class my_bitvector {
public:
  my_bitvector ( sdsl::bit_vector const& x,
		 sdsl::bit_vector const& y ) 
  : x(x),y(y) {
    rankx0 = sdsl::rank_support_v5<0>(&x); 
    size_ = x . size ();   
    block_cache_ . put ( 0, read_block ( 0 ) );
    block_cache_ . put ( 0, read_block ( 0 ) ); 
    word_cache_ . put ( 0, read_block ( 0 ) );
    word_cache_ . put ( 0, read_block ( 0 ) );  
  }
  uint64_t size ( void ) const {
    return size_;
  }
  bool operator [] ( int64_t i ) const {
    //return x[i] || y[rankx0(i)];
    uint64_t bit = i & 0x3F;
    uint64_t block = (i & 0xFFFFFFFFFFFFFFC0);
    if ( block_cache_ . count ( block ) == 0 ) block_cache_ . put ( block, read_block ( block ) );
    return ( block_cache_ . get ( block ) & ( ((uint64_t)1) << bit )) != 0;
  }
  uint64_t read_word(int64_t word_pos) const {
    if ( word_cache_ . count ( word_pos ) == 0 ) { 
      int64_t bit_pos = word_pos<<6;
      uint64_t word = read_block ( bit_pos );
      word_cache_ . put ( word_pos, word );
      return word;
    } else {
      return word_cache_ . get ( word_pos ) ;
    }
  }
  uint64_t read_block(int64_t bit_pos) const {
    uint64_t r = rankx0(bit_pos);
    uint64_t bit = 1;
    uint64_t end = std::min(bit_pos + 64, (int64_t)size_);
    uint64_t result = 0;
    while(bit_pos < end) {
      if (x[bit_pos++] || y[r++]) result |= bit;
      bit <<= 1;
    }
    return result;
  }
private:
  sdsl::bit_vector x;
  sdsl::bit_vector y;
  sdsl::rank_support_v5<0> rankx0;
  uint64_t size_;
  mutable TwoCache<int64_t, uint64_t> block_cache_;
  mutable TwoCache<int64_t, uint64_t> word_cache_;

};

template <class Functor>
void timetrial ( Functor const& f, uint64_t start, uint64_t end ) 
{
  auto starttime = std::chrono::high_resolution_clock::now();
  for ( uint64_t i = start; i <= end; ++ i ) f(i);
  auto endtime = std::chrono::high_resolution_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>(endtime - starttime);
  std::cout << "Measured time: " << dur.count()/(double)(end-start+1) << "ns per operation\n";
}

int main ( void ) {
  /// GENERATE RANDOM BITVECTORS
  uint64_t N = 1000000;
  sdsl::bit_vector x ( N );
  uint64_t M = N - sdsl::util::cnt_one_bits(x);
  sdsl::bit_vector y ( M );
  for ( uint64_t i = 0; i < N; ++ i ) x [ i ] = rand () % 2;
  for ( uint64_t i = 0; i < M; ++ i ) y [ i ] = rand () % 2;
  my_bitvector z ( x, y );

  /// USER BITVECTOR EXAMPLE
  typedef sdsl::user_bitvector<my_bitvector> custom_bitvector;
  custom_bitvector bv ( std::move(z) );
  sdsl::select_support_mcl<1,1,custom_bitvector> select ( &bv );

  /// USUAL SDSL::BITVECTOR EXAMPLE
  sdsl::bit_vector w ( N );
  for ( uint64_t i = 0; i < N; ++ i ) w[i] = bv[i];
  sdsl::select_support_mcl<> selectw ( &w );

  /// TEST TO MAKE SURE user_bitvector works with sdsl::util
  int64_t L = sdsl::util::cnt_one_bits(bv);

  /// TEST TO CHECK BITVECTOR CONSTRUCTED PROPERLY
  std::cout << "Testing access of user bit-vector.\n";
  timetrial ( [&] (uint64_t i) { 
    if ( w[i] != bv[i] ) { 
      throw std::logic_error("Implementation error. vectors not same!\n");
    } }, 0, N-1 );

  /// TEST TO MAKE SURE RESULTS MATCH
  std::cout << "Ensuring results of user bit-vector and sdsl::bit_vector are same.\n";
  timetrial ( [&] (uint64_t i) { 
      if ( select(i) != selectw(i) ) { 
        throw std::logic_error("Select not implemented properly.\n");
      } }, 1, L );
  std::cout << "Passed test.\n";

  /// FORWARD SEQUENTIAL SPEED TEST
  std::cout << "Sequential speed of select_support on user bit-vector.\n";
  timetrial ( [&] (uint64_t i) { select(i); }, 1, L );

  std::cout << "Sequential speed of select_support on sdsl::bit_vector\n";
  timetrial ( [&] (uint64_t i) { selectw(i); }, 1, L );

  /// BACKWARD SEQUENTIAL SPEED TEST
  std::cout << "Backwards sequential speed of select_support on user bit-vector\n";
  timetrial ( [&] (uint64_t i) { select(L - i + 1); }, 1, L );

  std::cout << "Backwards sequential speed of select_support on sdsl::bit_vector\n";
  timetrial ( [&] (uint64_t i) { selectw(L - i + 1); }, 1, L );

  /// RANDOM ACCESS SPEED TEST
  std::cout << "Random access speed of select_support on user bit-vector\n";
  timetrial ( [&] (uint64_t i) { select( 1 + rand()%L ); }, 1, L );

  std::cout << "Random access speed of select_support on sdsl::bit_vector\n";
  timetrial ( [&] (uint64_t i) { selectw( 1 + rand()%L ); }, 1, L );

  return 0;
}


