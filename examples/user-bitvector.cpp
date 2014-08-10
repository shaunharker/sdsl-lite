#include <iostream>
#include <utility>
#include <algorithm>
#include <chrono>

#include "sdsl/bits.hpp"
#include "sdsl/user_bitvector.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/rank_support_v5.hpp"

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

int main ( void ) {
  uint64_t N = 1000000;
  sdsl::bit_vector x ( N );
  sdsl::bit_vector y ( N );
  for ( uint64_t i = 0; i < N; ++ i ) {
    x [ i ] = rand () % 2;
    y [ i ] = rand () % 2;
  }
  my_bitvector z ( x, y );

#if 0
  std::cout << "x: ";
  for ( uint64_t i = 0; i < N; ++ i ) {
    std::cout << ( x[i] ? "1" : "0" );
  }
  std::cout << "\n";

  std::cout << "y: ";
  for ( uint64_t i = 0; i < N; ++ i ) {
    std::cout << ( y[i] ? "1" : "0" );
  }
  std::cout << "\n";

  sdsl::select_support_mcl<1,1> selecttest ( &x );

  std::cout << "z: ";
  for ( uint64_t i = 0; i < N; ++ i ) {
    std::cout << ( z[i] ? "1" : "0" );
    if ( z[i] ) ++ count;
  }
#endif
  int64_t count = 0;
  for ( uint64_t i = 0; i < N; ++ i ) {
    if ( z[i] ) ++ count;
  }
  //std::cout << "\n count = " << count << "\n";


  /// USER BITVECTOR EXAMPLE
  typedef sdsl::user_bitvector<my_bitvector> custom_bitvector;
  custom_bitvector bv ( std::move(z) );
  sdsl::select_support_mcl<1,1,custom_bitvector> select ( &bv );

  /// USUAL SDSL::BITVECTOR EXAMPLE
  sdsl::bit_vector w ( N );
  for ( uint64_t i = 0; i < N; ++ i ) {
    w[i] = bv[i];
  }
  sdsl::select_support_mcl<> selectw ( &w );

  /// TEST FOR CORRECTNESS
  for ( uint64_t i = 0; i < N; ++ i ) {
    if ( w[i] != bv[i] ) {
      std::cout << "Implementation error. vectors not same!\n";
      abort ();
    }
  }
  for ( uint64_t i = 1; i <= count; ++ i ) {
    if ( select(i) != selectw(i) ) {
      std::cout << "Implementation error. select(" << i 
        << ") = " << select(i) << " and selectw("<<i<<")=" << selectw(i) << "\n";
      std::cout << "count = " << count << "\n";
      std::cout << bv[select(i)] << " " << w[selectw(i)] << "\n";
      std::cout << bv[selectw(i)] << " " << w[select(i)] << "\n";

      abort ();
    }
  }
  std::cout << "Passed test.\n";

  /// FORWARD SEQUENTIAL SPEED TEST
  std::cout << "Sequential speed of select_support on user bit-vector.\n";
  auto start = std::chrono::high_resolution_clock::now();
  for ( uint64_t i = 1; i <= count; ++ i ) {
    select ( i );
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  std::cout << "Measured time: " << dur.count()/(double)N << "ns per select operation\n";

  std::cout << "Sequential speed of select_support on sdsl::bit_vector\n";
  start = std::chrono::high_resolution_clock::now();
  for ( uint64_t i = 1; i <= count; ++ i ) {
    selectw ( i );
  }
  end = std::chrono::high_resolution_clock::now();
  dur = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  std::cout << "Measured time: " << dur.count()/(double)N << "ns per select operation\n";


  /// BACKWARD SEQUENTIAL SPEED TEST
  std::cout << "Backwards sequential speed of select_support on user bit-vector\n";
  start = std::chrono::high_resolution_clock::now();
  for ( uint64_t i = count; i >= 1; -- i ) {
    select ( i );
  }
  end = std::chrono::high_resolution_clock::now();
  dur = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  std::cout << "Measured time: " << dur.count()/(double)N << "ns per select operation\n";

  std::cout << "Backwards sequential speed of select_support on sdsl::bit_vector\n";
  start = std::chrono::high_resolution_clock::now();
  for ( uint64_t i = count; i >= 1; -- i ) {
    selectw ( i );
  }
  end = std::chrono::high_resolution_clock::now();
  dur = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  std::cout << "Measured time: " << dur.count()/(double)N << "ns per select operation\n";

  /// RANDOM ACCESS SPEED TEST
  std::cout << "Random access speed of select_support on user bit-vector\n";
  start = std::chrono::high_resolution_clock::now();
  for ( uint64_t i = count; i >= 1; -- i ) {
    select ( 1 + rand()%count );
  }
  end = std::chrono::high_resolution_clock::now();
  dur = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  std::cout << "Measured time: " << dur.count()/(double)N << "ns per select operation\n";

  std::cout << "Random access speed of select_support on sdsl::bit_vector\n";
  start = std::chrono::high_resolution_clock::now();
  for ( uint64_t i = count; i >= 1; -- i ) {
    selectw ( 1 + rand()%count );
  }
  end = std::chrono::high_resolution_clock::now();
  dur = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  std::cout << "Measured time: " << dur.count()/(double)N << "ns\n";

  return 0;
}
