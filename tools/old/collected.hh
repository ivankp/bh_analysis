//---------------------------
// Developer: Ivan Pogrebnyak
//---------------------------

#ifndef collected_h
#define collected_h

#include <vector>
#include <memory>

template<class T, class Container=std::vector<std::unique_ptr<const collected<T>
struct collected {
  typedef T val_t;

  val_t val;

  // template <class... Args>
  // collected(Args&&... args): val(args...) { }

  static Container<std::unique_ptr<const collected<T,Container>>> all;

  template <class... Args>
  static void add(Args&&... args) {
    all.emplace_back( new collected<T,Container>(args...) );
  }
};

template<class T, class Container>
Container<std::unique_ptr<const collected<T,Container>>> collected<T,Container>::all;

#endif
