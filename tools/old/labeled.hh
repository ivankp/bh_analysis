//---------------------------
// Developer: Ivan Pogrebnyak
//---------------------------

#ifndef labeled_h
#define labeled_h

#include <string>
#include <functional>

template<class T, class Label=std::string, class Compare=std::less<Label>>
struct labeled {
  typedef T val_t;
  typedef Label label_t;

  val_t val;
  label_t label;

  // labeled(const val_t& val, const label_t& label)
  // : val(val), label(label) { }
  // virtual ~labeled() { }

  bool operator<(const labeled<T,Label,Compare>& other) const {
    return Compare(*this,other);
  }
};

#endif
