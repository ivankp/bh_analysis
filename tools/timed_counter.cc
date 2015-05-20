#include "timed_counter.hh"

#include <iostream>
#include <iomanip>

using namespace std;

timed_counter::timed_counter(bool newline)
: seconds(0), minutes(0), hours(0), last_time( time(0) ), newline(newline)
{ }

void timed_counter::prt(const num_t& ent) const noexcept {
  cout << setw(12) << ent << " | ";
  if (hours) {
    cout << setw(5) << hours << ':'
    << setfill('0') << setw(2) << minutes << ':'
    << setw(2) << seconds << setfill(' ');
  } else if (minutes) {
    cout << setw(2) << minutes << ':'
    << setfill('0') << setw(2) << seconds << setfill(' ');
  } else {
    cout << setw(2) << seconds <<'s';
  }
}

void timed_counter::operator()(const num_t& ent) noexcept {
  cur_time = time(0);
  if ( difftime(cur_time,last_time) > 1 ) {
    prt(ent);
    if (newline) cout << endl;
    else {
      cout.flush();
      if (hours)        for (char i=0;i<26;++i) cout << '\b';
      else if (minutes) for (char i=0;i<20;++i) cout << '\b';
      else              for (char i=0;i<18;++i) cout << '\b';
    }

    last_time = cur_time;
    ++seconds;
    if (seconds==60) {
      seconds = 0;
      ++minutes;
      if (minutes==60) {
        minutes = 0;
        ++hours;
      }
    }

  }
}
