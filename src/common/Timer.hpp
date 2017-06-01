#include <time.h>

class Timer {
public:
  Timer();
  // ms
  long get_time();

private:
  timespec start;
};





