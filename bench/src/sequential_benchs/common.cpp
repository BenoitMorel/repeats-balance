#include "common.h"


#include <stdarg.h>
#include <search.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <iostream>
#define RATE_CATS 4

  
Timer::Timer() {
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start); 
}

long Timer::get_time() {
  timespec end;
  timespec temp;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return (temp.tv_sec * 1000000000 + temp.tv_nsec) / 1000000; 
}




