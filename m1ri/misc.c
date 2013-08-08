#include <stdio.h>
#include <stdarg.h>
#include "misc.h"

void m1ri_die(const char *errormessage, ...) {
  va_list lst;
  va_start(lst, errormessage);
  vfprintf(stderr, errormessage, lst);
  va_end(lst);
  abort();
}
