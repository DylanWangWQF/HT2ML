#include "include/obl/obl_primitives.h"

namespace obl {

// TODO: define this as inline cause segment fault. Need to know why. Is it due
// to |inline| does not work well with |assembly| impl?
bool LessImplDouble(double x, double y) {
  bool result;
  __asm__ volatile(
      "movsd %1, %%xmm0;"
      "movsd %2, %%xmm1;"
      "comisd %%xmm1, %%xmm0;"
      "setb %0;"
      : "=r"(result)
      : "m"(x), "m"(y)
      : "cc");
  return result;
}

bool LessImplFloat(float x, float y) {
    bool result;
    __asm__ volatile(
            "movsd %1, %%xmm0;"
            "movsd %2, %%xmm1;"
            "comiss %%xmm1, %%xmm0;"
            "setb %0;"
            : "=r"(result)
            : "m"(x), "m"(y)
            : "cc");
    return result;
}

}  // namespace obl