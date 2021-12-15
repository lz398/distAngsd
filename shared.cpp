
#include <cassert>
#include "shared.h"

void my_bgzf_write(BGZF *fp, const void *data, size_t length){
  assert(bgzf_write(fp,data,length)==length);
}
