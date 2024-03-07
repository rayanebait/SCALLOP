#ifndef ERROR_H
#define ERROR_H

#include <errno.h>


void error_sys_set(const char*, const char*);
void error_set(const char*, const char*,...);
const char *error_sys_get(void);





#endif /*ERROR_H*/