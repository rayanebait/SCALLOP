#include <string.h>
#include <stdio.h>
#include <errno.h>
#include "error.h"

static char buf_err[256];

void
error_sys_set(const char *func, const char *sys)
{
    sprintf(buf_err, "%s(): %s(): [%d, %s]", func, sys, errno, strerror(errno));
}

void
error_set(const char* func, const char* fmt, ...)
{
    // memcpy(buf_err, fmt, );

}

const char*
error_sys_get(void)
{
    return buf_err;
}