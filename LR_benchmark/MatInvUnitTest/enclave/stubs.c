#include <openenclave/enclave.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

static const uint64_t _SEC_TO_MSEC = 1000UL;
static const uint64_t _MSEC_TO_USEC = 1000UL;
static const uint64_t _MSEC_TO_NSEC = 1000000UL;

uint64_t oe_get_time(void);

long __wrap_oe_SYS_clock_gettime_impl(long arg1, long arg2)
{
    clockid_t clock_id = (clockid_t)arg1;
    struct timespec* tp = (struct timespec*)arg2;
    int ret = -1;
    uint64_t msec;

    if (!tp)
        goto done;

    if (clock_id != CLOCK_REALTIME)
    {
        /* Only supporting CLOCK_REALTIME */
        // oe_assert("clock_gettime(): panic" == NULL);
        // goto done;
    }

    if ((msec = oe_get_time()) == (uint64_t)-1)
        goto done;

    tp->tv_sec = msec / _SEC_TO_MSEC;
    tp->tv_nsec = (msec % _SEC_TO_MSEC) * _MSEC_TO_NSEC;

    ret = 0;

done:

    return ret;
}