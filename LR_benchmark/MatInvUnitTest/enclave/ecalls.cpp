#include <openenclave/enclave.h>
#include <cstdio>
#include "hteva_t.h"
#include "dispatcher.h"

static ecall_dispatcher dispatcher;

int MatInv()
{
    return dispatcher.MatInv();
}

void close_encryptor()
{
    return dispatcher.close();
}

// Fix OE's locale implementation which returns NULL.
// Default locale is C.
#include <locale.h>
static char _locale[256] = "C";
extern "C" char* __wrap_setlocale(int category, const char* locale)
{
    OE_UNUSED(category);
    if (locale == NULL)
        return _locale;
    sprintf(_locale, "%s", locale);
    return _locale;
}