#include <openenclave/enclave.h>
#include <cstdio>
#include "fosdsc_t.h"
#include "dispatcher.h"

static ecall_dispatcher dispatcher;

int enclave_init(uint8_t* hecontext, size_t context_len)
{
    return dispatcher.enclave_init(hecontext, context_len);
}

int MatrixOperation(uint8_t* ectxt, size_t ectxt_len, uint8_t** octxt, size_t* octxt_len)
{
	return dispatcher.MatrixOperation(ectxt, ectxt_len, octxt, octxt_len);
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