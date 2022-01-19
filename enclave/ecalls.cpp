#include <openenclave/enclave.h>
#include <cstdio>
#include "fosdsc_t.h"
#include "dispatcher.h"

static ecall_dispatcher dispatcher;

int enclave_init(uint8_t* hecontext, size_t context_len, uint8_t** rawCtxt, size_t* rawCtxt_len, size_t* num_rawCtxt)
{
    return dispatcher.enclave_init(hecontext, context_len, rawCtxt, rawCtxt_len, num_rawCtxt);
}

// int CtxtTransform(uint8_t* ectxt, uint8_t** octxt, size_t ectxt_len, size_t* octxt_len)
// {
//     return dispatcher.CtxtTransform(ectxt, octxt, ectxt_len, octxt_len);
// }

// int multipleCtxtsTransform(uint8_t* hostCtxt[], size_t hostCtxt_length[], size_t num_resCtxt, uint8_t*** octxt, size_t** octxt_lengths, size_t* octxt_len)
// {
// 	return dispatcher.multipleCtxtsTransform(hostCtxt, hostCtxt_length, num_resCtxt, octxt, octxt_lengths, octxt_len);
// }

int multipleCtxtsTransform(uint8_t* ectxt, size_t ectxt_len, size_t num_ectxt, size_t batchIdx, uint8_t** octxt, size_t* octxt_len, size_t* num_octxt)
{
	return dispatcher.multipleCtxtsTransform(ectxt, ectxt_len, num_ectxt, batchIdx, octxt, octxt_len, num_octxt);
}

int HT_Classify(uint8_t* EvaCtxt, size_t EvaCtxt_len, size_t num_EvaCtxt)
{
    return dispatcher.HT_Classify(EvaCtxt, EvaCtxt_len, num_EvaCtxt);
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