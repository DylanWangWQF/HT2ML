#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/time.h>
#include <chrono>

#include <openenclave/host.h>

#include "hteva_u.h"

using namespace std;

oe_enclave_t* enclave = NULL;

int main(int argc, const char * argv[]) {
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    oe_result_t result;
    int ret = 0;
    const uint32_t flags = OE_ENCLAVE_FLAG_DEBUG_AUTO;

    /*---------------------------------------*/
    //  Create the enclave
    /*---------------------------------------*/
    cout << "Host: create enclave for image:" << argv[1] << endl;
    result = oe_create_hteva_enclave(argv[1], OE_ENCLAVE_TYPE_SGX, flags, NULL, 0, &enclave);
    if (result != OE_OK)
    {
        cerr << "oe_create_hteva_enclave() failed with " << argv[0] << " " << result << endl;
        ret = 1;
        goto exit;
    }

    /*-------------------------------------------------------------------*/
    //  Call into the enclave for classification
    /*-------------------------------------------------------------------*/
    result = MatInv(enclave, &ret);
    if (result != OE_OK){
        cerr << "Host: calling into enclave for classification failed. OE result = " << result << endl;
        ret = 1;
        goto exit;
    }

exit:
    cout << "Host: terminate the enclave" << endl;
    oe_terminate_enclave(enclave);
    cout << "Host: done  " << ((ret == 0) ? "succeeded" : "failed") << endl;
    return ret;
}
