#include <string.h>
#include <NTL/matrix.h>
#include <vector>
#include <chrono>

#include "dispatcher.h"
#include "trace.h"

//HTparams* params
ecall_dispatcher::ecall_dispatcher()
{
    // e_HTparams = params;
}

int ecall_dispatcher::enclave_init(uint8_t* hecontext, size_t context_len)
{
    int ret = 0;
    // if (context_len == 0)
    // {
    //     TRACE_ENCLAVE("passed hecontext is null string, context_len =0");
    //     ret = 1;
    //     goto exit;
    // }
    auto start= chrono::steady_clock::now();

    // TRACE_ENCLAVE("Enclave: receive the HE params");
    stringstream ess;
    string e_hecontext = string(hecontext, hecontext + context_len);
    // cout << "Enclave: check hecontext: " << hecontext << endl;
    // cout << "Enclave: check context_len, containing context, sk and pk: " << context_len << endl;
    ess << e_hecontext;
    // TRACE_ENCLAVE("Enclave: readfrom the string and reconstruct the context: ");
    e_context = Context::readPtrFrom(ess);
    // TRACE_ENCLAVE("Enclave: Setup Seckey by reading from ss!");
    activeSecKey = make_unique<SecKey>(SecKey::readFrom(ess, *e_context));
    // TRACE_ENCLAVE("Enclave: Setup PubKey by reading from ss!");
    activePubKey = make_unique<PubKey>(PubKey::readFrom(ess, *e_context));

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Enclave: Total Setup Time for HE inside enclave = " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;

exit:
    TRACE_ENCLAVE("Enclave: free memory in enclave_init()!");
    e_hecontext.shrink_to_fit();
    return ret;
}

int ecall_dispatcher::MatrixOperation(uint8_t* ectxt, size_t ectxt_len, uint8_t** octxt, size_t* octxt_len)
{
    auto start= chrono::steady_clock::now();
    stringstream css;

    // TRACE_ENCLAVE("Enclave: Receiving Ctxt start!");
    string etemp = string(ectxt, ectxt + ectxt_len);
    css << etemp;
    vector<long> cmsg1;
    int*** HostMatrix = new int**[2];
    for (int k = 0; k < 2; k++)
    {
        // initialize matrix using int*** format
        HostMatrix[k] = new int*[MatrixDim];
        for (int l = 0; l < MatrixDim; l++)
        {
            HostMatrix[k][l] = new int[MatrixDim];
        }

        // read two ctxts from string buffer
        cmsg1.clear();
        Ctxt ctemp(*activePubKey);
        ctemp.Ctxt::read(css);
        e_context->getEA().decrypt(ctemp, *activeSecKey, cmsg1);

        // assign the values in matrix
        int idx = 0;
        for(int i = 0; i < MatrixDim; ++i){
            for(int j = 0; j < MatrixDim; ++j){
                HostMatrix[k][i][j] = cmsg1[idx];
                idx++;
            }
        }
    }

    // HostMatrix[0] for inverse
    int** inverse_result = new int[MatrixDim];
    for (int l = 0; l < MatrixDim; l++)
    {
        inverse_result[l] = new int[MatrixDim];
    }

    // multiply the inverse by HostMatrix[1]
    int** final_result = new int[MatrixDim];
    for (int l = 0; l < MatrixDim; l++)
    {
        final_result[l] = new int[MatrixDim];
    }

    // send enclaveCtxt to server
    // TRACE_ENCLAVE("Enclave: Transforming Ctxt start!");
    vector<long> cmsg2(MatrixDim * MatrixDim);
    for(int i = 0; i < MatrixDim; ++i){
        for(int j = 0; j < MatrixDim; ++j){
            cmsg2[i * MatrixDim + j] = final_result[i][j];
        }
    }
    Ctxt ctemp(*activePubKey);
    e_context->getEA().encrypt(ctemp, *activePubKey, cmsg2);

    css.str(std::string());
    css.clear();
    ctemp.writeTo(css);
    string otemp = css.str();
    uint8_t* host_buf = (uint8_t*) oe_host_malloc(size_t(otemp.size() + 1));
    memcpy(host_buf, (uint8_t*)otemp.c_str(), otemp.size() + 1);
    *octxt = host_buf;
    *octxt_len = otemp.size();

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Enclave: Matrix Inverse Runtime inside enclave = " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;

    return 0;
}

void ecall_dispatcher::close()
{
    TRACE_ENCLAVE("Enclave: release context and HTparams!");
    // release context and keys
    delete e_context;
    TRACE_ENCLAVE("ecall_dispatcher::close");
}

int** ecall_dispatcher::matMul(int** A, int** B){
    int** C = new int*[MatrixDim];
    for(int i=0 ; i < MatrixDim; i++) C[i] = new int[MatrixDim];

    for(int i = 0; i < MatrixDim; i++)
    {
        for(int j = 0; j < MatrixDim; j++)
        {
            for(int k = 0; k < MatrixDim; k++)
            {
                C[i * MatrixDim + j] += A[ i * MatrixDim + k] * B[k * MatrixDim + j];
            }
        }
    }
    return C;
}

void ecall_dispatcher::LUP_Descomposition(int**& A, int**& L, int**& U, int*& P){
    int row = 0;
    for(int i = 0; i < MatrixDim; i++)
    {
        P[i] = i;
    }
    for(int i = 0; i < (MatrixDim - 1); i++)
    {
        int p = 0;
        for(int j = i; j < MatrixDim; j++)
        {
            if(abs(A[j][i]) > p)
            {
                p = abs(A[j][i]);
                row = j;
            }
        }
        if(0 == p)
        {
            cout<< "Unable to calculate matrix inverse due to singular matrix" <<endl;
            return;
        }

        // exchange P[i] and P[row]
        int tmp = P[i];
        P[i]  = P[row];
        P[row] = tmp;

        int tmp2 = 0;
        for(int j = 0; j < MatrixDim; j++)
        {
            // exchange A[i][j] and A[row][j]
            tmp2 =A[i][j];
            A[i][j] = A[row][j];
            A[row][j] = tmp2;
        }

        // LUP_Descomposition 
        int u = A[i][i], l = 0;
        for(int j = i + 1; j < MatrixDim; j++)
        {
            l = A[j*N+i]/u;
            A[j*N+i]=l;
            for(int k=i+1;k<N;k++)
            {
                A[j*N+k]=A[j*N+k]-A[i*N+k]*l;
            }
        }

    }

    //构造L和U
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<=i;j++)
        {
            if(i!=j)
            {
                L[i*N+j]=A[i*N+j];
            }
            else
            {
                L[i*N+j]=1;
            }
        }
        for(int k=i;k<N;k++)
        {
            U[i*N+k]=A[i*N+k];
        }
    }
}