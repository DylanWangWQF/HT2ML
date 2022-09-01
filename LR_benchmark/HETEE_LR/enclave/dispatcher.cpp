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
    // TRACE_ENCLAVE("Enclave: free memory in enclave_init()!");
    e_hecontext.shrink_to_fit();
    return ret;
}

int ecall_dispatcher::MatrixOperation(uint8_t* ectxt, size_t ectxt_len, uint8_t** octxt, size_t* octxt_len)
{
    // TRACE_ENCLAVE("Enclave: Receive Ctxt and calculate matrix inverse!");
    stringstream css;
    auto start= chrono::steady_clock::now();

    string etemp = string(ectxt, ectxt + ectxt_len);
    css << etemp;
    vector<long> cmsg1;
    int HostMatrix[2][MatrixDim][MatrixDim];
    for (int k = 0; k < 2; k++)
    {
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


    // Crop the matrices by removing zeros
    // TRACE_ENCLAVE("Enclave: Crop the matrices!");

    double *A = new double[numAttr*numAttr]();
    // int RealMatrix1[numAttr][numAttr]; // 6 * 6
    int RealMatrix2[numAttr][1];  // 6 * 1
    for (int i = 0; i < numAttr; i++)
    {
        for (int j = 0; j < numAttr; j++)
        {
            // RealMatrix1[i][j] = HostMatrix[0][i][j];
            A[i * numAttr + j] = HostMatrix[0][i][j];
            if (j == 0)
            {
                RealMatrix2[i][j] = HostMatrix[1][i][j];
            }
            
        }
    }

    // HostMatrix[0] for inverse
    // TRACE_ENCLAVE("Enclave: calculate the inverse!");
    auto inverse_start= chrono::steady_clock::now();

    // Adjoint
    // int inverse_result[numAttr][numAttr];
    // if (!inverse(RealMatrix1, inverse_result))
    // {
    //     cout << "Can't find its inverse" << endl;  
    // }

    // LUP Decomposition
    double *invOfA = new double[numAttr*numAttr]();
    invOfA = LUP_solve_inverse(A);
    

    auto inverse_end= chrono::steady_clock::now();
    auto inverse_diff = inverse_end - inverse_start;
    // double inverse_timeElapsed = chrono::duration <double, micro> (inverse_diff).count();
    // double inverse_timeElapsed = chrono::duration <double, nano> (inverse_diff).count();
    double inverse_timeElapsed = chrono::duration <double, milli> (inverse_diff).count();

    int inverse_result[numAttr][numAttr];
    for (int i = 0; i < numAttr; i++)
    {
        for (int j = 0; j < numAttr; j++)
        {
            inverse_result[i][j] = (int) invOfA[i * numAttr + j];
        }
    }

    // multiply the inverse by HostMatrix[1]
    // TRACE_ENCLAVE("Enclave: multiply the inverse!");
    int final_result[numAttr][1];
    for(int i = 0; i < numAttr; i++)
    {
        for(int j = 0; j < 1; j++)
        {
            for(int k = 0; k < numAttr; k++)
            {
                final_result[i][j] += inverse_result[i][k] * RealMatrix2[k][j];
            }
        }
    }

    // send enclaveCtxt to server
    // TRACE_ENCLAVE("Enclave: Transforming Ctxt start!");
    vector<long> cmsg2(MatrixDim * MatrixDim, 0);
    // for(int i = 0; i < MatrixDim; ++i){
    //     for(int j = 0; j < MatrixDim; ++j){
    //         cmsg2[i * MatrixDim + j] = final_result[i][j];
    //     }
    // }
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
    cout << "Enclave: Matrix inverse runtime inside enclave = " << inverse_timeElapsed << " ms" << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Enclave: Runtime inside enclave = " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;

    return 0;
}

int ecall_dispatcher::RefreshCtxt(uint8_t* ectxt, size_t ectxt_len, uint8_t** octxt, size_t* octxt_len)
{
    // TRACE_ENCLAVE("Enclave: Receive Ctxt and calculate matrix inverse!");
    stringstream css;
    auto start= chrono::steady_clock::now();

    string etemp = string(ectxt, ectxt + ectxt_len);
    css << etemp;
    vector<vector<long>> cmsg(2);
    int HostMatrix[2][MatrixDim][MatrixDim];
    for (int k = 0; k < 2; k++)
    {
        // read two ctxts from string buffer
        Ctxt ctemp(*activePubKey);
        ctemp.Ctxt::read(css);
        e_context->getEA().decrypt(ctemp, *activeSecKey, cmsg[k]);
    }

    css.str(std::string());
    css.clear();
    for (int k = 0; k < 2; k++)
    {
        // read two ctxts from string buffer
        Ctxt ctemp(*activePubKey);
        e_context->getEA().encrypt(ctemp, *activePubKey, cmsg[k]);
        ctemp.writeTo(css);
    }

    string otemp = css.str();
    uint8_t* host_buf = (uint8_t*) oe_host_malloc(size_t(otemp.size() + 1));
    memcpy(host_buf, (uint8_t*)otemp.c_str(), otemp.size() + 1);
    *octxt = host_buf;
    *octxt_len = otemp.size();

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Enclave: Runtime of refreshing the ctxt inside enclave = " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;

    return 0;
}

// Function to get cofactor of A[p][q] in temp[][]. n is current dimension of A[][]
void ecall_dispatcher::getCofactor(int A[numAttr][numAttr], int temp[numAttr][numAttr], int p, int q, int n){
    int i = 0, j = 0;

	// Looping for each element of the matrix
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			// Copying into temporary matrix only those element which are not in given row and column
			if (row != p && col != q)
			{
				temp[i][j++] = A[row][col];

				// Row is filled, so increase row index and reset col index
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

// Recursive function for finding determinant of matrix n is current dimension of A[][].
int ecall_dispatcher::determinant(int A[numAttr][numAttr], int n){
    int D = 0; // Initialize result

	// Base case : if matrix contains single element
	if (n == 1)
		return A[0][0];

	int temp[numAttr][numAttr]; // To store cofactors

	int sign = 1; // To store sign multiplier

	// Iterate for each element of first row
	for (int f = 0; f < n; f++)
	{
		// Getting Cofactor of A[0][f]
		getCofactor(A, temp, 0, f, n);
		D += sign * A[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

// Function to get adjoint of A[numAttr][numAttr] in adj[numAttr][numAttr].
void ecall_dispatcher::adjoint(int A[numAttr][numAttr],int adj[numAttr][numAttr])
{
	if (numAttr == 1)
	{
		adj[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][]
	int sign = 1, temp[numAttr][numAttr];

	for (int i=0; i<numAttr; i++)
	{
		for (int j=0; j<numAttr; j++)
		{
			// Get cofactor of A[i][j]
			getCofactor(A, temp, i, j, numAttr);

			// sign of adj[j][i] positive if sum of row and column indexes is even.
			sign = ((i+j)%2==0)? 1: -1;

			// Interchanging rows and columns to get the transpose of the cofactor matrix
			adj[j][i] = (sign)*(determinant(temp, numAttr-1));
		}
	}
}

// Function to calculate and store inverse, returns false if matrix is singular
bool ecall_dispatcher::inverse(int A[numAttr][numAttr], int inverse[numAttr][numAttr])
{
	// Find determinant of A[][]
	int det = determinant(A, numAttr);
	if (det == 0)
	{
		cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint
	int adj[numAttr][numAttr];
	adjoint(A, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i=0; i < numAttr; i++)
		for (int j=0; j<numAttr; j++)
			inverse[i][j] = adj[i][j]/det;

	return true;
}

void ecall_dispatcher::LUP_Descomposition(double A[numAttr*numAttr],double L[numAttr*numAttr],double U[numAttr*numAttr],int P[numAttr])
{
	int row = 0;
    for(int i = 0; i < numAttr;i++)
    {
        P[i]=i;
    }
    for(int i=0;i<numAttr-1;i++)
    {
        double p=0.0;
        for(int j=i;j<numAttr;j++)
        {
            if(abs(A[j*numAttr+i])>p)
            {
                p=abs(A[j*numAttr+i]);
                row=j;
            }
        }
        if(0==p)
        {
            cout<< "singular matrix, cannot find inverse" <<endl;
            return ;
        }

        int tmp=P[i];
        P[i]=P[row];
        P[row]=tmp;

        double tmp2=0.0;
        for(int j=0;j<numAttr;j++)
        {
            tmp2=A[i*numAttr+j];
            A[i*numAttr+j]=A[row*numAttr+j];
            A[row*numAttr+j]=tmp2;
        }

        double u=A[i*numAttr+i], l=0.0;
        for(int j=i+1;j<numAttr;j++)
        {
            l=A[j*numAttr+i]/u;
            A[j*numAttr+i]=l;
            for(int k=i+1;k<numAttr;k++)
            {
                A[j*numAttr+k]=A[j*numAttr+k]-A[i*numAttr+k]*l;
            }
        }

    }

    for(int i=0;i<numAttr;i++)
    {
        for(int j=0;j<=i;j++)
        {
            if(i!=j)
            {
                L[i*numAttr+j]=A[i*numAttr+j];
            }
            else
            {
                L[i*numAttr+j]=1;
            }
        }
        for(int k=i;k<numAttr;k++)
        {
            U[i*numAttr+k]=A[i*numAttr+k];
        }
    }
}

double* ecall_dispatcher::LUP_Solve(double L[numAttr*numAttr],double U[numAttr*numAttr],int P[numAttr],double b[numAttr])
{
	double *x=new double[numAttr]();
    double *y=new double[numAttr]();
	
    for(int i = 0;i < numAttr;i++)
    {
        y[i] = b[P[i]];
        for(int j = 0;j < i;j++)
        {
            y[i] = y[i] - L[i*numAttr+j]*y[j];
        }
    }
	
    for(int i = numAttr-1;i >= 0; i--)
    {
        x[i]=y[i];
        for(int j = numAttr-1;j > i;j--)
        {
            x[i] = x[i] - U[i*numAttr+j]*x[j];
        }
        x[i] /= U[i*numAttr+i];
    }
    return x;
}

int ecall_dispatcher::getNext(int i, int m, int n)
{
	return (i%n)*m + i/n;
}

int ecall_dispatcher::getPre(int i, int m, int n)
{
	return (i%m)*n + i/m;
}

void ecall_dispatcher::movedata(double *mtx, int i, int m, int n)
{
	double temp = mtx[i]; 
	int cur = i;   
	int pre = getPre(cur, m, n);
	while(pre != i)
	{
		mtx[cur] = mtx[pre];
		cur = pre;
		pre = getPre(cur, m, n);
	}
	mtx[cur] = temp;
}

void ecall_dispatcher::transpose(double *mtx, int m, int n)
{
	for(int i=0; i<m*n; ++i)
	{
		int next = getNext(i, m, n);
		while(next > i) 
			next = getNext(next, m, n);
		if(next == i)
			movedata(mtx, i, m, n);
	}
}

double* ecall_dispatcher::LUP_solve_inverse(double A[numAttr*numAttr])
{
	double* A_mirror = new double[numAttr*numAttr]();
    double* inv_A = new double[numAttr*numAttr]();
    double* inv_A_each = new double[numAttr]();
    double* b = new double[numAttr]();

    for(int i=0;i<numAttr;i++)
    {
        double *L=new double[numAttr*numAttr]();
        double *U=new double[numAttr*numAttr]();
        int *P=new int[numAttr]();

        for(int i=0;i<numAttr;i++)
        {
            b[i]=0;
        }
        b[i]=1;

        for(int i=0;i<numAttr*numAttr;i++)
        {
            A_mirror[i]=A[i];
        }

        LUP_Descomposition(A_mirror,L,U,P);

        inv_A_each=LUP_Solve (L,U,P,b);
        memcpy(inv_A+i*numAttr,inv_A_each,numAttr*sizeof(double));
    }
    transpose(inv_A,numAttr,numAttr);

    return inv_A;
}

double* ecall_dispatcher::mul(double A[numAttr*numAttr],double B[numAttr*numAttr])
{
	double *C=new double[numAttr*numAttr]{};
    for(int i=0;i<numAttr;i++)
    {
        for(int j=0;j<numAttr;j++)
        {
            for(int k=0;k<numAttr;k++)
            {
                C[i*numAttr+j] += A[i*numAttr+k]*B[k*numAttr+j];
            }
        }
    }

    for(int i=0;i<numAttr*numAttr;i++)
    {
        if(abs(C[i])<pow(10,-10))
        {
            C[i]=0;
        }
    }

    return C;
}

// void ecall_dispatcher::display(int A[MatrixDim][MatrixDim])
// {
// 	for (int i = 0; i < MatrixDim; i++)
// 	{
// 		for (int j = 0; j < MatrixDim; j++)
// 			cout << A[i][j] << " ";
// 		cout << endl;
// 	}
// }

// void ecall_dispatcher::mul(int A[MatrixDim][MatrixDim], int B[MatrixDim][MatrixDim], int C[MatrixDim][MatrixDim])
// {
//     for(int i = 0; i < MatrixDim; i++)
//     {
//         for(int j = 0; j < MatrixDim; j++)
//         {
//             for(int k = 0; k < MatrixDim; k++)
//             {
//                 C[i][j] += A[i][k] * B[k][j];
//             }
//         }
//     }
// }

void ecall_dispatcher::close()
{
    TRACE_ENCLAVE("Enclave: release context and HTparams!");
    // release context and keys
    delete e_context;
    TRACE_ENCLAVE("ecall_dispatcher::close");
}