#include <string.h>
#include <vector>
#include <chrono>
#include "dispatcher.h"
#include "trace.h"

//HTparams* params
ecall_dispatcher::ecall_dispatcher()
{
}

int ecall_dispatcher::MatInv()
{
	// int inputMat[numAttr][numAttr] = { {5, -2}, {1, 0}};

    // int inputMat[numAttr][numAttr] = { {5, -2, 2, 7}, {1, 0, 0, 3}, {-3, 1, 5, 0}, {3, -1, -9, 4}};

    // int inputMat[numAttr][numAttr] = 
    // {{6, 6, 9, 5, 6, 4},
	// {2, 9, 6, 3, 8, 1},
	// {2, 3, 3, 9, 1, 9},
	// {4, 0, 6, 2, 4, 7},
	// {7, 9, 4, 3, 2, 48},
	// {4, 3, 7, 3, 2, 0}};

 //    int inputMat[numAttr][numAttr] = 
 //    {{6, 6, 9, 5, 6, 4, 7, 4, 1, 2},
	// {2, 9, 6, 3, 8, 1, 5, 6, 3, 5},
	// {2, 3, 3, 9, 1, 9, 9, 0, 5, 4},
	// {4, 0, 6, 2, 4, 7, 4, 1, 3, 0},
	// {7, 9, 4, 3, 2, 4, 4, 7, 0, 8},
	// {4, 3, 7, 3, 2, 0, 7, 1, 5, 2},
	// {3, 7, 9, 6, 1, 1, 7, 9, 4, 4},
	// {9, 2, 2, 2, 9, 2, 2, 4, 5, 4},
	// {6, 7, 2, 3, 9, 4, 2, 7, 7, 8},
	// {8, 0, 8, 8, 4, 0, 5, 7, 1, 0}};

	// int inputMat[numAttr][numAttr];
	// srand((unsigned)time(0));
    // for(int i = 0; i< numAttr; i++)
    // {
	// 	for(int j = 0; j < numAttr; j++)
	// 	{
	// 		inputMat[i][j] = rand() % 100;
	// 	}
    // }

	// cout << "Generated input matrix: " << endl;
    // for (int i = 0; i < numAttr; i++)
    // {
    //     for (int j = 0; j < numAttr; j++)
    //     {
    //         cout << inputMat[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

	// LUP decomposition
	double *A = new double[numAttr*numAttr]();

    srand((unsigned)time(0));
    for(int i=0; i<numAttr ;i++)
    {
        for(int j=0; j<numAttr;j++)
        {
            A[i*numAttr+j]=rand()%100 *0.01;
        }
    }

    auto inverse_start= chrono::steady_clock::now();

    // int inverse_result[numAttr][numAttr];
    
    // if (!inverse(inputMat, inverse_result))
    // {
    //     cout << "Can't find its inverse" << endl;  
    // }
    // cout << "Final inverse: " << endl;
    // for (int i = 0; i < numAttr; i++)
    // {
    //     for (int j = 0; j < numAttr; j++)
    //     {
    //         cout << inverse_result[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

	double *invOfA = new double[numAttr*numAttr]();
    invOfA = LUP_solve_inverse(A);
	cout << "Final inverse: " << endl;
    for (int i = 0; i < numAttr; i++)
    {
        for (int j = 0; j < numAttr; j++)
        {
            cout << invOfA[i * numAttr + j] << " ";
        }
        cout << endl;
    }
    cout << endl;


    auto inverse_end= chrono::steady_clock::now();
    auto inverse_diff = inverse_end - inverse_start;
    // double inverse_timeElapsed = chrono::duration <double, milli> (inverse_diff).count() / (1000 * 3600);
	double inverse_timeElapsed = chrono::duration <double, milli> (inverse_diff).count() / 1000;

    cout << "------------------------------------------------------------------------" << endl;
    cout << "Enclave: Matrix Inverse Runtime = " << inverse_timeElapsed << " s" << endl;
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
    // cout << "Final adj: " << endl;
    // for (int i = 0; i < numAttr; i++)
    // {
    //     for (int j = 0; j < numAttr; j++)
    //     {
    //         cout << adj[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i=0; i < numAttr; i++)
		for (int j=0; j < numAttr; j++)
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

// void ecall_dispatcher::LUP_Descomposition(int A[numAttr][numAttr], int L[numAttr][numAttr], int U[numAttr][numAttr], int P[numAttr])
// {
// 	int row = 0;
//     for(int i = 0; i < numAttr; i++)
//     {
//         P[i]=i;
//     }

//     for(int i = 0; i < numAttr - 1; i++)
//     {
//         int p = 0;
//         for(int j = i; j < numAttr; j++)
//         {
//             if(abs(A[j][i]) > p)
//             {
//                 p = abs(A[j][i]);
//                 row = j;
//             }
//         }

//         if(0 == p)
//         {
//             cout<< "矩阵奇异，无法计算逆" <<endl;
//             return ;
//         }

//         int tmp = P[i];
//         P[i] = P[row];
//         P[row] = tmp;

//         int tmp2 = 0;
//         for(int j = 0; j < numAttr; j++)
//         {
//             tmp2 = A[i][j];
//             A[i][j] = A[row][j];
//             A[row][j] = tmp2;
//         }

//         int u = A[i][i], l = 0;
//         for(int j = i + 1;j < numAttr; j++)
//         {
//             l = A[j][i] / u;
//             A[j][i] = l;
//             for(int k = i + 1; k < numAttr; k++)
//             {
//                 A[j][k] = A[j][k] - A[i][k] * l;
//             }
//         }

//     }

//     for(int i = 0; i < numAttr; i++)
//     {
//         for(int j = 0; j <= i; j++)
//         {
//             if(i != j)
//             {
//                 L[i][j]=A[i][j];
//             }
//             else
//             {
//                 L[i][j]=1;
//             }
//         }
//         for(int k = i; k < numAttr; k++)
//         {
//             U[i][k] = A[i][k];
//         }
//     }
// }

// void ecall_dispatcher::LUP_Solve(int L[numAttr][numAttr], int U[numAttr][numAttr], int P[numAttr], double b[numAttr], int x[numAttr])
// {
// 	int y[numAttr];

//     for(int i = 0; i < numAttr; i++)
//     {
//         y[i] = b[P[i]];
//         for(int j = 0;j < i;j++)
//         {
//             y[i] = y[i] - L[i][j] * y[j];
//         }
//     }
	
//     for(int i = numAttr - 1; i >= 0; i--)
//     {
//         x[i] = y[i];
//         for(int j = numAttr - 1; j > i; j--)
//         {
//             x[i] = x[i] - U[i][j] * x[j];
//         }
//         x[i] /= U[i][i];
//     }
// }

void ecall_dispatcher::close()
{
    TRACE_ENCLAVE("Enclave: release context and HTparams!");
    TRACE_ENCLAVE("ecall_dispatcher::close");
}