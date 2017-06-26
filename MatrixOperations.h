//CMF_Project5
// LU decomposition,Matrix inversion, multiplication, etc.

#include <vector>
using namespace std;

// LU decomposition

// Structure to store upper and lower matrices
typedef struct {
	vector<vector<double> > L;
	vector<vector<double> > U;
}LUstruct;

LUstruct LU_Decomposition(vector<vector<double> > A) {
	int N = A.size();
	vector<vector<double> >  B(N, vector<double>(N));
	for (int i = 0; i <= N - 1; i++)
		for (int j = 0; j <= N - 1; j++)
			B[i][j] = A[i][j];

	for (int k = 0; k <= N - 2; k++) {
		for (int i = k + 1; i <= N - 1; i++)
			B[i][k] = B[i][k] / B[k][k];
		for (int j = k + 1; j <= N - 1; j++)
			for (int i = k + 1; i <= N - 1; i++)
				B[i][j] = B[i][j] - B[i][k] * B[k][j];
	}
	vector<vector<double> >  L(N, vector<double>(N));
	vector<vector<double> >  U(N, vector<double>(N));
	for (int i = 0; i <= N - 1; i++) {
		L[i][i] = 1.0;
		for (int j = 0; j <= N - 1; j++) {
			if (i>j)
				L[i][j] = B[i][j];
			else
				U[i][j] = B[i][j];
		}
	}
	LUstruct LU_Mat;
	LU_Mat.L = L;
	LU_Mat.U = U;
	return LU_Mat;
}

// Multiply two matrices: (m x k) %*% (k x n)=(m x n)
vector<vector<double> > Mult_Mat(vector<vector<double> >A, vector<vector<double> >B, int m, int k, int n) {
	vector<vector<double> >  C(m, vector<double>(n));
	for (int j = 0; j <= n - 1; j++)
		for (int i = 0; i <= m - 1; i++) {
			C[i][j] = 0;
			for (int l = 0; l <= k - 1; l++)
				C[i][j] += A[i][l] * B[l][j];
		}
	return C;
}



// Inverse of a matrix through LU decomposition
vector<vector<double> >  Inv_LU(vector<vector<double> > A) {
	LUstruct LU_Mat;
	LU_Mat = LU_Decomposition(A);
	vector<vector<double> >  L = LU_Mat.L;
	vector<vector<double> >  U = LU_Mat.U;

	// Inverse of the upper triangular matrix
	int N_U = U.size();
	vector<vector<double> >  Uinv(N_U, vector<double>(N_U));
	for (int j = N_U - 1; j >= 0; j--) {
		Uinv[j][j] = 1.0 / U[j][j];
		for (int i = j - 1; i >= 0; i--)
			for (int k = i + 1; k <= j; k++)
				Uinv[i][j] -= 1.0 / U[i][i] * U[i][k] * Uinv[k][j];
	}

	//  Inverse of a lower triangular matrix
		int N_L = L.size();
		vector<vector<double> >  Linv(N_L, vector<double>(N_L));
		for (int i = 0; i <= N_L - 1; i++) {
			Linv[i][i] = 1.0 / L[i][i];
			for (int j = i - 1; j >= 0; j--)
				for (int k = i - 1; k >= j; k--)
					Linv[i][j] -= 1.0 / L[i][i] * L[i][k] * Linv[k][j];
		}

	int m = Linv.size();
	int k = Linv[0].size();
	int n = Uinv[0].size();
	vector<vector<double> >  Ainv = Mult_Mat(Uinv, Linv, m, k, n);
	return Ainv;
}


// Transpose of a matrix
vector<vector<double> > Trans_Mat(vector<vector<double> >A) {
	int rows = A.size();
	int cols = A[0].size();
	vector<vector<double> >  B(cols, vector<double>(rows));
	for (int i = 0; i<cols; i++)
		for (int j = 0; j<rows; j++)
			B[i][j] = A[j][i];
	return B;
}

// Multiply a matrix by a vector
vector<double> Mult_Vec(vector<vector<double> > X, vector<double> Y) {
	int rows = X.size();
	int cols = X[0].size();
	vector<double> XY(rows);
	for (int l = 0; l<rows; l++) {
		XY[l] = 0.0;
		for (int c = 0; c<cols; c++)
			XY[l] += X[l][c] * Y[c];
	}
	return XY;
}

// Beta coefficients of Regression 
vector<double> Beta(vector<vector<double> > X, vector<double> Y) {
	int rows = X.size();
	int cols = X[0].size();
	vector<vector<double> > Xt = Trans_Mat(X);
	vector<vector<double> > XtX = Mult_Mat(Xt, X, cols, rows, cols);
	vector<vector<double> > XtXinv = Inv_LU(XtX);
	vector<double> XtY = Mult_Vec(Xt, Y);
	return Mult_Vec(XtXinv, XtY);
}