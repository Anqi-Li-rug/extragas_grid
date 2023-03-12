/*
 *
 * Template Numerical Toolkit (TNT)
 *
 * Mathematical and Computational Sciences Division
 * National Institute of Technology,
 * Gaithersburg, MD USA
 *
 *
 * This software was developed at the National Institute of Standards and
 * Technology (NIST) by employees of the Federal Government in the course
 * of their official duties. Pursuant to title 17 Section 105 of the
 * United States Code, this software is not subject to copyright protection
 * and is in the public domain. NIST assumes no responsibility whatsoever for
 * its use by other parties, and makes no guarantees, expressed or implied,
 * about its quality, reliability, or any other characteristic.
 *
 * May 2006: Modified by F. Fraternali
 */


#ifndef TNT_ARRAY2D_UTILS_H
#define TNT_ARRAY2D_UTILS_H

#include <cstdlib>
#include <cassert>

namespace ffarray {
	
	static int prev_i0=-1, prev_j0=-1;
	
	/*
	 INPUT / OUTPUT
	 */
	
	template <class T>
    std::ostream& operator<<(std::ostream &s, const Array2D<T> &A) {
		int M=A.dim1();
		int N=A.dim2();
		s << M << " " << N << "\n";
		for (int i=0; i<M; i++) {
			for (int j=0; j<N; j++) {
				s << A[i][j] << " ";
			}
			s << "\n";
		}
		return s;
	}
	
	template <class T>
    std::istream& operator>>(std::istream &s, Array2D<T> &A) {
		int M, N;
		s >> M >> N;
		Array2D<T> B(M,N);
		for (int i=0; i<M; i++)
			for (int j=0; j<N; j++) {
				s >>  B[i][j];
			}
		A = B;
		return s;
	}
	
	template <class T>
    void print(const Array2D<T> &A)
    // usage: print(A);
    {
		for (int i=0; i<A.dim1(); i++) {
			for (int j=0; j<A.dim2(); j++) {
				cout << A[A.dim1()-i-1][j] << " ";
			}
			cout << endl;
		}
    }
	
	template <class T>
    void print_nonzero(const Array2D<T> &A)
    {
		for (int i=0; i<A.dim1(); i++) 
			for (int j=0; j<A.dim2(); j++) 
				if (A[A.dim1()-i-1][j] != 0) 
					cout << i << " " << j << " " << A[A.dim1()-i-1][j] << endl;
    }
	
	/* 
     SUBARRAY
     
     Create a new view to a subarray defined by the boundaries
     [i0][i0] and [i1][j1].  The size of the subarray is
     (i1-i0) by (j1-j0).  If either of these lengths are zero
     or negative, the subarray view is null.
	 */
	
	template <class T>
    Array2D<T> subarray(const Array2D<T> &A, int i0, int i1, int j0, int j1) 
    {
		Array2D<T> subA;
		int m = i1-i0+1;
		int n = j1-j0+1;  
		if (m<1 || n<1)
			return subA;
		
		//      cout << m << " " << n << endl;
		//  subA.data_ = data_;
		subA.m_ = m;
		subA.n_ = n;
		subA.v_ = Array1D<T*>(m);
		T* p = &(A[0]) + i0 *  n + j0;
		for (int j=0; j<m; j++) {
			subA.v_[j] = p + j*n;
		}	
		return subA;
    }
	
	/*
	 SUBARRAY (2)
	 
	 Arguments:
	 Initial array,
	 initial vectors (x_i[i] and y_j[j] corresponds to A[i][j]),
	 coordinates of the centre of the subarray (with respect to vectors),
	 sizes of the subarray.
	 */
	template<class T> 
    Array2D<T> subarray(const Array2D<T> &A,
						Array1D<T> x_i, 
						Array1D<T> y_j, 
						T x, T y, 
						int sub_n, int sub_m)
    {
		unsigned int i0, j0, i, j;
		bool goon=true;
		if (x < x_i[0])
			if (y < y_j[0]) {
				// x and y out of grid
				i0=0; j0=0;
			}
			else {
				// only x out of grid
				i0=0, j=1;
				while ((goon) && j<y_j.dim()) {
					//	  for (int j=1; j<y_j.dim(); j++) {
					if (y_j[j-1] <= y && y_j[j] > y) {
						j0=j-sub_n/2; 
						if (j0 < 0) j0=0;
						goon=false;
					}
					j++;
				}
			}
			else {
				// only y out of grid
				if (y < y_j[0]) {
					j0=0, i=1;
					while ((goon) && i<x_i.dim()) {
						//	  for (int i=1; i<x_i.dim(); i++) {
						if (x_i[i-1] < x && x_i[i] >= x) {
							i0=i-sub_n/2; 
							if (i0 < 0) i0=0;
							goon=false;
						}
						i++;
					}
				}
				else {
					i=1;
					while ((goon) && i<x_i.dim()) {
						//for (int i=1; i<x_i.dim(); i++) {
						if (x_i[i-1] < x && x_i[i] >= x) {
							j=1;
							while ((goon) && j<y_j.dim()) {
								//		for (int j=1; j<y_j.dim() ;j++) {
								if (y_j[j-1] <= y && y_j[j] > y) {
									i0=i-sub_m/2;
									j0=j-sub_n/2;
									if (i0 < 0)
										i0=0;
									if (j0 < 0)
										j0=0;
									goon=false;
									prev_i0=i0;
									prev_j0=j0;
								}
								j++;
							}
						}
						i++;
					}
				}
			}
		// This seems necessary in caso i0 or j0 become larger than the size of the original matrix
		if (j0>A.dim2()) j0=A.dim2()-1;
		if (i0>A.dim1()) i0=A.dim1()-1;
		
		Array2D<T> subA(sub_m,sub_n);
		if (sub_m<1 || sub_n<1)
			return A;
		
		
		for (int i=0; i<sub_n; i++) {
			for (int j=0; j<sub_m; j++) {
//				cout << i << " " << j << " " << A.dim1() << " " << A.dim2() << " " << i0<< " " << j0 <<endl;
				subA[i][j]=A[i+i0][j+j0];
				//*(A[0]+A.dim1()*(j0+j+1)+i0+i+1);
			}
		}
		return subA;
    }
	
	
	/*
	 getValue
	 
	 Interpolates a matrix IN the POSITION (x,y) and return the value
	 */
	template<class T> 
    T getValue(Array2D<T> &A,
			   Array1D<T> x_i, 
			   Array1D<T> y_j, 
			   T x, T y,
			   int interp_degree)
    {
		T val;
		
		int sub_m, sub_n;
		if (interp_degree == 1) {sub_n=2, sub_m=2;}
		
		//cout << x_i << " " << y_j << endl;
		Array2D<T> subA=subarray(A, x_i, y_j, x ,y, sub_n, sub_m);
		Array1D<T> subx_i=subarray(x_i, x, sub_n);
		Array1D<T> suby_j=subarray(y_j, y, sub_m);
		
		// BI-LINEAR INTERPOLATION
		/*      cout << "x and y" << endl;
		 print(subx_i);
		 print(suby_j);
		 print(subA);
		 */
		
		/* TEMPORARY !!!!! */
		//	val=bilin(subx_i, suby_j, subA[0], x, y);
		if (interp_degree == 1) {
			float tt=(x-subx_i[0])/(subx_i[1]-subx_i[0]);
			float tu=(y-suby_j[0])/(suby_j[1]-suby_j[0]);
			val=(1-tt)*(1-tu)*subA[0][0]+tt*(1-tu)*subA[1][0]+tt*tu*subA[1][1]+(1-tt)*tu*subA[0][1];
		}
		else {
			if (interp_degree > 1) {
				// HIGHER ORDER INTERPOLATION (STILL PROBLEMS!!!
				//      NR::polin2(NRsubx_i, suby_j, subm_ij, 
				//		 x, y, &val, &val_err);
				//      T *p;
				//      splie2(p_subx_i, p_suby_j, p_subm_ij, 
				//	     interp_degree+1, interp_degree+1, 
				//	     p);
			}
			else {
				cout << "Error: degree of interpolation must be higher than 1" << endl;
			}
		}
		return val;
    }
	
	template<class T> 
    T getValueSubarray(Array2D<T> &subA,
					   Array1D<T> subx_i, 
					   Array1D<T> suby_j, 
					   T x, T y)
    {
		/* 
		 This version of getValue interpolates a matrix already reduced
		 to a SUBmatrix IN the POSITION (x,y) and returns the value
		 WITHOUT EXTRACTING A SUBMATRIX!
		 */
		T val, val_err;
		
		/* TEMPORARY !!!!! */
		//	val=bilin(subx_i, suby_j, subA[0], x, y);
		if (subx_i.dim() == 2 && suby_j.dim() == 2) {
			float tt=(x-subx_i[0])/(subx_i[1]-subx_i[0]);
			float tu=(y-suby_j[0])/(suby_j[1]-suby_j[0]);
			val=(1-tt)*(1-tu)*subA[0][0]+tt*(1-tu)*subA[1][0]+tt*tu*subA[1][1]+(1-tt)*tu*subA[0][1];
		}
		else {
			// HIGHER ORDER INTERPOLATION (STILL PROBLEMS!!!
			//      NR::polin2(NRsubx_i, suby_j, subm_ij, 
			//		 x, y, &val, &val_err);
			//      T *p;
			//      splie2(p_subx_i, p_suby_j, p_subm_ij, 
			//	     interp_degree+1, interp_degree+1, 
			//	     p);
		}
		return val;
    }
	
	template<class T> 
    Array2D<T> subarrayMem(const Array2D<T> &A,
						   Array1D<T> x_i, 
						   Array1D<T> y_j, 
						   T x, T y, 
						   int sub_n, int sub_m)
    /*
	 this version of subarray keeps in memory the last subarray used
	 and avoid to recalculate each time the subarray is it has not changed.
	 
	 useful in case of repeated operations with values of x and y that
	 do not change much
	 */
    {
		bool goon=true, flipping_i=true, flipping_j=true;
		int incr_i=1, incr_j=1;
		unsigned int i=1, i0, j=1, j0;
		if (prev_i0 > 0 && prev_i0 < x_i.dim()-sub_m)
			i=prev_i0;
		else 
			i=1;
		if (prev_j0 > 0 && prev_j0 < y_j.dim()-sub_n) 
			j=prev_j0;
		else 
			j=1;
		
		
		if (x < x_i[0])
			if (y < y_j[0]) {
				// x and y out of grid
				i0=0; j0=0;
			}
			else {
				// only x out of grid
				i0=0, j=1;
				while (goon) {
					if (y_j[j-1] <= y && y_j[j] > y) {
						j0=j-sub_n/2; 
						goon=false;
						prev_j0=j0;
					}
					j+=incr_j;
					if (flipping_j) {
						if (j>0 && j < y_j.dim())
							incr_j=-(abs(incr_j)+1)*(incr_j/abs(incr_j));
						else {
							if (j==0){
								j+=-(abs(incr_j)+1)*(incr_j/abs(incr_j));
								incr_j=1;
								flipping_j=false;
							}
							else {
								if (j==y_j.dim()){
									j+=-(abs(incr_j)+1)*(incr_j/abs(incr_j));
									incr_j=-1;
									flipping_j=false;
								}
							}
						}
					}
				}
			}
			else {
				// only y out of grid
				if (y < y_j[0]) {
					j0=0, i=1;
					while (goon) {
						if (x_i[i-1] <= x && x_i[i] > x) {
							i0=i-sub_n/2;
							if (i0 < 0)
								i0=0;
							goon=false;
							prev_i0=i0;
						}
						i+=incr_i;
						if (flipping_i) {
							if (i>0 && i < x_i.dim())
								incr_i=-(abs(incr_i)+1)*(incr_i/abs(incr_i));
							else {
								if (i==0){
									i+=-(abs(incr_i)+1)*(incr_i/abs(incr_i));
									incr_i=1;
									flipping_i=false;
								}
								else {
									if (i==x_i.dim()){
										i+=-(abs(incr_i)+1)*(incr_i/abs(incr_i));
										incr_i=-1;
										flipping_i=false;
									}
								}
							}
						}
					}
				}
				else {
					while (goon) {
						//	    cout << i << " " << prev_ii << " " << j << " " << prev_jj << endl;
						if (x_i[i-1] <= x && x_i[i] > x) {
							while (goon) {
								if (y_j[j-1] <= y && y_j[j] > y) {
									j0=j-sub_n/2;
									i0=i-sub_m/2;
									if (i0 < 0)
										i0=0;
									if (j0 < 0)
										j0=0;
									goon=false;
									prev_i0=i0;
									prev_j0=j0;
								}
								j+=incr_j;
								if (flipping_j) {
									if (j>0 && j < y_j.dim())
										incr_j=-(abs(incr_j)+1)*(incr_j/abs(incr_j));
									else {
										if (j==0){
											j+=-(abs(incr_j)+1)*(incr_j/abs(incr_j));
											incr_j=1;
											flipping_j=false;
										}
										else {
											if (j==y_j.dim()){
												j+=-(abs(incr_j)+1)*(incr_j/abs(incr_j));
												incr_j=-1;
												flipping_j=false;
											}
										}
									}
								}
							}
						}
						i+=incr_i;
						if (flipping_i) {
							if (i>0 && i < x_i.dim())
								incr_i=-(abs(incr_i)+1)*(incr_i/abs(incr_i));
							else {
								if (i==0){
									i+=-(abs(incr_i)+1)*(incr_i/abs(incr_i));
									incr_i=1;
									flipping_i=false;
								}
								else {
									if (i==x_i.dim()){
										i+=-(abs(incr_i)+1)*(incr_i/abs(incr_i));
										incr_i=-1;
										flipping_i=false;
									}
								}
							}
						}
					}
				}
			}
		Array2D<T> subA(sub_m,sub_n);
		if (sub_m<1 || sub_n<1)
			return A;
		
		for (int i=0; i<sub_n; i++) {
			for (int j=0; j<sub_m; j++) {
				subA[i][j]=A[i+i0][j+j0];
				//*(A[0]+A.dim1()*(j0+j+1)+i0+i+1);
			}
		}
		return subA;
    }
	
	
	/*
	 MATHEMATICAL OPERATORS
	 */
	
	template <class T>
    Array2D<T> operator+(const Array2D<T> &A, const Array2D<T> &B) {
		int m = A.dim1();
		int n = A.dim2();
		if (B.dim1() != m ||  B.dim2() != n )
			return Array2D<T>();
		else {
			Array2D<T> C(m,n);
			for (int i=0; i<m; i++) {
				for (int j=0; j<n; j++)
					C[i][j] = A[i][j] + B[i][j];
			}
			return C;
		}
	}
	
	template <class T>
    Array2D<T> operator-(const Array2D<T> &A, const Array2D<T> &B) {
		int m = A.dim1();
		int n = A.dim2();
		if (B.dim1() != m ||  B.dim2() != n )
			return Array2D<T>();
		else {
			Array2D<T> C(m,n);
			for (int i=0; i<m; i++) {
				for (int j=0; j<n; j++)
					C[i][j] = A[i][j] - B[i][j];
			}
			return C;
		}
	}
	
	template <class T>
    Array2D<T> operator*(const Array2D<T> &A, const Array2D<T> &B) {
		int m = A.dim1();
		int n = A.dim2();
		if (B.dim1() != m ||  B.dim2() != n )
			return Array2D<T>();
		else
		{
			Array2D<T> C(m,n);
			for (int i=0; i<m; i++)
			{
				for (int j=0; j<n; j++)
					C[i][j] = A[i][j] * B[i][j];
			}
			return C;
		}
	}
	
	template <class T>
    Array2D<T> operator/(const Array2D<T> &A, const Array2D<T> &B)
    {
		int m = A.dim1();
		int n = A.dim2();
		if (B.dim1() != m ||  B.dim2() != n )
			return Array2D<T>();
		else
		{
			Array2D<T> C(m,n);
			for (int i=0; i<m; i++)
			{
				for (int j=0; j<n; j++)
					C[i][j] = A[i][j] / B[i][j];
			}
			return C;
		}
    }
	
	template <class T>
    Array2D<T> operator/(const Array2D<T> &A, const T b)
    {
		int m = A.dim1();
		int n = A.dim2();
		Array2D<T> C(m,n);
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				C[i][j] = A[i][j] / b;
		return C;
    }
	
	template <class T>
    Array2D<T>&  operator+=(Array2D<T> &A, const Array2D<T> &B)
	{
		int m = A.dim1();
		int n = A.dim2();
		if (B.dim1() == m ||  B.dim2() == n )
		{
			for (int i=0; i<m; i++)
			{
				for (int j=0; j<n; j++)
					A[i][j] += B[i][j];
			}
		}
		return A;
	}
	
	template <class T>
    Array2D<T>&  operator-=(Array2D<T> &A, const Array2D<T> &B)
	{
		int m = A.dim1();
		int n = A.dim2();
		if (B.dim1() == m ||  B.dim2() == n )
		{
			for (int i=0; i<m; i++)
			{
				for (int j=0; j<n; j++)
					A[i][j] -= B[i][j];
			}
		}
		return A;
	}
	
	template <class T>
    Array2D<T>&  operator*=(Array2D<T> &A, const Array2D<T> &B)
	{
		int m = A.dim1();
		int n = A.dim2();
		if (B.dim1() == m ||  B.dim2() == n )
		{
			for (int i=0; i<m; i++)
			{
				for (int j=0; j<n; j++)
					A[i][j] *= B[i][j];
			}
		}
		return A;
	}
	
	template <class T>
    Array2D<T>&  operator/=(Array2D<T> &A, const Array2D<T> &B)
	{
		int m = A.dim1();
		int n = A.dim2();
		if (B.dim1() == m ||  B.dim2() == n )
		{
			for (int i=0; i<m; i++)
			{
				for (int j=0; j<n; j++)
					A[i][j] /= B[i][j];
			}
		}
		return A;
	}
	
	/**
     Matrix Multiply:  compute C = A*B, where C[i][j]
     is the dot-product of row i of A and column j of B.
     
     @param A an (m x n) array
     @param B an (n x k) array
     @return the (m x k) array A*B, or a null array (0x0)
     if the matrices are non-conformant (i.e. the number
     of columns of A are different than the number of rows of B.)
	 */
	template <class T>
    Array2D<T> matmult(const Array2D<T> &A, const Array2D<T> &B)
    {
		if (A.dim2() != B.dim1())
			return Array2D<T>();
		int M = A.dim1();
		int N = A.dim2();
		int L = B.dim2();
		Array2D<T> C(M,L);
		for (int i=0; i<M; i++)
			for (int j=0; j<L; j++)
			{
				T sum = 0;
				for (int k=0; k<N; k++)
					sum += A[i][k] * B [k][j];
				C[i][j] = sum;
			}
		return C;
    }
	
	// STATISTICS
	template <class T>
    T total(const Array2D<T> &A)
    {
		T tot=0;
		for (int j=0; j<A.dim2(); j++)
			for (int i=0; i<A.dim1(); i++)
				tot+=A[i][j];
		return tot;
    }
	
	
	
	template <class T>
    Array2D<T> smooth2D(const Array2D<T> & A, float ratio)
    {
		int x_sizeKernel, y_sizeKernel;
		/*      Array2D<double> kernel.makeKernel(1., 
		 1., 
		 ratio, 
		 ratio,
		 x_sizeKernel,
		 y_sizeKernel);
		 
		 int m=A.dim1(), n=A.dim2();
		 
		 Array2D<T> smooA(m, n);
		 
		 for (int j=y_sizeKernel/2;j<n-y_sizeKernel/2;j++)
		 for (int i=x_sizeKernel/2;i<m-x_sizeKernel/2;i++)
		 smooA[i][j]=0.;
		 
		 for (int j=y_sizeKernel/2;j<n-y_sizeKernel/2;j++)
		 for (int i=x_sizeKernel/2;i<m-x_sizeKernel/2;i++)
		 for (int kj=0;kj<y_sizeKernel;kj++)
		 for (int ki=0;ki<x_sizeKernel;ki++){
		 smooA[i][j]+=kernel[ki][kj]*A[i-x0+ki][j-y0+kj];
		 }
		 return smooA;
		 */
    }
	
	// STATISTICS
	
	template <class T>
    T max(const Array2D<T> &A)
    {
		int m = A.dim1();
		int n = A.dim2();
		T maxArray=A[0][0];
		for (int i=1; i<m; i++)
			for (int j=1; j<n; j++)
				if (A[i][j] > maxArray)
					maxArray=A[i][j];
		return maxArray;
    }
	
	template <class T>
    T min(const Array2D<T> &A)
    {
		int m = A.dim1();
		int n = A.dim2();
		T minArray=A[0][0];
		for (int i=1; i<m; i++)
			for (int j=1; j<n; j++)
				if (A[i][j] < minArray)
					minArray=A[i][j];
		return minArray;
    }
	
	template <class T>
    Array2D<T> histogram(const Array1D<T> &A, T min, T max, T step)
    {
		int nbins=(int)((max-min)/step+0.5)+1, j;
		Array2D<T> histo(nbins,2);
		bool found;
		
		for (int i=0; i<A.dim1(); i++){
			if (A[i] > min && A[i] < max){
				found=false;
				j=0;
				while (found == false) {
					if (A[i] > min && A[i] <= min+step*(j+1)){
						histo[j][1]=histo[j][1]+1;
						found=true;
					}
					j++;
				}
			}
		}
		for (int k=0; k<nbins; k++)
			histo[k][0]=min+k*step+step/2.;
		return histo;
    }
	
	template <class T>
    Array2D<T> histogram(const Array1D<T> &A, int nbins)
    {
		int j;
		Array2D<T> histo(nbins,2);
		for (int i=0; i<nbins; i++)
			histo[i][1]=0; // this is necessary to inizialize histo (in case it has been used before by another call)
		
		T Amin=min(A);
		T Amax=max(A);
		T step=(Amax-Amin)/(nbins-3);
		Amin=Amin-3*step/2.;
		Amax=Amax+3*step/2.;
		bool found;
		
		for (int i=0; i<A.dim1(); i++){
			if (A[i] > Amin && A[i] < Amax){
				found=false;
				j=0;
				while (found == false) {
					if (A[i] > Amin && A[i] <= Amin+step*(j+1)){
						histo[j][1]=histo[j][1]+1;
						found=true;
					}
					j++;
				}
			}
		}
		for (int k=0; k<nbins; k++)
			histo[k][0]=Amin+k*step+step/2.;
		return histo;
    }
	
	template <class T>
    Array2D<T> histogram_no0(const Array1D<T> &A, int nbins)
    {
		int j=0, dimno0=0;
		Array2D<T> histo(nbins,2);
		for (int i=0; i<nbins; i++)
			histo[i][1]=0; // this is necessary to inizialize histo (in case it has been used before by another call)
		
		for (int i=0; i<A.dim1(); i++)
			if (A[i] != 0)
				dimno0++;
		
		Array1D<T> Ano0(dimno0);
		
		for (int i=0; i<A.dim1(); i++)
			if (A[i] != 0) {
				Ano0[j]=A[i];
				j++;
			}
		
		T Amin=min(Ano0);
		T Amax=max(Ano0);
		T step=(Amax-Amin)/(nbins-3);
		Amin=Amin-3*step/2.;
		Amax=Amax+3*step/2.;
		bool found;
		
		for (int i=0; i<dimno0; i++){
			if (Ano0[i] > Amin && Ano0[i] < Amax){
				found=false;
				j=0;
				while (found == false) {
					if (Ano0[i] > Amin && Ano0[i] <= Amin+step*(j+1)){
						histo[j][1]=histo[j][1]+1;
						found=true;
					}
					j++;
				}
			}
		}
		for (int k=0; k<nbins; k++)
			histo[k][0]=Amin+k*step+step/2.;
		return histo;
    }
	
	
	template <class T>
    Array2D<T> histogram_no0(const Array1D<T> &A, T Amin, T Amax, T step)
    {
		int nbins=(int)((Amax-Amin)/step+0.5)+1;
		
		int j=0, dimno0=0;
		Array2D<T> histo(nbins,2);
		for (int i=0; i<nbins; i++)
			histo[i][1]=0; // this is necessary to inizialize histo (in case it has been used before by another call)
		
		for (int i=0; i<A.dim1(); i++)
			if (A[i] != 0)
				dimno0++;
		
		Array1D<T> Ano0(dimno0);
		
		for (int i=0; i<A.dim1(); i++)
			if (A[i] != 0) {
				Ano0[j]=A[i];
				j++;
			}
		
		bool found;
		
		for (int i=0; i<dimno0; i++){
			if (Ano0[i] > Amin && Ano0[i] < Amax){
				found=false;
				j=0;
				while (found == false) {
					if (Ano0[i] > Amin && Ano0[i] <= Amin+step*(j+1)){
						histo[j][1]=histo[j][1]+1;
						found=true;
					}
					j++;
				}
			}
		}
		for (int k=0; k<nbins; k++)
			histo[k][0]=Amin+k*step+step/2.;
		return histo;
    }
	
}

#endif

