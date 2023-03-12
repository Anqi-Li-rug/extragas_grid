/*
  array3d : Utilities
*/

#ifndef TNT_ARRAY3D_UTILS_H
#define TNT_ARRAY3D_UTILS_H

#include <cstdlib>
#include <cassert>

namespace ffarray {

  template <class T>
    std::ostream& operator<<(std::ostream &s, const Array3D<T> &A)
    {
      int M=A.dim1();
      int N=A.dim2();
      int L=A.dim3();
      
      s << M << " " << N << " " << L << "\n";
      
      for (int i=0; i<M; i++)
	{
	  for (int j=0; j<N; j++)
	    {
	      for (int k=0; k<L; k++)
            	s << A[i][j][k] << " ";
			s << "\n";
	    }
	  s << "\n";
	}
      
      
      return s;
    }
  
  template <class T>
    std::istream& operator>>(std::istream &s, Array3D<T> &A)
    {
      
      int M, N, L;
      
      s >> M >> N >> L;
      
      Array3D<T> B(M,N,L);
      
      for (int i=0; i<M; i++)
        for (int j=0; j<N; j++)
	  for (int k=0; k<L; k++)
	    s >>  B[i][j][k];
      
      A = B;
      return s;
    }
  
  template <class T>
    void print(const Array3D<T> &A)
    {
      for (int k=0; k<A.dim3(); k++) {
	for (int i=0; i<A.dim1(); i++) {
	  for (int j=0; j<A.dim2(); j++) {
	    cout << A[A.dim1()-i-1][j][k] << " ";
	  }
	  cout << endl;
	}
	cout << endl;
      }
    }
  
  
  template <class T>
    Array3D<T> operator+(const Array3D<T> &A, const Array3D<T> &B)
    {
      int m = A.dim1();
      int n = A.dim2();
      int p = A.dim3();
      
      if (B.dim1() != m ||  B.dim2() != n || B.dim3() != p )
	return Array3D<T>();
      
      else
	{
	  Array3D<T> C(m,n,p);
	  
	  for (int i=0; i<m; i++)
	    for (int j=0; j<n; j++)
	      for (int k=0; k<p; k++)
		C[i][j][k] = A[i][j][k] + B[i][j][k];
	  
	  return C;
	}
    }
  
  
  template <class T>
    Array3D<T> operator-(const Array3D<T> &A, const Array3D<T> &B)
    {
      int m = A.dim1();
      int n = A.dim2();
      int p = A.dim3();
      
      if (B.dim1() != m ||  B.dim2() != n || B.dim3() != p )
	return Array3D<T>();
      
      else
	{
	  Array3D<T> C(m,n,p);
	  
	  for (int i=0; i<m; i++)
	    for (int j=0; j<n; j++)
	      for (int k=0; k<p; k++)
		C[i][j][k] = A[i][j][k] - B[i][j][k];
	  
	  return C;
	}
    }
    
  template <class T>
    Array3D<T> operator*(const Array3D<T> &A, const Array3D<T> &B)
    {
      int m = A.dim1();
      int n = A.dim2();
      int p = A.dim3();
      
      if (B.dim1() != m ||  B.dim2() != n || B.dim3() != p )
	return Array3D<T>();
      
      else
	{
	  Array3D<T> C(m,n,p);
	  
	  for (int i=0; i<m; i++)
	    for (int j=0; j<n; j++)
	      for (int k=0; k<p; k++)
		C[i][j][k] = A[i][j][k] * B[i][j][k];
	  
	  return C;
	}
    }
    
  template <class T>
    Array3D<T> operator/(const Array3D<T> &A, const Array3D<T> &B)
    {
      int m = A.dim1();
      int n = A.dim2();
      int p = A.dim3();
      
      if (B.dim1() != m ||  B.dim2() != n || B.dim3() != p )
	return Array3D<T>();
      
      else
	{
	  Array3D<T> C(m,n,p);
	  
	  for (int i=0; i<m; i++)
	    for (int j=0; j<n; j++)
	      for (int k=0; k<p; k++)
		C[i][j][k] = A[i][j][k] / B[i][j][k];
	  
	  return C;
	}
    }
  
  template <class T>
    Array3D<T>& operator+=(Array3D<T> &A, const Array3D<T> &B) 
  {
    int m = A.dim1();
    int n = A.dim2();
    int p = A.dim3();
    
    if (B.dim1() == m &&  B.dim2() == n && B.dim3() == p )
      {
	for (int i=0; i<m; i++)
	  for (int j=0; j<n; j++)
	    for (int k=0; k<p; k++)
	      A[i][j][k] += B[i][j][k];
      }
    
    return A;
  }

  template <class T>
    Array3D<T>& operator-=(Array3D<T> &A, const Array3D<T> &B)
  {
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() == m &&  B.dim2() == n && B.dim3() == p )
	{
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				for (int k=0; k<p; k++)
					A[i][j][k] -= B[i][j][k];
	}

	return A;
  }

  template <class T>
    Array3D<T>& operator*=(Array3D<T> &A, const Array3D<T> &B)
  {
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() == m &&  B.dim2() == n && B.dim3() == p )
	{
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				for (int k=0; k<p; k++)
					A[i][j][k] *= B[i][j][k];
	}

	return A;
  }


  template <class T>
  Array3D<T>& operator/=(Array3D<T> &A, const Array3D<T> &B)
  {
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() == m &&  B.dim2() == n && B.dim3() == p )
	{
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				for (int k=0; k<p; k++)
					A[i][j][k] /= B[i][j][k];
	}

	return A;
  }

  template <class T>
    T max(const Array3D<T> &A)
    {
      int m = A.dim1();
      int n = A.dim2();
      int p = A.dim3();
      T maxArray=A[0][0][0];
      for (int i=1; i<m; i++)
	for (int j=1; j<n; j++)
	  for (int k=1; k<p; k++)
	    if (A[i][j][k] > maxArray)
	      maxArray=A[i][j][k];
      return maxArray;
    }

  template <class T>
    T min(const Array3D<T> &A)
    {
      int m = A.dim1();
      int n = A.dim2();
      int p = A.dim3();
      T minArray=A[0][0][0];
      for (int i=1; i<m; i++)
	for (int j=1; j<n; j++)
	  for (int k=1; k<p; k++)
	    if (A[i][j][k] < minArray)
	      minArray=A[i][j][k];
      return minArray;
    }

  // STATISTICS
  template <class T>
    T total(const Array3D<T> &A)
    {
      T tot=0;
      for (int k=0; k<A.dim3(); k++)
	for (int j=0; j<A.dim2(); j++)
	  for (int i=0; i<A.dim1(); i++)
	    tot+=A[i][j][k];
      return tot;
    }

} // namespace array

#endif
