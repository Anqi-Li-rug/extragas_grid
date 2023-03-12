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
 */

#ifndef TNT_ARRAY1D_UTILS_H
#define TNT_ARRAY1D_UTILS_H

#include <cstdlib>
#include <cassert>

namespace ffarray {
  
  template <class T>
    std::ostream& operator<<(std::ostream &s, const Array1D<T> &A)
    {
      int N=A.dim1();
#ifdef TNT_DEBUG
      s << "addr: " << (void *) &A[0] << "\n";
#endif
      s << N << "\n";
      for (int j=0; j<N; j++) {
	s << A[j] << "\n";
      }
      s << "\n";
      return s;
    }
  
  template <class T>
    std::istream& operator>>(std::istream &s, Array1D<T> &A)
    {
      int N;
      s >> N;
      Array1D<T> B(N);
      for (int i=0; i<N; i++)
	s >> B[i];
      A = B;
      return s;
    }
  
  template <class T>
    void print(const Array1D<T> &v)
    {
      for (int i=0; i<v.dim1(); i++) 
	cout << v[i] << " ";
      cout << endl;
    }
  
  /* 
     SUBARRAYS
  */
  static int prev_i;

  template<class T> 
    Array1D<T> subarrayMem(Array1D<T> v, double x, int sub_n)
    { 
      bool goon=true, flipping=true;
      int incr_i=1;
      unsigned int i=1, i0;
     if (prev_i > 0 && prev_i < v.dim()) 
	i=prev_i;
      else 
	i=1;
      if (x < v[0])
	i0=0;
      else {
	//	cout << i << " " << incr_i << " " << prev_i << endl;
	while (goon) {
	  //	cout << x << " " << v[i-1] << " " << v[i] << endl;
	  if (v[i-1] <= x && v[i] > x) {
	    i0=i-sub_n/2;
	    if (i0 < 0)
	      i0=0;
	    goon=false;
	    prev_i=i;
	  }
	  i+=incr_i;
	  if (flipping) {
	    if (i>0 && i < v.dim())
	      incr_i=-(abs(incr_i)+1)*(incr_i/abs(incr_i));
	    else {
	      if (i==0){
		i+=-(abs(incr_i)+1)*(incr_i/abs(incr_i));
		incr_i=1;
		flipping=false;
	      }
	      else {
		if (i==v.dim()){
		  i+=-(abs(incr_i)+1)*(incr_i/abs(incr_i));
		  incr_i=-1;
		  flipping=false;
		}
	      }
	    }
	  }
	}
      }
      Array1D<T> subv(sub_n);
      if (sub_n<1)
	return v;
      
      for (int i=0; i<sub_n; i++) {
	subv[i]=v[i0+i];
      }
      return subv;
    }

  template<class T> 
    Array1D<T> subarray(Array1D<T> v, double x, int sub_n)
    { 
      int i0=1;
      int i=1;
      bool goon=true;
      if (x < v[0])
	i0=0;
      else {
	while ((goon) && i<v.dim()) {
	  if (v[i-1] < x && v[i] >= x) {
	    i0=i-sub_n/2; 
	    goon=false;
	  }
	  i++;
	}
      }
      Array1D<T> subv(sub_n);
      if (sub_n<1)
	return v;
      
      for (int i=0; i<sub_n; i++) {
	subv[i]=v[i0+i];
      }
      return subv;
    }

  /*
    Mathematical operators
  */
  
  template <class T>
    Array1D<T> operator+(const Array1D<T> &A, const Array1D<T> &B)
    {
      int n = A.dim1();
      if (B.dim1() != n )
	return Array1D<T>();
      else {
	Array1D<T> C(n);
	for (int i=0; i<n; i++) {
	  C[i] = A[i] + B[i];
	}
	return C;
      }
    }
  
  template <class T>
    Array1D<T> operator-(const Array1D<T> &A, const Array1D<T> &B)
    {
      int n = A.dim1();
      if (B.dim1() != n )
	return Array1D<T>();
      else {
	Array1D<T> C(n);
	for (int i=0; i<n; i++) {
	  C[i] = A[i] - B[i];
	}
	return C;
      }
    }
  
  template <class T>
    Array1D<T> operator*(const Array1D<T> &A, const Array1D<T> &B)
    {
      int n = A.dim1();
      if (B.dim1() != n )
	return Array1D<T>();
      else {
	Array1D<T> C(n);
	for (int i=0; i<n; i++) {
	  C[i] = A[i] * B[i];
	}
	return C;
      }
    }

  template <class T>
    Array1D<T> operator*(const Array1D<T> &A, const float b)
    {
      int n = A.dim1();
      Array1D<T> C(n);
      for (int i=0; i<n; i++) {
	C[i] = A[i] * b;
      }
      return C;
    }
  
  template <class T>
    Array1D<T> operator/(const Array1D<T> &A, const Array1D<T> &B)
    {
      int n = A.dim1();
      if (B.dim1() != n )
	return Array1D<T>();
      else {
	Array1D<T> C(n);
	for (int i=0; i<n; i++)
	  {
	    C[i] = A[i] / B[i];
	  }
	return C;
      }
    }
  
  template <class T>
    Array1D<T> operator/(const Array1D<T> &A, const float b)
    {
      int n = A.dim1();
      Array1D<T> C(n);
      for (int i=0; i<n; i++) {
	C[i] = A[i] / b;
      }
      return C;
    }

  template <class T>
    Array1D<T>&  operator+=(Array1D<T> &A, const Array1D<T> &B)
  {
    int n = A.dim1();
    if (B.dim1() == n) {
      for (int i=0; i<n; i++)
	{
	  A[i] += B[i];
	}
    }
    return A;
  }
  
  template <class T>
    Array1D<T>&  operator-=(Array1D<T> &A, const Array1D<T> &B)
  {
    int n = A.dim1();
    if (B.dim1() == n) {
      for (int i=0; i<n; i++)
	{
	  A[i] -= B[i];
	}
    }
    return A;
  }
  
  template <class T>
    Array1D<T>&  operator*=(Array1D<T> &A, const Array1D<T> &B)
  {
    int n = A.dim1();
    if (B.dim1() == n) {
      for (int i=0; i<n; i++) {
	A[i] *= B[i];
      }
    }
    return A;
  }

  template <class T>
    Array1D<T>&  operator*=(Array1D<T> &A, const float b)
  {
    int n = A.dim1();
    Array1D<T> C(n);
    for (int i=0; i<n; i++) {
      C[i] = A[i] * b;
    }
    return C;
  }
  
  template <class T>
    Array1D<T>&  operator/=(Array1D<T> &A, const Array1D<T> &B)
  {
    int n = A.dim1();
    if (B.dim1() == n) {
      for (int i=0; i<n; i++) {
	A[i] /= B[i];
      }
    }
    return A;
  }

  // STATISTICS
  
  template <class T>
    T max(const Array1D<T> &A)
    {
      int n = A.dim1();
      T maxArray=A[0];
      for (int i=1; i<n; i++)
	if (A[i] > maxArray)
	  maxArray=A[i];
      return maxArray;
    }

  template <class T>
    T min(const Array1D<T> &A)
    {
      int n = A.dim1();
      T minArray=A[0];
      for (int i=1; i<n; i++)
	if (A[i] < minArray)
	  minArray=A[i];
      return minArray;
    }

  template <class T>
    T total(const Array1D<T> &A)
    {
      T tot=0;
      for (int i=0; i<A.dim1(); i++)
          tot+=A[i];
      return tot;
    }

  template <class T>
    T mean(const Array1D<T> &A)
    {
      T mean=0;
      for (int i=0; i<A.dim1(); i++)
          mean+=A[i];
      return mean/A.dim1();
    }

  template <class T>
    T mean_no0(const Array1D<T> &A)
    {
      T mean=0;
      int n=0;
      for (int i=0; i<A.dim1(); i++) {
	if (A[i] != 0) {
          mean+=A[i];
	  n++;
	}
      }
      //      print(A);
      //      cout << mean << " " << n << endl;
      return mean/n;
    }

  template <class T>
    T variance(const Array1D<T> &A)
    {
      T variance=0;
      T ave=mean(A);
      for (int i=0; i<A.dim1(); i++) {
	variance+=(A[i]-ave)*(A[i]-ave);
      }
      return variance/A.dim1();
    }

  template <class T>
    T variance_no0(const Array1D<T> &A)
    {
      T variance=0;
      T ave=mean_no0(A);
      int n=0;
      for (int i=0; i<A.dim1(); i++) {
	if (A[i] != 0) {
	  variance+=fabs((A[i]-ave)*(A[i]-ave));
	  n++;
	}
      }
      return variance/n;
    }


}

#endif
  
  
