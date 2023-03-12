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
*
* N.B. I tried to separate the header from the definitions, it didn't work.
* It complains about undefined symbols. I think the problem is related to
* the presence of templates, I read this somewhere...
*
*/

#ifndef TNT_ARRAY2D_H
#define TNT_ARRAY2D_H

#include <cstdlib>
#include <cmath>
#include <iostream>
#ifdef TNT_BOUNDS_CHECK
#include <assert.h>
#endif

#include "array1d.h"

extern const double PI;//=3.14159265359;

namespace ffarray {
  
  const float bcorr=2.354820045; // FWHM = 2.35 sigma
  const float bcorr2=1.067; // for Gaussian convolution
  
  template <class T>
    class Array2D 
    {
    private:
      Array1D<T> data_;
      Array1D<T*> v_;
      int m_;
      int n_;
      
    public:
      typedef T value_type;
      Array2D();
      Array2D(int m, int n);
      Array2D(int m, int n,  T *a);
      Array2D(int m, int n, const T &a);
      inline Array2D(const Array2D &A);
      inline operator T**();
      inline operator const T**();
      inline Array2D & operator=(const T &a);
      inline Array2D & operator=(const Array2D &A);
      inline Array2D & ref(const Array2D &A);
      Array2D copy() const;
      Array2D & inject(const Array2D & A);
      inline T* operator[](int i);
      inline const T* operator[](int i) const;
      inline int dim1() const;
      inline int dim2() const;
      ~Array2D();
      
      inline int ref_count();
      inline int ref_count_data();
      inline int ref_count_dim1();
      Array2D subarray(int i0, int i1, int j0, int j1);
      //  Array2D subarray(Array1D<T> , Array1D<T>, T, T, int, int); 

      //      Array2D<double> makeKernel(float, float, float, float, int&, int& );
      Array2D<double> makeKernel(float, float, float, float);
      int x_sizeKernel, y_sizeKernel, x0, y0, sizeKernel;
      Array2D smooth2D(float, float, float, float);
      Array2D smooth2D(float, float);
      Array2D smooth2D(float);
      Array2D smooth2D_keepsize(float);

      Array2D smooth1D(char, float, float);
      Array2D smooth1D_noZero(char, float, float);
    };


  template <class T>
  Array2D<T>::Array2D() 
    : data_(), v_(), m_(0), n_(0) {} 
  
  template <class T>
  Array2D<T>::Array2D(const Array2D<T> &A) 
    : data_(A.data_), v_(A.v_), m_(A.m_), n_(A.n_) {}
  
  template <class T>
  Array2D<T>::Array2D(int m, int n) 
    : data_(m*n), v_(m), m_(m), n_(n)
  {
    if (m>0 && n>0) {
      T* p = &(data_[0]);
      for (int i=0; i<m; i++) {
	v_[i] = p;
	p += n;
      }
    }
  }

  template <class T>
  Array2D<T>::Array2D(int m, int n, const T &val) 
    : data_(m*n), v_(m), m_(m), n_(n) 
  {
    if (m>0 && n>0) {
      data_ = val;
      T* p  = &(data_[0]);
      for (int i=0; i<m; i++) {
	v_[i] = p;
	p += n;
      }
    }
  }
  
  template <class T>
  Array2D<T>::Array2D(int m, int n, T *a) 
    : data_(m*n, a), v_(m), m_(m), n_(n)
  {
    if (m>0 && n>0) {
      T* p = &(data_[0]);
      for (int i=0; i<m; i++) {
	v_[i] = p;
	p += n;
      }
    }
  }
  
  /*
    OVERLOADED OPERATORS
  */
  template <class T>
  inline T* Array2D<T>::operator[](int i) 
  { 
#ifdef TNT_BOUNDS_CHECK
    assert(i >= 0);
    assert(i < m_);
#endif
    return v_[i]; 
  }
  
  template <class T>
  inline const T* Array2D<T>::operator[](int i) const
  { 
#ifdef TNT_BOUNDS_CHECK
    assert(i >= 0);
    assert(i < m_);
#endif
    return v_[i]; 
  }
  
  template <class T>
  Array2D<T> & Array2D<T>::operator=(const T &a)
  {
    /* non-optimzied, but will work with subarrays in future verions */
    for (int i=0; i<m_; i++)
      for (int j=0; j<n_; j++)
	v_[i][j] = a;
    return *this;
  }
  
  template <class T>
  inline Array2D<T>::operator T**()
  {
    return &(v_[0]);
  }
  
  template <class T>
  inline Array2D<T>::operator const T**()
  {
    return &(v_[0]);
  }
  
  /*
    COPY and INJECT CONSTRUCTORS
  */
  template <class T>
  Array2D<T> Array2D<T>::copy() const
  {
    Array2D A(m_, n_);
    for (int i=0; i<m_; i++)
      for (int j=0; j<n_; j++)
	A[i][j] = v_[i][j];
    return A;
  }
  
  template <class T>
  Array2D<T> & Array2D<T>::inject(const Array2D &A)
  {
    if (A.m_ == m_ &&  A.n_ == n_) {
      for (int i=0; i<m_; i++)
	for (int j=0; j<n_; j++)
	  v_[i][j] = A[i][j];
    }
    return *this;
  }
  
  template <class T>
  Array2D<T> & Array2D<T>::ref(const Array2D<T> &A)
  {
    if (this != &A) {
      v_ = A.v_;
      data_ = A.data_;
      m_ = A.m_;
      n_ = A.n_;
    }
    return *this;
  }
  
  template <class T>
  Array2D<T> & Array2D<T>::operator=(const Array2D<T> &A)
  {
    return ref(A);
  }
  
  template <class T>
  inline int Array2D<T>::dim1() const { return m_; }
  
  template <class T>
  inline int Array2D<T>::dim2() const { return n_; }
  
  template <class T>
  Array2D<T>::~Array2D() {}
  
  /*
    SUBARRAY
    
    Create a new view to a subarray defined by the boundaries
    [i0][i0] and [i1][j1].  The size of the subarray is
    (i1-i0) by (j1-j0).  If either of these lengths are zero
    or negative, the subarray view is null.
  */
  template <class T>
  Array2D<T> Array2D<T>::subarray(int i0, int i1, int j0, int j1) 
  {
    Array2D<T> A;
    int m = i1-i0+1;
    int n = j1-j0+1;  
    /* if either length is zero or negative, this is an invalide
       subarray. return a null view.
    */
    if (m<1 || n<1)
      return A;
    
    cout << m << " " << n << endl;
    A.data_ = data_;
    A.m_ = m;
    A.n_ = n;
    A.v_ = Array1D<T*>(m);
    T* p = &(data_[0]) + i0 *  n_ + j0;
    for (int j=0; j<m; j++) {
      A.v_[j] = p + j*n_;
    }	
    return A;
  }
  
  template <class T>
  inline int Array2D<T>::ref_count()
  {
    return ref_count_data();
  }
  
  template <class T>
  inline int Array2D<T>::ref_count_data()
  {
    return data_.ref_count();
  }
  
  template <class T>
  inline int Array2D<T>::ref_count_dim1()
  {
    return v_.ref_count();
  }
  
  /*
    SMOOTHING ARRAYS
  */
  template <class T>
  Array2D<double> Array2D<T>::makeKernel(float old_x, 
					 float old_y, 
					 float new_x, 
					 float new_y)
    //					   int& x_sizeKernel,
    //					   int& y_sizeKernel)
  {
    // NB. IT DOES NOT WORK YET FOR ELLIPTICAL SMOOTHING: USE ROUND!!!!
    
    float sig_x=(new_x/old_x/bcorr/bcorr2);
    float sig_y=(new_y/old_y/bcorr/bcorr2);
    
    x_sizeKernel=(int)(10.*sig_x+0.5);
    y_sizeKernel=(int)(10.*sig_y+0.5);
    
    if (x_sizeKernel%2 == 0)
      x_sizeKernel++;
    if (y_sizeKernel%2 == 0)
      y_sizeKernel++;
    
    Array2D<double> kernel(x_sizeKernel, y_sizeKernel);
    
    //cout << sig_x << " " << x_sizeKernel << " " << sig_y << " " << y_sizeKernel << endl;
    
    x0=x_sizeKernel/2;
    y0=y_sizeKernel/2;
    
    for (int i=0;i<x_sizeKernel;i++) {
      for (int j=0;j<y_sizeKernel;j++) {
	kernel[i][j]=1/2./PI/sig_x/sig_y*exp(-(i-x0)*(i-x0)/2./sig_x/sig_x)*exp(-(j-y0)*(j-y0)/2./sig_y/sig_y);
      }
    }
    return kernel;
  }
  
  template <class T>
  Array2D<T> Array2D<T>::smooth2D(float old_x, 
				  float old_y, 
				  float new_x, 
				  float new_y)
  {
    Array2D<double> kernel=makeKernel(old_x, old_y, new_x, new_y);
    
    Array2D<T> smooA(m_, n_, 0.);
    
    for (int j=y_sizeKernel/2;j<n_-y_sizeKernel/2;j++)
      for (int i=x_sizeKernel/2;i<m_-x_sizeKernel/2;i++)
	for (int kj=0;kj<y_sizeKernel;kj++)
	  for (int ki=0;ki<x_sizeKernel;ki++){
	    if (v_[i-x0+ki][j-y0+kj] > 0)
	      smooA[i][j]+=kernel[ki][kj]*v_[i-x0+ki][j-y0+kj];
	  }
    return smooA;
  }
  
  template <class T>
    Array2D<T> Array2D<T>::smooth2D(float old_beam, 
				  float new_beam)
    {
      Array2D<double> kernel=makeKernel(old_beam, old_beam, new_beam, new_beam);
      
      Array2D<T> smooA(m_, n_, 0.);
      
      for (int j=y_sizeKernel/2;j<n_-y_sizeKernel/2;j++)
	for (int i=x_sizeKernel/2;i<m_-x_sizeKernel/2;i++)
	  for (int kj=0;kj<y_sizeKernel;kj++)
	    for (int ki=0;ki<x_sizeKernel;ki++){
	      smooA[i][j]+=kernel[ki][kj]*v_[i-x0+ki][j-y0+kj];
	    }
      return smooA;
    }
  
  template <class T>
    Array2D<T> Array2D<T>::smooth2D(float ratio)
    {
      Array2D<double> kernel=makeKernel(1., 1., ratio, ratio);
      
      Array2D<T> smooA(m_, n_, 0.);
    
      for (int j=y_sizeKernel/2;j<n_-y_sizeKernel/2;j++)
	for (int i=x_sizeKernel/2;i<m_-x_sizeKernel/2;i++) {
	  for (int kj=0;kj<y_sizeKernel;kj++)
	    for (int ki=0;ki<x_sizeKernel;ki++){
	      smooA[i][j]+=kernel[ki][kj]*v_[i-x0+ki][j-y0+kj];
	    }
	}
    
      return smooA;
  }
  
  template <class T>
    Array2D<T> Array2D<T>::smooth2D_keepsize(float ratio)
    {
      Array2D<double> kernel=makeKernel(1.f, 1.f, ratio, ratio);
      
      Array2D<T> smooA(m_, n_);
      
         double renormK;
   renormK=0;

      // Bottom row
      for (int j=0; j<y_sizeKernel/2; j++)
	for (int i=x_sizeKernel/2;i<m_-x_sizeKernel/2;i++) {
	  renormK=0;
	  for (int kj=y_sizeKernel/2-j; kj<y_sizeKernel; kj++)
	    for (int ki=0;ki<x_sizeKernel;ki++)
	      renormK+=kernel[ki][kj];
	  for (int kj=y_sizeKernel/2-j; kj<y_sizeKernel; kj++)
	    for (int ki=0;ki<x_sizeKernel;ki++) 
	      smooA[i][j]+=kernel[ki][kj]*v_[i-x0+ki][j-y0+kj]/renormK;
	}
      
      // Top row
      for (int j=n_-y_sizeKernel/2; j<n_; j++)
	for (int i=x_sizeKernel/2;i<m_-x_sizeKernel/2;i++) {
	  renormK=0;
	  for (int kj=0; kj<y_sizeKernel/2+(n_-j); kj++)
	    for (int ki=0;ki<x_sizeKernel;ki++)
	      renormK+=kernel[ki][kj];
	  for (int kj=0; kj<y_sizeKernel/2+(n_-j); kj++)
	    for (int ki=0;ki<x_sizeKernel;ki++) 
	      smooA[i][j]+=kernel[ki][kj]*v_[i-x0+ki][j-y0+kj]/renormK;
	}

      // Left column
      for (int i=0; i<x_sizeKernel/2; i++) 
	for (int j=y_sizeKernel/2; j<n_-y_sizeKernel/2; j++) {
	  renormK=0;
	  for (int ki=x_sizeKernel/2-i; ki<x_sizeKernel; ki++)
	    for (int kj=0; kj<y_sizeKernel; kj++)
	      renormK+=kernel[ki][kj];
	  for (int ki=x_sizeKernel/2-i; ki<x_sizeKernel; ki++)
	    for (int kj=0; kj<y_sizeKernel; kj++)
	      smooA[i][j]+=kernel[ki][kj]*v_[i-x0+ki][j-y0+kj]/renormK;
	}
      
      // Right column
      for (int j=y_sizeKernel/2;j<n_-y_sizeKernel/2;j++)
	for (int i=m_-x_sizeKernel/2; i<m_; i++) {
	  renormK=0;
	  for (int kj=0;kj<y_sizeKernel;kj++)
	    for (int ki=0; ki<x_sizeKernel/2+(m_-i); ki++)
	      renormK+=kernel[ki][kj];
	  for (int kj=0;kj<y_sizeKernel;kj++)
	    for (int ki=0; ki<x_sizeKernel/2+(m_-i); ki++)
	      smooA[i][j]+=kernel[ki][kj]*v_[i-x0+ki][j-y0+kj]/renormK;
	}
      
      for (int i=0; i<x_sizeKernel/2; i++)
	for (int j=0; j<n_; j++)
	  smooA[i][j]=v_[i][j];
      for (int i=m_-x_sizeKernel/2; i<m_; i++)
	for (int j=0; j<n_; j++)
	  smooA[i][j]=v_[i][j];

      for (int j=y_sizeKernel/2;j<n_-y_sizeKernel/2;j++)
	for (int i=x_sizeKernel/2;i<m_-x_sizeKernel/2;i++)
	  for (int kj=0;kj<y_sizeKernel;kj++)
	    for (int ki=0;ki<x_sizeKernel;ki++){
	      smooA[i][j]+=kernel[ki][kj]*v_[i-x0+ki][j-y0+kj];
	    }

      /*
      for (int j=0; j<y_sizeKernel/2;j++) {
	for (int i=0; i<x_sizeKernel/2;i++)
	  //	  for (int kj=0;kj<j;kj++)
	  //	    for (int ki=0;ki<i;ki++)
	  //	      smooA[i][j]+=kernel[ki][kj]*v_[i+ki][j+kj];
	  smooA[i][j]=v_[i][j];
	for (int i=m_-x_sizeKernel/2; i<m_; i++)
	  //	  for (int kj=0;kj<j;kj++)
	  //	    for (int ki=0;ki<(m_-i);ki++)
	  //	      smooA[i][j]+=kernel[ki][kj]*v_[i+ki][j+kj];
	  smooA[i][j]=v_[i][j];
      }
      for (int j=n_-y_sizeKernel/2; j=n_; j++) {
	for (int i=0; i<x_sizeKernel/2;i++)
	  //	  for (int kj=0;kj<(n_-j);kj++)
	  //	    for (int ki=0;ki<i;ki++)
	  //	      smooA[i][j]+=kernel[ki][kj]*v_[i+ki][j+kj];
	  smooA[i][j]=v_[i][j];
	for (int i=m_-x_sizeKernel/2; i<m_; i++)
	  //	  for (int kj=0;kj<(n_-j);kj++)
	  //	    for (int ki=0;ki<(m_-i);ki++)
	  //	      smooA[i][j]+=kernel[ki][kj]*v_[i+ki][j+kj];
	  smooA[i][j]=v_[i][j];
      }
      */
      
      return smooA;
    }
  
  template <class T>
    Array2D<T> Array2D<T>::smooth1D(char axis,
				    float old_size, 
				    float new_size) 
    {
      float sig=(new_size/old_size/bcorr/bcorr2);

      sizeKernel=(int)(6.*sig+0.5);
      
      if (sizeKernel%2 == 0)
	sizeKernel++;
      
      Array1D<double> kernel(sizeKernel);
    
      x0=sizeKernel/2;
    
      for (int i=0;i<sizeKernel;i++) {
	kernel[i]=1/sqrt(2.*PI)/sig*exp(-(i-x0)*(i-x0)/2./sig/sig);
      }

      Array2D<T> smooA(m_, n_, 0.);
    
      if (axis == 'y') {
	for (int j=sizeKernel/2;j<n_-sizeKernel/2;j++)
	  for (int i=0;i<m_;i++)
	    for (int ki=0;ki<sizeKernel;ki++)
	      smooA[i][j]+=kernel[ki]*v_[i][j-x0+ki];
      }
      else {
	if (axis == 'x') {
	  for (int i=sizeKernel/2;i<m_-sizeKernel/2;i++)
	    for (int j=0;j<n_;j++)
	      for (int ki=0;ki<sizeKernel;ki++)
		smooA[i][j]+=kernel[ki]*v_[i-x0+ki][j];
	}
      }
      return smooA;
    }

  template <class T>
    Array2D<T> Array2D<T>::smooth1D_noZero(char axis,
					   float old_size, 
					   float new_size) 
    {
      float sig=(new_size/old_size/bcorr/bcorr2);
      sizeKernel=(int)(6.*sig+0.5);
      if (sizeKernel%2 == 0)
	sizeKernel++;
      Array1D<double> kernel(sizeKernel);
      x0=sizeKernel/2;
    
      for (int i=0;i<sizeKernel;i++) {
	kernel[i]=1/sqrt(2.*PI)/sig*exp(-(i-x0)*(i-x0)/2./sig/sig);
      }

      Array2D<T> smooA(m_, n_, 0.);
      //      Array2D<T> A(m_, n_);

      //      for (int j=0;j<n_;j++)
      //	for (int i=0;i<m_;i++)
      //	  A[i][j]=v_[i][j];
      
      int nZero=sizeKernel;
      double renormK;
      if (axis == 'y') {
	for (int j=sizeKernel/2;j<n_-sizeKernel/2;j++)
	  for (int i=0;i<m_;i++) {
	    renormK=0;
	    for (int ki=0;ki<sizeKernel;ki++)
	      if (v_[i][j-x0+ki] == 0)
		renormK+=kernel[ki];
	    for (int ki=0;ki<sizeKernel;ki++)
	      if (v_[i][j-x0+ki] != 0)
		if (renormK != 0)
		  smooA[i][j]+=kernel[ki]*v_[i][j-x0+ki]/renormK;
		else
		  smooA[i][j]+=kernel[ki]*v_[i][j-x0+ki];
	  }
      }

      if (axis == 'x') {
	for (int i=sizeKernel/2;i<m_-sizeKernel/2;i++)
	  for (int j=0;j<n_;j++) {
	    renormK=0;
	    for (int ki=0;ki<sizeKernel;ki++) {
	      if (v_[i-x0+ki][j] == 0) {
		renormK+=kernel[ki];
	      }
	    }
	    renormK=1-renormK;
	    for (int ki=0;ki<sizeKernel;ki++) {
	      if (v_[i-x0+ki][j] != 0) {
		smooA[i][j]+=kernel[ki]*v_[i-x0+ki][j]/renormK;
	      }
	    }
	  }
	}
      return smooA;
    }

}
#endif
