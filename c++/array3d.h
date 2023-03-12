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



#ifndef TNT_ARRAY3D_H
#define TNT_ARRAY3D_H

#include <cstdlib>
#include <iostream>
#ifdef TNT_BOUNDS_CHECK
#include <assert.h>
#endif

#include "array1d.h"
#include "array2d.h"

extern const double PI;
extern bool verbose;

namespace ffarray
{

  template <class T>
    class Array3D 
    {
    private:
      Array1D<T> data_;
      Array2D<T*> v_;
      int m_;
      int n_;
      int g_;
        
    public:      
      typedef         T   value_type;      
      Array3D();
      Array3D(int m, int n, int g);
      Array3D(int m, int n, int g,  T val);
      Array3D(int m, int n, int g, T *a);
      
      inline operator T***();
      inline operator const T***();
      inline Array3D(const Array3D &A);
      inline Array3D & operator=(const T &a);
      inline Array3D & operator=(const Array3D &A);
      inline Array3D & ref(const Array3D &A);
      Array3D copy() const;
      Array3D & inject(const Array3D & A);
      
      inline T** operator[](int i);
      inline const T* const * operator[](int i) const;
      inline int dim1() const;
      inline int dim2() const;
      inline int dim3() const;
      ~Array3D();
      
      inline int ref_count(){ return data_.ref_count(); }
      Array3D subarray(int i0, int i1, int j0, int j1, 
		       int k0, int k1);
      Array3D smooth2D(float, float, float, float);
      Array3D smooth2D(int );
      Array3D smoothMW(float, float, float, float);
      void Hanning(int );
    };
  
  /*
CONSTRUCTORS
  */

  template <class T>
    Array3D<T>::Array3D() : data_(), v_(), m_(0), n_(0) {}
  
  template <class T>
    Array3D<T>::Array3D(const Array3D<T> &A) : data_(A.data_), 
    v_(A.v_), m_(A.m_), n_(A.n_), g_(A.g_)
    {
    }
    
  template <class T>
    Array3D<T>::Array3D(int m, int n, int g) : data_(m*n*g), v_(m,n),
    m_(m), n_(n), g_(g)
    {
      
      if (m>0 && n>0 && g>0)
	{
	  T* p = & (data_[0]);
	  int ng = n_*g_;
	  
	  for (int i=0; i<m_; i++)
	    {	
	      T* ping = p+ i*ng;
	      for (int j=0; j<n; j++)
		v_[i][j] = ping + j*g_;
	    }
	}
    }
    
  template <class T>
    Array3D<T>::Array3D(int m, int n, int g, T val) : data_(m*n*g, val), 
    v_(m,n), m_(m), n_(n), g_(g)
    {
      if (m>0 && n>0 && g>0){
	T* p = & (data_[0]);
	int ng = n_*g_;
	
	for (int i=0; i<m_; i++) {	
	  T* ping = p+ i*ng;
	  for (int j=0; j<n; j++)
	    v_[i][j] = ping + j*g_;
	}
      }
    }
    
  template <class T>
    Array3D<T>::Array3D(int m, int n, int g, T* a) : 
		data_(m*n*g, a), v_(m,n), m_(m), n_(n), g_(g)
    {
      if (m>0 && n>0 && g>0) {
	T* p = & (data_[0]);
	int ng = n_*g_;
	for (int i=0; i<m_; i++) {	
	  T* ping = p+ i*ng;
	  for (int j=0; j<n; j++)
	    v_[i][j] = ping + j*g_;
	}
      }
    }

  /*
OPERATORS
  */  
  
  template <class T>
    inline T** Array3D<T>::operator[](int i) 
    { 
#ifdef TNT_BOUNDS_CHECK
      assert(i >= 0);
      assert(i < m_);
#endif
      return v_[i]; 
    }
  
  template <class T>
    inline const T* const * Array3D<T>::operator[](int i) const 
    { return v_[i]; }
  
  template <class T>
  Array3D<T> & Array3D<T>::operator=(const T &a)
  {
    for (int i=0; i<m_; i++)
      for (int j=0; j<n_; j++)
	for (int k=0; k<g_; k++)
	  v_[i][j][k] = a;
    return *this;
  }

  template <class T>
    Array3D<T> Array3D<T>::copy() const
    {
      Array3D A(m_, n_, g_);
      for (int i=0; i<m_; i++)
	for (int j=0; j<n_; j++)
	  for (int k=0; k<g_; k++)
	    A.v_[i][j][k] = v_[i][j][k];
      return A;
    }

  template <class T>
    Array3D<T> & Array3D<T>::inject(const Array3D &A)
    {
      if (A.m_ == m_ &&  A.n_ == n_ && A.g_ == g_)
	for (int i=0; i<m_; i++)
	  for (int j=0; j<n_; j++)
	    for (int k=0; k<g_; k++)
	      v_[i][j][k] = A.v_[i][j][k];
      return *this;
    }
  
  template <class T>
    Array3D<T> & Array3D<T>::ref(const Array3D<T> &A)
    {
      if (this != &A) {
	m_ = A.m_;
	n_ = A.n_;
	g_ = A.g_;
	v_ = A.v_;
	data_ = A.data_;
      }
      return *this;
    }

  template <class T>
    Array3D<T> & Array3D<T>::operator=(const Array3D<T> &A)
  {
    return ref(A);
  }

  template <class T>
    inline int Array3D<T>::dim1() const { return m_; }
  
  template <class T>
    inline int Array3D<T>::dim2() const { return n_; }
  
  template <class T>
    inline int Array3D<T>::dim3() const { return g_; }
  
  template <class T>
    Array3D<T>::~Array3D() 
    {
      // **************** DELETE ARRAYS??? ***************
      //      for (int i=0;i<x_size;i++)
      //	delete [] array[i];
      //      delete [] array;
    }
  
  template <class T>
    inline Array3D<T>::operator T***()
    {
      return v_;
    }

  template <class T>
    inline Array3D<T>::operator const T***()
    {
      return v_;
    }

  /*
SUBARRAY
  */
  
  template <class T>
    Array3D<T> Array3D<T>::subarray(int i0, int i1, int j0,
				    int j1, int k0, int k1)
    {
      /* check that ranges are valid. */
      if (!( 0 <= i0 && i0 <= i1 && i1 < m_ &&
	     0 <= j0 && j0 <= j1 && j1 < n_ &&
	     0 <= k0 && k0 <= k1 && k1 < g_))
	return Array3D<T>();  /* null array */

      Array3D<T> A;

      A.data_ = data_;
      A.m_ = i1-i0+1;
      A.n_ = j1-j0+1;
      A.g_ = k1-k0+1;
      A.v_ = Array2D<T*>(A.m_,A.n_);
      T* p = &(data_[0]) + i0*n_*g_ + j0*g_ + k0; 

      for (int i=0; i<A.m_; i++) {
	T* ping = p + i*n_*g_;
	for (int j=0; j<A.n_; j++)
	  A.v_[i][j] = ping + j*g_ ;
      }
      return A;
    }
  
  template <class T>
    Array3D<T> Array3D<T>::smooth2D(float old_x, 
				    float old_y, 
				    float new_x, 
				    float new_y)
    {
      //      Array2D<double> makeKernel(float, float, float, float);

      Array2D<double> kern;
      Array2D<double> kernel=kern.makeKernel(old_x, old_y, new_x, new_y);

      Array3D<T> smooCube(m_, n_, g_);
      
      int x_sizeKernel=kern.x_sizeKernel;
      int y_sizeKernel=kern.y_sizeKernel;
      
      if(verbose) cout << "Smoothing final cube" << endl;
      for (int k=0; k<g_; k++){
	for (int j=y_sizeKernel/2;j<n_-y_sizeKernel/2;j++)
	  for (int i=x_sizeKernel/2;i<m_-x_sizeKernel/2;i++)
	    for (int kj=0;kj<y_sizeKernel;kj++)
	      for (int ki=0;ki<x_sizeKernel;ki++) {
		if (v_[i-kern.x0+ki][j-kern.y0+kj][k] > 0)
		  smooCube[i][j][k]+=
		    kernel[ki][kj]*v_[i-kern.x0+ki][j-kern.y0+kj][k];
	      }
      }
      return smooCube;
    }

  template <class T>
    Array3D<T> Array3D<T>::smooth2D(int ratio)
    {

      Array2D<double> kern;
      Array2D<double> kernel=kern.makeKernel(1, 1, ratio, ratio);

      Array3D<T> smooCube(m_, n_, g_);
      
      int x_sizeKernel=kern.x_sizeKernel;
      int y_sizeKernel=kern.y_sizeKernel;
      
      cout << "Smoothing final cube" << endl;
      for (int k=0; k<g_; k++){
	for (int j=y_sizeKernel/2;j<n_-y_sizeKernel/2;j++)
	  for (int i=x_sizeKernel/2;i<m_-x_sizeKernel/2;i++)
	    for (int kj=0;kj<y_sizeKernel;kj++)
	      for (int ki=0;ki<x_sizeKernel;ki++) {
		if (v_[i-kern.x0+ki][j-kern.y0+kj][k] > 0)
		  smooCube[i][j][k]+=
		    kernel[ki][kj]*v_[i-kern.x0+ki][j-kern.y0+kj][k];
	      }
      }
      return smooCube;
    }

  template <class T>
    Array3D<T> Array3D<T>::smoothMW(float old_x, 
				    float old_y, 
				    float new_x, 
				    float new_y)
    {
      // PB: It assumes that the old beam is = to the pixel size

      const float bcorr=2.354820045; // FWHM = 2.35 sigma
      const float bcorr2=1.067; // for Gaussian convolution
      int x0, y0, ii, jj;
      int n_sigma = 6;

      Array3D<T> smooCube(m_, n_, g_);

      float sig_y=fabs(new_y/old_y/bcorr/bcorr2), sig_x;

      unsigned long y_sizeKernel=(int)(n_sigma*sig_y+0.5), x_sizeKernel;
      if (y_sizeKernel%2 == 0)
	y_sizeKernel++;

      y0=y_sizeKernel/2;

      float b, bmax; // latitude

      cout << "Smoothing final cube" << endl;

      //      for (int k=0; k<g_; k++){
      for (int k=0; k<g_; k++){

	cout << "channel : " << k <<endl;
	//	for (int j=y_sizeKernel/2;j<n_-y_sizeKernel/2;j++) {
	for (int j=0;j<n_;j++) {

	  // calculating latitude correction
	  b = (float)(j-n_/2)/float(n_)*PI;
	  //	  std::cout << b << endl;
	  bmax = abs(b)+sig_y*n_sigma/2.*PI/180.;
	  if (bmax > PI/2.)
	    bmax = PI/2.;

	  sig_x=fabs((new_x/old_x/bcorr/bcorr2)/cos(bmax));
	  x_sizeKernel=(int)(n_sigma*sig_x+0.5);
	  
	  //	  cout << b << " " << bmax << " " << sig_x << " " << x_sizeKernel << endl;

	  if (x_sizeKernel > m_) {
	    x_sizeKernel = m_;
	  }

	  //	  cout << j << " " << x_sizeKernel << " " << intensityCorr << endl;

	  x0=x_sizeKernel/2;
	  
	  //	  Array2D<double> kernel(x_sizeKernel, y_sizeKernel);
	  
	  //	  for (int i=x_sizeKernel/2;i<m_-x_sizeKernel/2;i++) {
	  for (int i=0;i<m_;i++) {
	    
	    for (int kj=0;kj<y_sizeKernel;kj++)
	      
	      for (int ki=0;ki<x_sizeKernel;ki++) {
		
		b = (float)(j-n_/2)/float(n_)*PI;

		sig_x=fabs(new_x/old_x/bcorr/bcorr2)/cos(b);

		ii = i-x0+ki;
		jj = j-y0+kj;

		if (ii < 0)
		  ii = ii+m_;
		if (ii >= m_)
		  ii = ii-m_;
		if (jj < 0)
		  jj = jj+n_;
		if (jj >= n_)
		  jj = jj-n_;
		//		cout << j << " " << i << " " << ki << " " << x_sizeKernel << " " << ii << endl;

		if (v_[ii][jj][k] != 0.0) {

		  // correction for edge points (sphere)
		  //		  cout << i << " " << ki << " " << m_ << endl;
	      
		  smooCube[i][j][k]+=
		    1/2./PI/sig_x/sig_y*exp(-(ki-x0)*(ki-x0)/2./sig_x/sig_x)
		    *exp(-(kj-y0)*(kj-y0)/2./sig_y/sig_y)
		    *v_[ii][jj][k];
		  
		}

	      }
	  }
	}
      }
      
      return smooCube;
    }
  
  template <class T>
    void Array3D<T>::Hanning(int n_smooth)
    {
      //      array = new Type*[x_size];
      //      for (int i=0;i<x_size;i++)
      //	array[i]=new Type[y_size];
  
      if(verbose) cout << "Smoothing array in velocity..."  <<endl;

      Array3D<T> temp (m_, n_, g_);
	  for (int k=0; k<g_; k++) 
	  {
		  for (int i=0; i<m_; i++) 
		  {
			  for (int j=0; j<n_; j++) temp[i][j][k]=v_[i][j][k];
		  }
	  }	
	  if (n_smooth == 3) 
	  {
		  for (int k=1; k<g_-1; k++) 
		  {
			  for (int i=0; i<m_; i++) 
			  {
				  for (int j=0; j<n_; j++) v_[i][j][k]=.5*temp[i][j][k]+.25*temp[i][j][k-1]+.25*temp[i][j][k+1];
			  }
		  }
      }
      else {
		  cout << "Hanning smoothing for n>3 not implemented yet..." << endl;
      }
     }
} 

#endif

