/*
 *  arrays3d.h
 *  TestCfits
 *
 *  Created by Antonino Marasco on 11/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

// ary3d.h ///////////////////////////////////////////////////////////
// 2D/3D array classes 
//  double3d , float3d , int3d , short3d , char3d
//  double2d , float2d , int2d , short2d , char2d
//
//  Usage:
//   int x=5,y=4,z=6,i,j,k;                int x=6,y=5,i,j;
//
//   // declaration                        // declaration
//   double3d a(x,y,z);                    float2d a(x,y);
//
//   // assignment                         // assignment
//   for(i=0;i<x;i++)                      for(i=0;i<x;i++) 
//     for(j=0;j<y;j++)                      for(j=0;j<y;j++)
//       for(k=0;i<z;k++)                      a(i,j) = i*2+j;
//         a(i,j,i) = i*j+k;
//
//   // reference                          // reference
//   for(i=0;i<x;i++)                      for(i=0;i<x;i++)
//     for(j=0;j<y;j++)                      for(j=0;j<y;j++)
//       for(k=0;k<z;k++)                      cout << a(i,j) << endl;
//         cout << a(i,j,k) << endl;
//
// NOTE1:
//   It is recommended to precompile this header file. Precompiled version
//  is faster and more reliable.
//     $  makecint -mk Makefile -dl ary3d.dll -H ary3d.h
//     $  make -f Makefile
//
// NOTE2:
//   Only simple run-time array index check is done. Run-time array index
//  check can be turned off by undefining INDEXCHECK macro.
//
// NOTE3:
//   a[i][j][k] notation is available if INDEXOPR macro is defined.
//  This notation is less efficient compared to a(i,j,k). Not recommended.
//
/////////////////////////////////////////////////////////////////////
#include <iostream>

#ifndef ARRAY2D
#define ARRAY2D

//#define INDEXCHECK //boundary check! Sometimes useful, but you loose speed.

#define NEGATIVE_ALLOWED//my version, useful for smoothing an all-sky map. Don't use it if you are a noob.

/////////////////////////////////////////////////////////////////////
// 2 dimentional array class template
/////////////////////////////////////////////////////////////////////

using namespace std;

namespace myarray{
template<class T> class array2d {
 public:
  // constructors
  array2d(int jx,int kx) { 
    maxj=jx; maxk=kx;
    isnew=1;
    p = new T[maxj*maxk];
  }
  array2d(array2d &obj) {copy(obj);}
//#ifdef INDEXOPR
  array2d() { p=NULL; }
  array2d(T* pin,int jx,int kx) {init(pin,jx,kx);}
  void init(T* pin,int jx,int kx) { 
    maxj=jx; maxk=kx;
    isnew=0;
    p = pin;
  } 
//#endif

  // destructor
  ~array2d() {if(isnew&&p) delete[] p;}

  // operator overloading
  array2d& operator=(array2d& obj) {
    if(p==obj.p) return(*this);
    if(isnew&&p) delete[] p;
    copy(obj);
    return(*this);
  }
	
  T& operator()(int j,int k) {
#ifdef NEGATIVE_ALLOWED //my version, to allow for negative indexes
	  myj = j;
	  myk = k;
	  if(j<0) {myj += maxj - 1;}
	  else if(j>=maxj) {myj += 1 - maxj;}
	  if(k<0)
	  {
		  myk = -k;
		  if(myj<=maxj/2) {myj += maxj/2;}
		  else{myj -= maxj/2;}
	  }
	  else if(k>=maxk)
	  {
		  myk = 2*maxk - k - 2;
		  if(myj<=maxj/2) {myj += maxj/2;}
		  else{myj -= maxj/2;}
	  }
	  return( *(p + myj + myk*maxj) ); 
#endif
#ifdef INDEXCHECK
    if(j<0||maxj<=j||k<0||maxk<=k) {
      cerr << "Bad index ("<<j<<","<<k<<") > ("<<maxj<<","<<maxk<<")"<< endl;
      return(*p);
    }
#endif
	  return( *(p + j + k*maxj) );
  }

//#ifdef INDEXOPR
  T* operator[](int k) {return(p+k*maxj);}
//#endif

  //friend bool operator==(array2d<T>& a,array2d<T>& b);
  //friend bool operator!=(array2d<T>& a,array2d<T>& b);
  //friend ostream& operator<<(ostream& ost,array2d<T>& a);

	
	
 private:
  T* p;
  int maxj,maxk;
  int myj, myk;
  int isnew;

  // utility function
  void copy(array2d& obj) {
    maxj=obj.maxj; maxk=obj.maxk;
    if(isnew) {
      isnew=1;
      p = new T[maxj*maxk];
      memcpy(p,obj.p,sizeof(T)*maxj*maxk);
    }
    else {
      isnew=0;
      p = obj.p;
    }
  }
};

template<class T>
bool operator==(array2d<T>& a,array2d<T>& b) {
  if(a.maxj!=b.maxj || a.maxk!=b.maxk) return(false);
  int i,max=a.maxj*a.maxk;
  for(i=0;i<max;i++) {
    if(a.p[i]!=b.p[i]) return(false);
  }
  return(true);
}

template<class T>
bool operator!=(array2d<T>& a,array2d<T>& b) {
  return(!(a==b));
}

template<class T>
ostream& operator<<(ostream& ost,array2d<T>& a) {
  int j,k;
  ost << '(' ;
  for(j=0;j<a.maxj;j++) {
    if(j) ost << ',';
    ost << '(' ;
    for(k=0;k<a.maxk;k++) {
      if(k) ost << ',';
      ost << a(j,k);
    }
    ost << ')' ;
  }
  ost << ')' ;
  return(ost);
}
#endif

#ifndef ARRAY3D
#define ARRAY3D
/////////////////////////////////////////////////////////////////////
// 3 dimentional array class template
/////////////////////////////////////////////////////////////////////
template<class T> class array3d {
 public:
  // constructor
  array3d(int ix,int jx,int kx) { 
    maxi=ix; maxj=jx; maxk=kx;
    isnew=1;
    p = new T[maxi*maxj*maxk];
  }
  array3d(array3d& obj) {copy(obj);}
//#ifdef INDEXOPR
  array3d() {p=NULL;}
  array3d(T* pin,int ix,int jx,int kx) {init(pin,ix,jx,kx);}
  void init(T* pin,int ix,int jx,int kx) { 
    maxi=ix; maxj=jx; maxk=kx;
    isnew=0;
    p = pin;
  }
//#endif

  // destructor
  ~array3d() {if(isnew&&p) delete[] p;}

  // operator overloading
  array3d& operator=(array3d& obj) {
    if(p==obj.p) return(*this);
    if(isnew&&p) delete[] p;
    copy(obj);
    return(*this);
  }
	
  T& operator()(int i,int j,int k) {
#ifdef NEGATIVE_ALLOWED  
	  myi = i;
	  myj = j;
	  myk = k;
	  if(i<0) {myi += maxi - 1;}
	  else if(i>=maxi) {myi += 1 - maxi;}
	  if(k<0) {myk += maxk - 1;}
	  else if(k>=maxk) {myk += 1 - maxk;}
	  if(j<0)
	  {
		  myj = -j;
		  if(myi<=maxi/2) {myi += maxi/2;}
		  else{myi -= maxi/2;}
	  }
	  else if(j>=maxj)
	  {
		  myj = 2*maxj - j - 2;
		  if(myi<=maxi/2) {myi += maxi/2;}
		  else{myi -= maxi/2;}
	  }
	  return( *(p + myi + maxi*myj + maxi*maxj*myk) );
#endif
#ifdef INDEXCHECK
    if(i<0||maxi<=i||j<0||maxj<=j||k<0||maxk<=k) 
	{
      cerr << "Bad index ("<<i<<","<<j<<","<<k<<") > (";
      cerr <<maxi<<","<<maxj<<","<<maxk<<")"<< endl;
      return(*p);
    }
#endif
	  return( *(p + i + maxi*j + maxi*maxj*k) );
  }

//#ifdef INDEXOPR
  array2d<T>& operator[](int k) {
    subary.init(p+maxi*maxj*k,maxi,maxj);//NOT REALLY SURE ABOUT THIS
    return(subary);
  }
//sizes	
	int xsize(){return maxi;}
	int ysize(){return maxj;}
	int zsize(){return maxk;}
	
//#endif
  //friend bool operator==(array3d<T>& a,array3d<T>& b);
  //friend bool operator!=(array3d<T>& a,array3d<T>& b);
  //friend ostream& operator<<(ostream& ost,array3d<T>& a);

 private:
  T* p;
  int maxi,maxj,maxk;
  int myi, myj, myk;
  int isnew;
//#ifdef INDEXOPR
  array2d<T> subary;
//#endif

  // utility function
  void copy(array3d& obj) {
    maxi=obj.maxi; maxj=obj.maxj; maxk=obj.maxk;
    if(isnew) {
      isnew=1;
      p = new T[maxi*maxj*maxk];
      memcpy(p,obj.p,sizeof(T)*maxi*maxj*maxk);
    }
    else {
      isnew=0;
      p = obj.p;
    }
  }
};

template<class T>
bool operator==(array3d<T>& a,array3d<T>& b) {
  if(a.maxi!=b.maxi || a.maxj!=b.maxj || a.maxk!=b.maxk) return(false);
  int i,max=a.maxi*a.maxj*a.maxk;
  for(i=0;i<max;i++) {
    if(a.p[i]!=b.p[i]) return(false);
  }
  return(true);
}
/*
template<class T>
bool operator=(array3d<T>& a,array3d<T>& b) {
	a.maxi=b.maxi;a.maxj=b.maxj;a.maxk=b.maxk;
	A.p = new T[maxi*maxj*maxk];
	memcpy(A.p,B.p,sizeof(T)*maxi*maxj*maxk);
}
*/	
template<class T>
bool operator!=(array3d<T>& a,array3d<T>& b) {
  return(!(a==b));
}

template<class T>
ostream& operator<<(ostream& ost,array3d<T>& a) {
  int i,j,k;
  ost << '(' ;
  for(i=0;i<a.maxi;i++) {
    if(i) ost << ',';
    ost << '(' ;
    for(j=0;j<a.maxj;j++) {
      if(j) ost << ',';
      ost << '(' ;
      for(k=0;k<a.maxk;k++) {
        if(k) ost << ',';
        ost << a(i,j,k);
      }
      ost << ')' ;
    }
    ost << ')' ;
  }
  ost << ')' ;
  return(ost);
}
/////////////////////////////////////////////////////////////////////
// instantiation of class template as typedef
/////////////////////////////////////////////////////////////////////
typedef array3d<double> double3d;
typedef array2d<double> double2d;
typedef array3d<float>  float3d;
typedef array2d<float>  float2d;
typedef array3d<int>    int3d;
typedef array2d<int>    int2d;
typedef array3d<short>  short3d;
typedef array2d<short>  short2d;
typedef array3d<char>   char3d;
typedef array2d<char>   char2d;
}

#endif