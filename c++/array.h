/*
*
* Template Numerical Toolkit (TNT): Linear Algebra Module
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
*
* Modified by F. Fraternali for astrophyisical purposes.
*
*/


#ifndef ARRAY_H
#define ARRAY_H



//---------------------------------------------------------------------
// Define this macro if you want  TNT to track some of the out-of-bounds
// indexing. This can encur a small run-time overhead, but is recommended 
// while developing code.  It can be turned off for production runs.
// 
//       #define TNT_BOUNDS_CHECK
//---------------------------------------------------------------------
//

#define TNT_BOUNDS_CHECK

//#include "array_version.h"
//#include "array_math_utils.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "array1d_extra.h"
#include "array2d_extra.h"
#include "array3d_extra.h"

//#include "tnt_sparse_vector.h"
//#include "tnt_sparse_matrix.h"

//#include "tnt_stopwatch.h"
//#include "tnt_subscript.h"
//#include "tnt_vector.h"
//#include "tnt_matrix.h"


#endif
