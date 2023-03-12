// Put an assert to check if x is NULL, this is to catch
// program "logic" errors early. Even though delete works
// fine with NULL by using assert you are actually catching
// "bad code" very early

// Defining Zap using templates
// Use zap instead of delete as this will be very clean

#ifndef ZAP_H
#define ZAP_H

#include <assert.h>

template <class T>
inline void zap(T & x)
{
	{assert(x != NULL);}
	delete x;
	x = NULL;
}

// In C++ the reason there are 2 forms of the delete operator is - because
// there is no way for C++ to tell the difference between a pointer to
// an object and a pointer to an array of objects. The delete operator
// relies on the programmer using "[]" to tell the two apart.
// Hence, we need to define zaparr function below.
// To delete array of pointers
template <class T>
inline void zaparr(T & x)
{
	{assert(x != NULL);}
	delete [] x;
	x = NULL;
}

#endif