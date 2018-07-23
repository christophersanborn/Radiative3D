// array.hpp
//
// This file defines three pseudo arrays:
//
// Array::Pair
// Array::Triple
// Array::Quad
//
#ifndef ARRAY_H
#define ARRAY_H
//
#include <cstdlib>
#include <iostream>

//________________________________________________________________________
//************************************************************************
// NAMESPACE: Array
// PURPOSE:
//
//  Allow construction of pseudo arrays which can be initialized in the
//  constructor initialization list.
//
//________________________________________________________________________
//************************************************************************

namespace Array {

//////
// CLASS:   ::::  Pair ::::
//
//   Defines an array of two objects of 
//   type T using the template format.
//
template<class T>
class Pair
{
protected:

  // :::::::::::::::::::::::::::::::::
  // ::: Member Data  (Pair Class) :::
  // :::::::::::::::::::::::::::::::::

  T mA;
  T mB;


public:

  // ::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Pair Class) :::
  // ::::::::::::::::::::::::::::::::::

  Pair(T a, T b):
    mA(a), mB(b){};

  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (Pair Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  T & operator[](int index){
    switch(index){
    case (0):
      return mA;
      break;
    case (1):
      return mB;
      break;
    default:
      std::cerr << "Incorrect Pair index" << std::endl;
      exit(1);
      break;
    }
  }

 const T & operator[](int index) const{
    switch(index){
    case (0):
      return mA;
      break;
    case (1):
      return mB;
      break;
    default:
      std::cerr << "Incorrect Pair index" << std::endl;
      exit(1);
      break;
    }
  }
};


//////
// CLASS:   ::::  Triple ::::
//
//   Defines an array of three objects of 
//   type T using the template format.
//
template<class T>
class Triple
{
protected:

  // :::::::::::::::::::::::::::::::::::
  // ::: Member Data  (Triple Class) :::
  // :::::::::::::::::::::::::::::::::::

  T mA;
  T mB;
  T mC;

public:

  // ::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Triple Class) :::
  // ::::::::::::::::::::::::::::::::::::

  Triple(T a, T b, T c):
    mA(a), mB(b), mC(c){};

  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (Triple Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::

  T & operator[](int index){
    switch(index){
    case (0):
      return mA;
      break;
    case (1):
      return mB;
      break;
    case(2):
      return mC;
      break;
    default:
      std::cerr << "Incorrect Triple index" << std::endl;
      exit(1);
      break;
    }
  }

 const T & operator[](int index) const{
    switch(index){
    case (0):
      return mA;
      break;
    case (1):
      return mB;
      break;
    case(2):
      return mC;
      break;
    default:
      std::cerr << "Incorrect Triple index" << std::endl;
      exit(1);
      break;
    }
  }
};


//////
// CLASS:   ::::  Quad ::::
//
//   Defines an array of four objects of 
//   type T using the template format.
//
template<class T>
class Quad
{
protected:

  // :::::::::::::::::::::::::::::::::
  // ::: Member Data  (Quad Class) :::
  // :::::::::::::::::::::::::::::::::

  T mA;
  T mB;
  T mC;
  T mD;

public:

  // ::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Quad Class) :::
  // ::::::::::::::::::::::::::::::::::

  Quad(T a, T b, T c, T d):
    mA(a), mB(b), mC(c), mD(d){};

  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (Quad Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  T & operator[](int index){
    switch(index){
    case (0):
      return mA;
      break;
    case (1):
      return mB;
      break;
    case(2):
      return mC;
      break;
    case(3):
      return mD;
      break;
    default:
      std::cerr << "Incorrect Quad index" << std::endl;
      exit(1);
      break;
    }
  }


  const T & operator[](int index) const{
    switch(index){
    case (0):
      return mA;
      break;
    case (1):
      return mB;
      break;
    case(2):
      return mC;
      break;
    case(3):
      return mD;
      break;
    default:
      std::cerr << "Incorrect Quad index" << std::endl;
      exit(1);
      break;
    }
  }

};


};// End namespace Array
///
//

#endif //#ifndef ARRAY_H_
