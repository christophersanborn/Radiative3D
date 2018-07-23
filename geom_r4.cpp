// geom_r4.cpp
//
#include "geom_r4.hpp"  /* picks up base and r3 as well */
#include <iostream>
#include <cmath>
using namespace R4;
using namespace std;

// IN THIS FILE:
// Class Implementations for:
//
//   o  Classes in the R4:: namespace
//

Column Matrix::SolveAXB(Column B) const
{
  int numRows = 4; 
 
  AugmentedRow row[4];
 
  row[0] = AugmentedRow(GetRow(0),B.x(0));
  row[1] = AugmentedRow(GetRow(1),B.x(1));
  row[2] = AugmentedRow(GetRow(2),B.x(2));
  row[3] = AugmentedRow(GetRow(3),B.x(3));
 
  //Upper Triangular
  for(int i = 0; i < numRows; i++)
    {
      //Swap Rows
      AugmentedRow largest = row[i];
      int largestNumber = i;
      for(int k = i; k < numRows; k++){
        if(abs(row[k].x(i)) > abs(row[i].x(i))){
          largest = row[k];
          largestNumber = k;
        }
      }

      if(largestNumber != i){
        row[largestNumber] = row[i];
        row[i] = largest;
      }

      if(row[i].x(i)!= 0){
        row[i] *= 1/row[i].x(i);
        
        for(int j = 1; j < numRows - i; j++)
          {
            row[i+j] -= (row[i].ScaledBy(row[i+j].x(i)));
          }

      }

    }
 
    
  //Reduced Echelon 
  for(int j = numRows-1; j >= 1; j--)
    {
      for(int i = 0; i <= j-1; i++)
        {
          row[i] -= (row[j].ScaledBy(row[i].x(j)));
        }
    }

  return Column(row[0].x(4), row[1].x(4), row[2].x(4), row[3].x(4));
}
