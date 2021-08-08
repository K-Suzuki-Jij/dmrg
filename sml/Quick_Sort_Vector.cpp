#include "SML.hpp"

//Ascending order
template <typename T_Target, typename T_1>
void Quick_Sort_Vector(std::vector<T_Target> &Array_Target, std::vector<T_1> &Array1, long left, long right) {
   
   if (right - left <= 1) {
      return;
   }
   
   long     pivot_index = (left + right)/2;
   T_Target pivot       = Array_Target[pivot_index];
   
   std::swap(Array_Target[pivot_index], Array_Target[right - 1]);
   std::swap(Array1[pivot_index]      , Array1[right - 1]      );

   
   long index = left;
   for (long i = left; i < right - 1; i++) {
      if (Array_Target[i] < pivot) {
         std::swap(Array_Target[index], Array_Target[i]);
         std::swap(Array1[index]      , Array1[i]      );
         index++;
      }
   }
  
   std::swap(Array_Target[index], Array_Target[right - 1]);
   std::swap(Array1[index]      , Array1[right - 1]      );

   Quick_Sort_Vector(Array_Target, Array1, left     , index);
   Quick_Sort_Vector(Array_Target, Array1, index + 1, right);
   
}

template void Quick_Sort_Vector<double, int>    (std::vector<double> &, std::vector<int>    &, long left, long right);
template void Quick_Sort_Vector<double, long>   (std::vector<double> &, std::vector<long>   &, long left, long right);
template void Quick_Sort_Vector<double, double> (std::vector<double> &, std::vector<double> &, long left, long right);
template void Quick_Sort_Vector<int, int>       (std::vector<int>    &, std::vector<int>    &, long left, long right);
template void Quick_Sort_Vector<int, long>      (std::vector<int>    &, std::vector<long>   &, long left, long right);
template void Quick_Sort_Vector<int, double>    (std::vector<int>    &, std::vector<double> &, long left, long right);
template void Quick_Sort_Vector<long, int>      (std::vector<long>   &, std::vector<int>    &, long left, long right);
template void Quick_Sort_Vector<long, long>     (std::vector<long>   &, std::vector<long>   &, long left, long right);
template void Quick_Sort_Vector<long, double>   (std::vector<long>   &, std::vector<double> &, long left, long right);


template <typename T_Target, typename T_1, typename T_2>
void Quick_Sort_Vector(std::vector<T_Target> &Array_Target, std::vector<T_1> &Array1, std::vector<T_2> &Array2, long left, long right) {
   
   if (right - left <= 1) {
      return;
   }
   
   long     pivot_index = (left + right)/2;
   T_Target pivot       = Array_Target[pivot_index];
   
   std::swap(Array_Target[pivot_index], Array_Target[right - 1]);
   std::swap(Array1[pivot_index]      , Array1[right - 1]      );
   std::swap(Array2[pivot_index]      , Array2[right - 1]      );

   
   long index = left;
   for (long i = left; i < right - 1; i++) {
      if (Array_Target[i] < pivot) {
         std::swap(Array_Target[index], Array_Target[i]);
         std::swap(Array1[index]      , Array1[i]      );
         std::swap(Array2[index]      , Array2[i]      );
         index++;
      }
   }
  
   std::swap(Array_Target[index], Array_Target[right - 1]);
   std::swap(Array1[index]      , Array1[right - 1]      );
   std::swap(Array2[index]      , Array2[right - 1]      );

   Quick_Sort_Vector(Array_Target, Array1, Array2, left     , index);
   Quick_Sort_Vector(Array_Target, Array1, Array2, index + 1, right);
   
}

template void Quick_Sort_Vector<double, int, int> (std::vector<double> &, std::vector<int> &, std::vector<int> &, long left, long right);
template void Quick_Sort_Vector<int, int, int>    (std::vector<int>    &, std::vector<int> &, std::vector<int> &, long left, long right);

