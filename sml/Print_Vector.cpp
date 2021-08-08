#include <iostream>
#include <iomanip>
#include <ios>
#include "SML.hpp"

template <typename T>
void Print_Vector(std::vector<T> &Vec, std::string Name) {
   
   std::cout << std::fixed << std::setprecision(15);
   for (int i = 0; i < (int)Vec.size(); i++) {
      std::cout << Name << "[";
      std::cout << std::left << std::noshowpos << std::setw(2) << i << "]=";
      std::cout << std::showpos << Vec[i] << "\n";
   }
   
}

template void Print_Vector<double>        (std::vector<double>        &, std::string Name);
template void Print_Vector<long>          (std::vector<long>          &, std::string Name);
template void Print_Vector<unsigned long> (std::vector<unsigned long> &, std::string Name);
template void Print_Vector<int>           (std::vector<int>           &, std::string Name);
template void Print_Vector<unsigned int>  (std::vector<unsigned int>  &, std::string Name);
template void Print_Vector<short>         (std::vector<short>         &, std::string Name);


template <typename T>
void Print_Vector(std::vector<std::vector<T>> &Vec, std::string Name) {
   
   std::cout << std::fixed << std::setprecision(15);
   for (int i = 0; i < (int)Vec.size(); i++) {
      for (int j = 0; j < (int)Vec[i].size(); j++) {
         std::cout << Name << "[";
         std::cout << std::left << std::noshowpos << std::setw(2) << i << "][" << std::setw(2) << j;
         std::cout << "]=" << std::showpos << Vec[i][j] << "\n";
      }
   }
   
}

template void Print_Vector<double>        (std::vector<std::vector<double>>        &, std::string Name);
template void Print_Vector<long>          (std::vector<std::vector<long>>          &, std::string Name);
template void Print_Vector<unsigned long> (std::vector<std::vector<unsigned long>> &, std::string Name);
template void Print_Vector<int>           (std::vector<std::vector<int>>           &, std::string Name);
template void Print_Vector<unsigned int>  (std::vector<std::vector<unsigned int>>  &, std::string Name);
template void Print_Vector<short>         (std::vector<std::vector<short>>         &, std::string Name);
