#ifndef ED_hpp
#define ED_hpp

#include <vector>
#include <unordered_map>
#include "SML.hpp"

struct A_Basis_Set {
   std::vector<double>           Val;
   std::vector<long>             Basis;
   std::vector<long>             Site_Constant;
   std::vector<int>              Local_Basis;
   std::unordered_map<long, int> Check_Basis;
   double zero_precision;
};

struct ED_Time {
   double bases;
   double ham;
   double diag;
   double onsite_exp;
   double intersite_exp;
   double total;
};

void ED_Expectation_Intersite_Q1(CRS &M, std::vector<double> &Vec_In, std::vector<long> &Bases_In, std::vector<long> &Bases_Out, std::vector<double> &Out, int site_ref, int tot_site, int p_threads);
void ED_Expectation_Intersite_Q1(CRS &M, std::vector<double> &Vec_In, std::vector<long> &Bases_In, std::vector<long> &Bases_Out_1, std::vector<long> &Bases_Out_2, std::vector<double> &Out, int site_ref, int tot_site, int p_threads);
void ED_Expectation_Onsite(CRS &M, std::vector<double> &Vec, std::vector<long> &Bases, std::vector<double> &Out, int tot_site, int p_threads);
int  ED_Find_Local_Basis(long basis, int site, int dim_onsite);
void ED_Make_Intersite_Element(long basis, int site1, int site2, CRS &M1, CRS &M2, double coeef, int sign, A_Basis_Set &A_Set);
void ED_Make_Onsite_Element(long basis, int site, CRS &M, double coeef, A_Basis_Set &A_Set);
void ED_Matrix_Vector_Product_Q1(CRS &M, std::vector<double> &Vec_In, std::vector<long> &Bases_In, std::vector<double> &Vec_Out, std::vector<long> &Bases_Out, int site, int p_threads);
void ED_Output_Time(ED_Time &Time, std::string file_name);
void ED_Make_Zero_Element(long basis, A_Basis_Set &A_Set);

#endif /* ED_hpp */
