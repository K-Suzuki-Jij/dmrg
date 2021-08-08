#ifndef Model_1D_XXZ_hpp
#define Model_1D_XXZ_hpp

#include <string>
#include "SML.hpp"

struct Model_1D_XXZ {
   
   int p_thread; //Number of Parallel Threads
   
   std::string BC; //Boundary Condition
   std::string Mat_Type;
   
   int tot_sz;      //Total Sz of the Target Sector
   int num_of_qn;
   int spin;        //Twice the Magnitude of the Spin
   int system_size; //System Size
   int site_cf_ref;
   int dim_onsite;  //dimension of the Local Hilbert Space
   int dim_target;  //dimension of the Target Hilbert Space
   double zero_precision;
   
   std::vector<int> Param_Tot_Sz;
   std::vector<double> Param_J_xy, Param_J_z, Param_D_z, Param_h_z;
   double J_xy,J_z,D_z,h_z;
   
   CRS Ham_On;
   CRS Sx_On, iSy_On, Sz_On, Sp_On, Sm_On;
   CRS SxSx_On, SySy_On, SzSz_On;
   
   std::vector<long> Site_Constant;
   
   void Check_Parameters();
   void Set_Onsite_Op();
   void Output_Parameters(std::string file_name);
   
   int  Find_Dim_Onsite();
   int  Find_Dim_Target(int target_sz);
   int  Find_Site_Sz(int basis);
   
   void Get_Ham_On (CRS &M);
   void Get_Sx_On  (CRS &M, double coeef);
   void Get_iSy_On (CRS &M, double coeef);
   void Get_Sz_On  (CRS &M, double coeef);
   void Get_Sp_On  (CRS &M, double coeef);
   void Get_Sm_On  (CRS &M, double coeef);
   void Get_SxSx_On(CRS &M, double coeef);
   void Get_SySy_On(CRS &M, double coeef);
   void Get_SzSz_On(CRS &M, double coeef);
   
   void Output_Onsite_Values(std::vector<double> &Val, std::string file_name);
   void Output_Intersite_Values(std::vector<double> &Val, std::string file_name);
   void Output_Average_Values(double val, std::string file_name);
   void Output_Average_Values(std::vector<double> &Val, std::string file_name);
   
};

#endif /* Model_1D_XXZ_hpp */
