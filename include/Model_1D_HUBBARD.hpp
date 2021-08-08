#ifndef Model_1D_HUBBARD_hpp
#define Model_1D_HUBBARD_hpp

#include <string>
#include "SML.hpp"

struct Model_1D_HUBBARD {
   
   int p_thread; //Number of Parallel Threads
   
   std::string BC; //Boundary Condition
   std::string Mat_Type;
   
   int tot_sz;      //Total Sz in the Target Sector
   int tot_ele;     //# of total electron in the Target Sector
   int num_of_qn;   //Number of Conserved Quantum Number, it must be 2 (total spin and total electron)
   int system_size; //System Size
   int site_cf_ref;
   int dim_onsite;  //dimension of the Local Hilbert Space
   int dim_target;  //dimension of the Target Hilbert Space
   double zero_precision;
   
   std::vector<int> Param_Tot_Ele, Param_Tot_Sz;
   std::vector<double> Param_t, Param_J_xy, Param_J_z, Param_U, Param_V, Param_h_z, Param_mu;
   double t, J_xy, J_z, U, V, h_z, mu;
   
   CRS Ham_On;
   CRS Sx_On, iSy_On, Sz_On, Sp_On, Sm_On, SxSx_On, SySy_On, SzSz_On;
   CRS CUp_On, CUp_D_On, CDown_On, CDown_D_On, NCUp_On, NCDown_On, NC_On;
   
   std::vector<long> Site_Constant;
   
   void Check_Parameters();
   void Set_Onsite_Op();
   
   int Find_Dim_Onsite() {return 4;};
   
   int Find_Site_Sz (int basis);
   int Find_Site_Ele(int basis);
   
   void Get_Ham_On    (CRS &M);
   void Get_Sx_On     (CRS &M, double coeef);
   void Get_iSy_On    (CRS &M, double coeef);
   void Get_Sz_On     (CRS &M, double coeef);
   void Get_Sp_On     (CRS &M, double coeef);
   void Get_Sm_On     (CRS &M, double coeef);
   void Get_SxSx_On   (CRS &M, double coeef);
   void Get_SySy_On   (CRS &M, double coeef);
   void Get_SzSz_On   (CRS &M, double coeef);
   void Get_CUp_On    (CRS &M, double coeef);
   void Get_CUp_D_On  (CRS &M, double coeef);
   void Get_CDown_On  (CRS &M, double coeef);
   void Get_CDown_D_On(CRS &M, double coeef);
   void Get_NCUp_On   (CRS &M, double coeef);
   void Get_NCDown_On (CRS &M, double coeef);
   void Get_NC_On     (CRS &M, double coeef);

};

#endif /* Model_1D_HUBBARD_hpp */
