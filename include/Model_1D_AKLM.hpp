#ifndef Model_1D_XXZ_hpp
#define Model_1D_XXZ_hpp

#include <string>
#include "SML.hpp"

struct Model_1D_AKLM {
   
   int p_thread; //Number of Parallel Threads
   
   std::string BC; //Boundary Condition
   std::string Mat_Type;
   
   int tot_sz;      //Total Sz in the Target Sector
   int lspin;    //Twice the Magnitude of the Spin
   int tot_ele;     //# of total electron in the Target Sector
   int num_of_qn;   //Number of Conserved Quantum Number, it must be 2 (total spin and total electron)
   int system_size; //System Size
   int site_cf_ref;
   int dim_lspin;
   int dim_ele;
   int dim_onsite;  //dimension of the Local Hilbert Space
   int dim_target;  //dimension of the Target Hilbert Space
   double zero_precision;
   
   std::vector<int> Param_Tot_Ele, Param_Tot_Sz;
   
   std::vector<double> Param_t, Param_J_xy, Param_J_z, Param_I_xy, Param_I_z, Param_K_xy, Param_K_z;
   std::vector<double> Param_U, Param_V, Param_D_z, Param_hc_z, Param_hl_z, Param_mu;
   double t, J_xy, J_z, K_xy, K_z, I_xy, I_z, U, V, D_z, hc_z, hl_z, mu;
   
   CRS Ham_On, SCSL_On;
   CRS SxL_On, iSyL_On, SzL_On, SpL_On, SmL_On, SxLSxL_On, SyLSyL_On, SzLSzL_On;
   CRS SxC_On, iSyC_On, SzC_On, SpC_On, SmC_On, SxCSxC_On, SyCSyC_On, SzCSzC_On;
   CRS CUp_On, CUp_D_On, CDown_On, CDown_D_On, NCUp_On, NCDown_On, NC_On;
   
   std::vector<long> Site_Constant;
   
   void Check_Parameters();
   void Set_Onsite_Op();
   
   int Find_Basis_Ele  (int basis) { return basis/dim_lspin; };
   int Find_Basis_Lspin(int basis) { return basis%dim_lspin; };
   int Find_Dim_Lspin  (         ) { return lspin + 1;       };
   int Find_Dim_Ele    (         ) { return 4;               };
   int Find_Dim_Onsite (         ) { return (lspin + 1)*4;   };
   
   int Find_Site_Sz (int basis);
   int Find_Site_Ele(int basis);
      
   void Get_Ham_On    (CRS &M);
   void Get_SCSL_On   (CRS &M, double coeef);
   void Get_SxL_On    (CRS &M, double coeef);
   void Get_iSyL_On   (CRS &M, double coeef);
   void Get_SzL_On    (CRS &M, double coeef);
   void Get_SpL_On    (CRS &M, double coeef);
   void Get_SmL_On    (CRS &M, double coeef);
   void Get_SxLSxL_On (CRS &M, double coeef);
   void Get_SyLSyL_On (CRS &M, double coeef);
   void Get_SzLSzL_On (CRS &M, double coeef);
   void Get_SxC_On    (CRS &M, double coeef);
   void Get_iSyC_On   (CRS &M, double coeef);
   void Get_SzC_On    (CRS &M, double coeef);
   void Get_SpC_On    (CRS &M, double coeef);
   void Get_SmC_On    (CRS &M, double coeef);
   void Get_SxCSxC_On (CRS &M, double coeef);
   void Get_SyCSyC_On (CRS &M, double coeef);
   void Get_SzCSzC_On (CRS &M, double coeef);
   void Get_CUp_On    (CRS &M, double coeef);
   void Get_CUp_D_On  (CRS &M, double coeef);
   void Get_CDown_On  (CRS &M, double coeef);
   void Get_CDown_D_On(CRS &M, double coeef);
   void Get_NCUp_On   (CRS &M, double coeef);
   void Get_NCDown_On (CRS &M, double coeef);
   void Get_NC_On     (CRS &M, double coeef);

};

#endif /* Model_1D_XXZ_hpp */
