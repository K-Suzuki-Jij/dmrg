#ifndef Model_1D_TAKLM_hpp
#define Model_1D_TAKLM_hpp

#include <string>
#include "SML.hpp"

struct Model_1D_TAKLM {
   
   int p_thread;    //Number of Parallel Threads
   
   std::string BC;  //Boundary Condition
   
   int tot_sz;        //Total Sz in the Target Sector
   int tot_ele_1;     //# of total electron in an orbital 1
   int tot_ele_2;     //# of total electron in an orbital 2
   int lspin;         //Twice the Magnitude of the Spin
   int num_of_qn;     //Number of Conserved Quantum Number, it must be 3 (total spin and total electron)
   int system_size;   //System Size
   int site_cf_ref;
   int dim_lspin;
   int dim_ele;     //set as 4
   int dim_onsite;  //dimension of the Local Hilbert Space
   int dim_target;  //dimension of the Target Hilbert Space
   double zero_precision;
   
   std::vector<int> Param_Tot_Ele_1, Param_Tot_Ele_2, Param_Tot_Sz;
   
   std::vector<double> Param_t_1, Param_t_2, Param_J_xy_1, Param_J_z_1, Param_J_xy_2, Param_J_z_2, Param_K_xy_1, Param_K_z_1, Param_K_xy_2, Param_K_z_2;
   std::vector<double> Param_I_xy, Param_I_z, Param_U_1, Param_V_1, Param_U_2, Param_V_2, Param_D_z, Param_hc_z_1, Param_hc_z_2, Param_hl_z, Param_mu_1, Param_mu_2;
   double t_1, t_2, J_xy_1, J_z_1, J_xy_2, J_z_2, K_xy_1, K_z_1, K_xy_2, K_z_2, I_xy, I_z, U_1, V_1, U_2, V_2, D_z, hc_z_1, hc_z_2, hl_z, mu_1, mu_2;
   
   CRS Ham_On, SC_1SL_On, SC_2SL_On;
   CRS SxL_On, iSyL_On, SzL_On, SpL_On, SmL_On, SxLSxL_On, SyLSyL_On, SzLSzL_On;
   CRS SxC_1_On, iSyC_1_On, SzC_1_On, SpC_1_On, SmC_1_On, SxC_1SxC_1_On, SyC_1SyC_1_On, SzC_1SzC_1_On;
   CRS SxC_2_On, iSyC_2_On, SzC_2_On, SpC_2_On, SmC_2_On, SxC_2SxC_2_On, SyC_2SyC_2_On, SzC_2SzC_2_On;
   CRS CUp_1_On, CUp_1_D_On, CDown_1_On, CDown_1_D_On, NCUp_1_On, NCDown_1_On, NC_1_On;
   CRS CUp_2_On, CUp_2_D_On, CDown_2_On, CDown_2_D_On, NCUp_2_On, NCDown_2_On, NC_2_On;
   CRS NC_On;
   
   std::vector<long> Site_Constant;
   
   void Check_Parameters();
   void Set_Onsite_Op();
   
   int Find_Basis_Ele_1(int basis) { return basis/dim_lspin/dim_ele;};
   int Find_Basis_Ele_2(int basis) { return (basis/dim_lspin)%dim_ele;};
   int Find_Basis_Lspin(int basis) { return basis%dim_lspin;};
   int Find_Dim_Lspin  (         ) { return lspin + 1;};
   int Find_Dim_Ele    (         ) { return 4;};
   int Find_Dim_Onsite (         ) { return (lspin + 1)*4*4;};
   
   int Find_Site_Sz   (int basis);
   int Find_Site_Ele_1(int basis);
   int Find_Site_Ele_2(int basis);
   int Find_Site_Tot_Ele(int basis) {return Find_Site_Ele_1(basis) + Find_Site_Ele_2(basis);};
   
   void Get_Ham_On       (CRS &M);
   void Get_SC_1SL_On    (CRS &M, double coeef);
   void Get_SC_2SL_On    (CRS &M, double coeef);
   void Get_SxL_On       (CRS &M, double coeef);
   void Get_iSyL_On      (CRS &M, double coeef);
   void Get_SzL_On       (CRS &M, double coeef);
   void Get_SpL_On       (CRS &M, double coeef);
   void Get_SmL_On       (CRS &M, double coeef);
   void Get_SxLSxL_On    (CRS &M, double coeef);
   void Get_SyLSyL_On    (CRS &M, double coeef);
   void Get_SzLSzL_On    (CRS &M, double coeef);
   void Get_SxC_1_On     (CRS &M, double coeef);
   void Get_iSyC_1_On    (CRS &M, double coeef);
   void Get_SzC_1_On     (CRS &M, double coeef);
   void Get_SpC_1_On     (CRS &M, double coeef);
   void Get_SmC_1_On     (CRS &M, double coeef);
   void Get_SxC_1SxC_1_On(CRS &M, double coeef);
   void Get_SyC_1SyC_1_On(CRS &M, double coeef);
   void Get_SzC_1SzC_1_On(CRS &M, double coeef);
   void Get_CUp_1_On     (CRS &M, double coeef);
   void Get_CUp_1_D_On   (CRS &M, double coeef);
   void Get_CDown_1_On   (CRS &M, double coeef);
   void Get_CDown_1_D_On (CRS &M, double coeef);
   void Get_NCUp_1_On    (CRS &M, double coeef);
   void Get_NCDown_1_On  (CRS &M, double coeef);
   void Get_NC_1_On      (CRS &M, double coeef);
   void Get_SxC_2_On     (CRS &M, double coeef);
   void Get_iSyC_2_On    (CRS &M, double coeef);
   void Get_SzC_2_On     (CRS &M, double coeef);
   void Get_SpC_2_On     (CRS &M, double coeef);
   void Get_SmC_2_On     (CRS &M, double coeef);
   void Get_SxC_2SxC_2_On(CRS &M, double coeef);
   void Get_SyC_2SyC_2_On(CRS &M, double coeef);
   void Get_SzC_2SzC_2_On(CRS &M, double coeef);
   void Get_CUp_2_On     (CRS &M, double coeef);
   void Get_CUp_2_D_On   (CRS &M, double coeef);
   void Get_CDown_2_On   (CRS &M, double coeef);
   void Get_CDown_2_D_On (CRS &M, double coeef);
   void Get_NCUp_2_On    (CRS &M, double coeef);
   void Get_NCDown_2_On  (CRS &M, double coeef);
   void Get_NC_2_On      (CRS &M, double coeef);
   void Get_NC_On        (CRS &M, double coeef);

};

#endif /* Model_1D_TAKLM_hpp */
