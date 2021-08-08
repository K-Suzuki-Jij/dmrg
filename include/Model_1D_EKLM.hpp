#ifndef Model_1D_EKLM_hpp
#define Model_1D_EKLM_hpp

#include <cmath>
#include <string>
#include "SML.hpp"


struct Model_1D_EKLM {
   
   std::string BC;
   
   int p_threads;
   int num_ele_orbit;
   int num_lspin_orbit;
   int system_size;
   int site_cf_ref;
   int num_of_qn;
   int dim_onsite;
   int dim_ele_orbit = 4;
   double zero_precision;
   
   int tot_sz;
   std::vector<int> Tot_Ele;
   std::vector<int> Magnitude_2LSpin;
   
   std::vector<double> Ele_Ele_t_site, Ele_Ele_V_site, LSpin_LSpin_Jxy_site, LSpin_LSpin_Jz_site;
   double ele_mu, ele_hz, ele_U, lspin_Dz, lspin_hz, ele_lspin_Jxy_orbit, ele_lspin_Jz_orbit, ele_ele_V_orbit, lspin_lspin_Jxy_orbit, lspin_lspin_Jz_orbit;
   
   CRS Ham_On, NC_Tot_On;
   std::vector<CRS> CUp_On, CUp_D_On, CDown_On, CDown_D_On, NC_Up_On, NC_Down_On, NC_On;
   std::vector<CRS> SpC_On, SmC_On, SzC_On, SxC_On, iSyC_On;
   std::vector<CRS> SpL_On, SmL_On, SzL_On, SxL_On, iSyL_On;
   std::vector<std::vector<CRS>> SzLSzL_On, SCSL_On;
   
   void Get_Ham_On       (CRS &M);
   void Get_CUp_On       (CRS &M, int target_ele_orbit  , double coeef);
   void Get_CUp_D_On     (CRS &M, int target_ele_orbit  , double coeef);
   void Get_CDown_On     (CRS &M, int target_ele_orbit  , double coeef);
   void Get_CDown_D_On   (CRS &M, int target_ele_orbit  , double coeef);
   void Get_NC_Up_On     (CRS &M, int target_ele_orbit  , double coeef);
   void Get_NC_Down_On   (CRS &M, int target_ele_orbit  , double coeef);
   void Get_NC_On        (CRS &M, int target_ele_orbit  , double coeef);
   void Get_SpC_On       (CRS &M, int target_ele_orbit  , double coeef);
   void Get_SmC_On       (CRS &M, int target_ele_orbit  , double coeef);
   void Get_SzC_On       (CRS &M, int target_ele_orbit  , double coeef);
   void Get_SxC_On       (CRS &M, int target_ele_orbit  , double coeef);
   void Get_iSyC_On      (CRS &M, int target_ele_orbit  , double coeef);
   void Get_SpL_On       (CRS &M, int target_lspin_orbit, double coeef);
   void Get_SmL_On       (CRS &M, int target_lspin_orbit, double coeef);
   void Get_SzL_On       (CRS &M, int target_lspin_orbit, double coeef);
   void Get_SxL_On       (CRS &M, int target_lspin_orbit, double coeef);
   void Get_iSyL_On      (CRS &M, int target_lspin_orbit, double coeef);
   void Get_NC_Tot_On    (CRS &M, double coeef);
   void Get_SzC_Tot_On   (CRS &M, double coeef);
   void Get_SzL_Tot_On   (CRS &M, double coeef);
   void Get_SzLSzL_Tot_On(CRS &M, double coeef);
   void Get_SzLSzL_On    (CRS &M, int target_lspin_orbit_1, int target_lspin_orbit_2, double coeef);
   void Get_SCSL_On      (CRS &M, int target_ele_orbit    , int target_lspin_orbit  , double coeef);
   void Get_Onsite_Operator();
   
   int Find_Dim_Ele   ();
   int Find_Dim_LSpin ();
   int Find_Dim_Onsite();
   int Find_Basis_Ele  (int basis) { return basis%Find_Dim_Ele(); }
   int Find_Basis_LSpin(int basis) { return basis/Find_Dim_Ele(); }
   int Find_Basis_Ele  (int basis, int target_ele_orbit  );
   int Find_Basis_LSpin(int basis, int target_lspin_orbit);
   int Find_Site_Sz(int basis);
   int Find_Site_Ele(int basis, int target_ele_orbit);
   
};

#endif /* Model_1D_EKLM_hpp */
