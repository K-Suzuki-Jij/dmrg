#ifndef DMRG_hpp
#define DMRG_hpp

#include <string>
#include "SML.hpp"

struct DMRG_Param {
   
   std::string Initial_Guess;
   std::string Enviro_Copy;

   int max_dim_system;
   int now_sweep;
   int tot_sweep;
   int param_now_iter;
   int param_tot_iter;
   int renorm_now_iter;
   int renorm_tot_iter;
   
};

struct DMRG_Basis_LLLRRRRL {
  
   std::vector<short> LL, LR, RR, RL;
   std::vector<int> Inv;
   int base_LRRRRL, base_RRRL, base_RL;
   
   std::vector<int> QN;
   
   //Initial Guess
   std::vector<short> LL_Bef, LR_Bef, RR_Bef, RL_Bef;
   std::vector<int> Inv_Bef;
   
};

struct DMRG_Basis_LLLR {
  
   std::vector<short> LL, LR;
   std::vector<int> Inv, Count_Enviro;
   std::vector<std::vector<int>> QN_LLLR;
   
};

struct DMRG_Basis_LRRL {
  
   std::vector<short> LR, RL;
   std::vector<int> Inv;
   
};

struct DMRG_Basis_RRRL {
  
   std::vector<short> RR, RL;
   std::vector<int> Inv;
   
};

struct DMRG_Basis_LLRR {
   
   std::vector<short> LL, RR;
   std::vector<int> Inv;
   
};

struct DMRG_Basis_LLRL {
   
   std::vector<short> LL, RL;
   std::vector<int> Inv;
   
};

struct DMRG_Basis_LRRR {
   
   std::vector<short> LR, RR;
   std::vector<int> Inv;
   
};

struct DMRG_Basis_Stored {
  
   std::vector<std::vector<short>> LL_LLLR_Stored, LR_LLLR_Stored;
   std::vector<std::vector<int>>   Inv_LLLR_Stored;
   std::vector<std::vector<std::vector<int>>> QN_LL_LL_Stored;
   int qn_LL_LL_stored_ele_start, qn_LL_LL_stored_ele_end;
   int qn_LL_LL_stored_parity_start, qn_LL_LL_stored_parity_end;
   std::vector<CCS> Trans_Mat_Stored;
   
};

struct DMRG_Block_Information {
  
   int LL_site, RR_site;
   int dim_onsite, dim_LL, dim_RR;
   
};


struct DMRG_Ground_State {
   
   double val, error, tr_error, sum_tr_val;
   std::vector<double> Vector, Vector_Guess;
   
};

struct DMRG_Basis {
   
   DMRG_Basis_LLLRRRRL LLLRRRRL;
   DMRG_Basis_LLLR     LLLR;
   DMRG_Basis_LLRR     LLRR;
   DMRG_Basis_LLRL     LLRL;
   DMRG_Basis_LRRR     LRRR;
   DMRG_Basis_LRRL     LRRL;
   DMRG_Basis_RRRL     RRRL;

};

struct DMRG_Block_Hamiltonian {
  
   CRS LLLR, LLRR, LLRL, LRRR, LRRL, RRRL;
   std::vector<int> Ele_LL, Ele_RR, Ele_On;
   
};

struct DMRG_T_Mat {
   CCS Trans_Mat;
   //std::vector<int> QN1, QN2, QN3;
   std::vector<std::vector<int>> QN;
   std::vector<double> Eig_Val_DM;
};

struct DMRG_Onsite_Basis {
   int LL, LR, RR, RL;
   int LLLR, LLRR, LLRL, LRRR, LRRL, RRRL;
   int row;
};

struct DMRG_A_Basis_Set {
   int LL_site, RR_site;
   int dim_RR;
   int base_LRRRRL, base_RRRL, base_RL;
   double zero_precision;
   int elem_num;
   int c_obc;
   std::vector<int> Inv;
   std::vector<double> Val;
   std::vector<int> Check_Basis;
};

struct DMRG_Time {
   double diag;
   double inv_iter;
   double make_ham;
   double make_basis;
   double make_dm_mat;
   double total;
};

void DMRG_Clear_Check_Basis(DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv);
void DMRG_Construct_Superblock(DMRG_Basis &Basis, const DMRG_Basis_Stored &Basis_System, const DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, std::string Initial_Guess_Flag, DMRG_Time &Time, int p_threads);
void DMRG_Expectation_Intersite(std::vector<double> &Out, int site_ref, std::vector<CRS> &M_CF_LL, std::vector<CRS> &M_LL, CRS &M_On, std::vector<CRS> &M_RR, std::vector<double> &Vec, DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, int LL_site, int RR_site, int p_threads);
void DMRG_Expectation_Intersite(std::vector<double> &Out, int site_ref, std::vector<CRS> &M_CF_LL, std::vector<CRS> &M_LL, CRS &M_On, std::vector<CRS> &M_RR, std::vector<double> &Vec, DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, std::vector<DMRG_Basis_LLLRRRRL> &W_Basis_LLLRRRRL, int LL_site, int RR_site, int p_threads);
void DMRG_Expectation_Onsite(std::vector<double> &Out, const std::vector<CRS> &M_LL, const CRS &M_On, const std::vector<CRS> &M_RR, const std::vector<double> &Vec, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, int LL_site, int RR_site, int p_threads);
void DMRG_Extend_Mat_LL_LR_To_LLLR(const CRS &M_LL, const CRS &M_LR, CRS &Out_LLLR, const std::vector<int> &Ele_LL, std::string Sign_Flag, const DMRG_Basis_LLLR &Basis_LLLR);
int  DMRG_Find_Renormalized_Dim(std::vector<std::vector<double>> &Value_DM, std::vector<int> &Block_Sorted, std::vector<int> &Basis_Sorted, int max_dim_system);
void DMRG_Get_A_Basis_Set(std::vector<DMRG_A_Basis_Set> &A_Basis, int num, int p_threads);
void DMRG_Get_Density_Matrix(std::vector<std::vector<double>> &M, int block, const std::vector<int> &Dim_Block, const std::vector<double> &GS_Vec, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const DMRG_Basis_LLLR &Basis_LLLR);
void DMRG_Get_Trans_Matrix(DMRG_T_Mat &T_Mat, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const DMRG_Basis_LLLR &Basis_LLLR, DMRG_Ground_State &GS, int num_of_qn, int max_dim_system, int p_threads, DMRG_Time &Time);
void DMRG_Make_Elem_Onsite_LLLR_RRRL(const std::string Mat_Type, int basis_change, int basis_no_change, DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv, const CRS &M, double coeff);
void DMRG_Make_Elem_Onsite_LLLRRRRL(const std::string Mat_Type, int basis_change, const DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, const std::vector<short> &Basis_1, const std::vector<short> &Basis_2, const std::vector<int> &Inv_LLLRRRRL, const CRS &M, double coeff, const std::vector<int> &Ele_RR, const std::vector<int> &Ele_On, const std::string Sign_Flag);
void DMRG_Make_Elem_Intersite_Block(std::string Mat_Type, int basis_change_1, int basis_change_2, DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv, const CRS &M_1, const CRS &M_2, double coeff, const std::vector<int> &Ele_LL_or_LR_or_RR, std::string Sign_Flag);
void DMRG_Make_Elem_Zero_LLLRRRRL(const DMRG_Onsite_Basis &Basis, DMRG_A_Basis_Set &A_Basis, const std::vector<int> &Inv_LLLRRRRL);
void DMRG_Make_Elem_Ham_LLLRRRRL(const DMRG_Onsite_Basis &Basis_Onsite, DMRG_A_Basis_Set &A_Basis, const DMRG_Basis &Basis, const DMRG_Block_Hamiltonian &Ham, std::string Sign);
void DMRG_Get_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, const DMRG_Block_Hamiltonian &Block_Ham, const DMRG_Basis &Basis, const DMRG_Block_Information &Block, std::string Sign_Flag, int p_threads);
void DMRG_Matrix_Vector_Product(const std::string Mat_Type, const CRS &M, const std::vector<double> &V_In, std::vector<double> &V_Out, const std::vector<int> &Inv, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const int p_threads);
void DMRG_Output_Time(DMRG_Time &Time);
void DMRG_Print_Status(int num_of_qn, std::string BC, const DMRG_Ground_State &GS, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const DMRG_Basis_Stored &Basis_System, const DMRG_Basis_Stored &Basis_Enviro, const DMRG_Block_Information &Block,const  DMRG_Param &Dmrg_Param, const DMRG_Time &Time);
void DMRG_Renormalize_Matrix(std::string Mat_Type, const CRS &M, const CCS &Trans_Mat, CRS &M_Out_LL, const DMRG_Block_Information &Block, std::string Sign_Flag, const DMRG_Basis_LLLR &Basis_LLLR, const std::vector<int> &Ele_LL, int p_threads);
void DMRG_Renormalize_Matrix_LLLR(const CRS &M_LL, const CRS &M_LR, const CCS &Trans_Mat, CRS &M_Out_LL, const std::vector<int> &Ele_LL, std::string Sign_Flag, const DMRG_Basis_LLLR &Basis_LLLR);
void DMRG_Sort_Bases(DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, std::vector<std::vector<int>> &QN_LLLR_LLLRRRRL, int left, int right);
void DMRG_Sort_Density_Matrix_Bases(std::vector<std::vector<double>> &Value_DM, std::vector<int> &Block_Sorted, std::vector<int> &Basis_Sorted, std::vector<std::vector<int>> &QN_Block_Sorted, long left, long right);
void DMRG_Store_Trans_Matrix(CCS &Trans_Mat, DMRG_Basis_Stored &Basis_System, DMRG_Basis_LLLR &Basis_LLLR, DMRG_Param &Dmrg_Param, int LL_site, int system_size);
void DMRG_Transform_Matrix_One(const CRS &M_On, std::vector<CRS> &M_Out, const DMRG_Basis_Stored &Basis_System, const DMRG_Block_Information &Block, int LL_site_now, int p_threads);
void DMRG_Transform_Matrix_Two(const CRS &M_On, std::vector<CRS> &M_Out, const DMRG_Basis_Stored &Basis_System, const DMRG_Block_Information &Block, std::vector<int> &Ele_LL, std::string Sign_Flag, int LL_site_now, int cf_origin, int p_threads);
void DMRG_Initial_Guess(DMRG_Ground_State &GS, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const DMRG_Basis_LLLR &Basis_LLLR, const DMRG_Basis_Stored &Basis_System, const DMRG_Basis_Stored &Basis_Enviro, const DMRG_Block_Information &Block, std::string Initial_Guess_Flag, int p_threads);
void DMRG_Get_Inv_LLLRRRRL(DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const DMRG_Block_Information &Block, int p_threads);
void DMRG_Diagonalize_Ham_LLLRRRRL(CRS &Ham_LLLRRRRL, DMRG_Ground_State &GS, DMRG_Basis &Basis, DMRG_Block_Information &Block, DMRG_Param &Dmrg_Param, Diag_Param &Diag_Param, DMRG_Time &Time, int p_threads);
#endif /* DMRG_hpp */
