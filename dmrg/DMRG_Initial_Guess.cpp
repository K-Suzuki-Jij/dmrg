#include <iostream>
#include "DMRG.hpp"

void DMRG_Initial_Guess(DMRG_Ground_State &GS, const DMRG_Basis_LLLRRRRL &Basis_LLLRRRRL, const DMRG_Basis_LLLR &Basis_LLLR, const DMRG_Basis_Stored &Basis_System, const DMRG_Basis_Stored &Basis_Enviro, const DMRG_Block_Information &Block, std::string Initial_Guess_Flag, int p_threads) {
   
   int LL_site = Block.LL_site;
   int RR_site = Block.RR_site;
   
   if (LL_site <= 0 || LL_site == RR_site || Initial_Guess_Flag != "Yes") {
      GS.Vector_Guess.resize(0);
      return;
   }
   
   int LL_site_bef      = LL_site - 1;
   int dim_onsite       = Block.dim_onsite;
   int dim_LL           = Block.dim_LL;
   int dim_RR           = Block.dim_RR;
   int dim_RR_bef       = (int)Basis_Enviro.QN_LL_LL_Stored[RR_site + 1][0].size();
   int dim_LLLR         = (int)Basis_LLLR.LL.size();
   int dim_LLLRRRRL     = (int)Basis_LLLRRRRL.LL.size();
   int dim_LLLRRRRL_bef = (int)Basis_LLLRRRRL.LL_Bef.size();
   int base_LRRRRL      = Basis_LLLRRRRL.base_LRRRRL;
   int base_RRRL        = Basis_LLLRRRRL.base_RRRL;
   int base_RL          = Basis_LLLRRRRL.base_RL;
   int base_LRRRRL_bef  = dim_onsite*dim_RR_bef*dim_onsite;
   int base_RRRL_bef    = dim_RR_bef*dim_onsite;
   
   if (dim_LL != Basis_System.Trans_Mat_Stored[LL_site_bef].col_dim) {
      std::cout << "Error in DMRG_Initial_Guess at 2" << std::endl;
      std::exit(0);
   }
   if (dim_RR_bef != Basis_Enviro.Trans_Mat_Stored[RR_site].col_dim) {
      std::cout << "Error in DMRG_Initial_Guess at 3" << std::endl;
      std::exit(0);
   }
   
   std::vector<int> Ele_LR(dim_onsite, 0), Ele_RR(dim_RR, 0);
   
   for (int i = Basis_System.qn_LL_LL_stored_ele_start; i < Basis_System.qn_LL_LL_stored_ele_end; i++) {
      for (int j = 0; j < dim_onsite; j++) {
         Ele_LR[j] += Basis_System.QN_LL_LL_Stored[0][i][j];
      }
      for (int j = 0; j < dim_RR; j++) {
         Ele_RR[j] += Basis_System.QN_LL_LL_Stored[RR_site][i][j];
      }
   }
   
   std::vector<std::vector<double>> Temp_V(dim_LLLR, std::vector<double>(dim_RR_bef, 0.0));

   CCS Trans;
   Matrix_Transpose(Basis_System.Trans_Mat_Stored[LL_site_bef], Trans);

   //Fix me: Use OpenMp
   for (int i = 0; i < dim_LLLRRRRL_bef; i++) {
      int LL_bef       = Basis_LLLRRRRL.LL_Bef[i];
      int LR_bef       = Basis_LLLRRRRL.LR_Bef[i];
      int RR_bef       = Basis_LLLRRRRL.RR_Bef[i];
      int LR_new       = Basis_LLLRRRRL.RL_Bef[i];//LR_new = RL_bef
      int inv_LLLR_bef = Basis_System.Inv_LLLR_Stored[LL_site_bef][LL_bef*dim_onsite + LR_bef];
      long LLLRRRRL    = (long)LL_bef*base_LRRRRL_bef + LR_bef*base_RRRL_bef + RR_bef*base_RL + LR_new;
      int  inv_bef     = Basis_LLLRRRRL.Inv_Bef[LLLRRRRL];
      if (inv_LLLR_bef >= 0 && inv_bef >= 0) {
         double val = GS.Vector[inv_bef];
         for (long j = Trans.Col[inv_LLLR_bef]; j < Trans.Col[inv_LLLR_bef + 1]; j++) {
            int inv_LLLR_new = Basis_LLLR.Inv[Trans.Row[j]*dim_onsite + LR_new];
            if (inv_LLLR_new >= 0) {
               Temp_V[inv_LLLR_new][RR_bef] += val*Trans.Val[j];
            }
         }
      }
   }
   
   GS.Vector_Guess.resize(dim_LLLRRRRL);

   Matrix_Transpose(Basis_Enviro.Trans_Mat_Stored[RR_site], Trans);
   
#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      int LL_new       = Basis_LLLRRRRL.LL[i];
      int LR_new       = Basis_LLLRRRRL.LR[i];
      int RR_new       = Basis_LLLRRRRL.RR[i];
      int RL_new       = Basis_LLLRRRRL.RL[i];
      int inv_LLLR_new = Basis_LLLR.Inv[LL_new*dim_onsite + LR_new];
      int inv_RRRL_new = Basis_Enviro.Inv_LLLR_Stored[RR_site][RR_new*dim_onsite + RL_new];
      double val = 0.0;
            
      int sign = 1;
      
      if (Ele_LR[LR_new]%2 == 1 && (Ele_RR[RR_new] + Ele_LR[RL_new])%2 == 0) {
         sign = -1;
      }
      
      if (inv_LLLR_new >= 0 && inv_RRRL_new >= 0) {
         for (long j = Trans.Col[inv_RRRL_new]; j < Trans.Col[inv_RRRL_new + 1]; j++) {
            val += Temp_V[inv_LLLR_new][Trans.Row[j]]*Trans.Val[j]*sign;
         }
      }
      
      int inv_new = Basis_LLLRRRRL.Inv[(long)LL_new*base_LRRRRL + LR_new*base_RRRL + RR_new*base_RL + RL_new];
      
      if (inv_new < 0) {
         std::cout << "Error in DMRG_Initial_Guess at 3" << std::endl;
         std::exit(0);
      }
      GS.Vector_Guess[inv_new] = val;
   }

   Normalize(GS.Vector_Guess, p_threads);

}
