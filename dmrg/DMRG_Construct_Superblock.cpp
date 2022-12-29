#include <iostream>
#include <vector>
#include <omp.h>
#include <climits>
#include "DMRG.hpp"

void DMRG_Construct_Superblock(DMRG_Basis &Basis, const DMRG_Basis_Stored &Basis_System, const DMRG_Basis_Stored &Basis_Enviro, DMRG_Block_Information &Block, std::string Initial_Guess_Flag, DMRG_Time &Time, int p_threads) {
   
   double start = omp_get_wtime();
   
   int LL_site    = Block.LL_site;
   int RR_site    = Block.RR_site;
   int dim_LL     = (int)Basis_System.QN_LL_LL_Stored[LL_site][0].size();
   int dim_RR     = (int)Basis_Enviro.QN_LL_LL_Stored[RR_site][0].size();
   int dim_onsite = (int)Basis_System.QN_LL_LL_Stored[   0   ][0].size();
   int num_of_qn  = (int)Basis_System.QN_LL_LL_Stored[0].size();
   
   if (dim_LL <= 0 || dim_RR <= 0 || dim_onsite <= 0) {
      std::cout << "Error in DMRG_Construct_Superblock at 1" << std::endl;
      std::exit(0);
   }
   
   Block.dim_LL     = dim_LL;
   Block.dim_RR     = dim_RR;
   Block.dim_onsite = dim_onsite;
   
   std::vector<int> Target_QN(num_of_qn);
   
   for (int qn = 0; qn < num_of_qn; qn++) {
      Target_QN[qn] = Basis.LLLRRRRL.QN[qn];
   }
   
   //For Initial Guess of GS Vector
   if (Initial_Guess_Flag == "Yes") {
      Basis.LLLRRRRL.LL_Bef  = Basis.LLLRRRRL.LL;
      Basis.LLLRRRRL.LR_Bef  = Basis.LLLRRRRL.LR;
      Basis.LLLRRRRL.RR_Bef  = Basis.LLLRRRRL.RR;
      Basis.LLLRRRRL.RL_Bef  = Basis.LLLRRRRL.RL;
      Basis.LLLRRRRL.Inv_Bef = Basis.LLLRRRRL.Inv;
   }
   
   std::vector<std::vector<int>> QN_LLLR_LLLRRRRL(num_of_qn);
   
   Basis.LLLRRRRL.LL.clear();
   Basis.LLLRRRRL.LR.clear();
   Basis.LLLRRRRL.RR.clear();
   Basis.LLLRRRRL.RL.clear();
   
   long count = 0;
   
#pragma omp parallel for num_threads (p_threads)
   for (int LL = 0; LL < dim_LL; LL++) {
      for (int LR = 0; LR < dim_onsite; LR++) {
         for (int RR = 0; RR < dim_RR; RR++) {
            for (int RL = 0; RL < dim_onsite; RL++) {
                              
               int flag = 0;
               for (int i = 0; i < num_of_qn; i++) {
                  int qn = Basis_System.QN_LL_LL_Stored[LL_site][i][LL] +
                  Basis_System.QN_LL_LL_Stored[   0   ][i][LR] +
                  Basis_Enviro.QN_LL_LL_Stored[RR_site][i][RR] +
                  Basis_Enviro.QN_LL_LL_Stored[   0   ][i][RL];
                  
                  if (Basis_System.qn_LL_LL_stored_parity_start <= i && i < Basis_System.qn_LL_LL_stored_parity_end) {
                     qn = qn%2;
                  }
                  
                  if (qn != Target_QN[i]) {
                     flag = 1;
                     break;
                  }
               }
               
               if (flag == 0) {
#pragma omp critical
                  {
                     Basis.LLLRRRRL.LL.push_back(LL);
                     Basis.LLLRRRRL.LR.push_back(LR);
                     Basis.LLLRRRRL.RR.push_back(RR);
                     Basis.LLLRRRRL.RL.push_back(RL);
                     for (int i = 0; i < num_of_qn; i++) {
                        int qn = Basis_System.QN_LL_LL_Stored[LL_site][i][LL] + Basis_System.QN_LL_LL_Stored[0][i][LR];
                        if (Basis_System.qn_LL_LL_stored_parity_start <= i && i < Basis_System.qn_LL_LL_stored_parity_end) {
                           qn = qn%2;
                        }
                        QN_LLLR_LLLRRRRL[i].push_back(qn);
                     }
                     count++;
                  }
               }
            }
         }
      }
   }
      
   if (count != (int)Basis.LLLRRRRL.LL.size() || INT_MAX <= count || count <= 0) {
      std::cout << "Error in DMRG_Construct_Superblock at 2" << std::endl;
      std::cout << "dim_LLLRRRRL=" << count << std::endl;
      std::exit(0);
   }
   
   int dim_LLLRRRRL = (int)Basis.LLLRRRRL.LL.size();
   DMRG_Sort_Bases(Basis.LLLRRRRL, QN_LLLR_LLLRRRRL, 0, dim_LLLRRRRL);
   
   DMRG_Get_Inv_LLLRRRRL(Basis.LLLRRRRL, Block, p_threads);
   
   Basis.LLLR.Inv.resize(dim_LL*dim_onsite);
   Basis.LLRL.Inv.resize(dim_LL*dim_onsite);
   Basis.LLRR.Inv.resize(dim_LL*dim_RR);
   Basis.LRRR.Inv.resize(dim_RR*dim_onsite);
   Basis.LRRL.Inv.resize(dim_onsite*dim_onsite);
   Basis.RRRL.Inv.resize(dim_RR*dim_onsite);
   
#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < dim_LL*dim_onsite; i++) {
      Basis.LLLR.Inv[i] = -1;
      Basis.LLRL.Inv[i] = -1;
   }

#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < dim_RR*dim_onsite; i++) {
      Basis.LRRR.Inv[i] = -1;
      Basis.RRRL.Inv[i] = -1;
   }
   
#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < dim_LL*dim_RR; i++) {
      Basis.LLRR.Inv[i] = -1;
   }
   
#pragma omp parallel for num_threads (p_threads)
   for (int i = 0; i < dim_onsite*dim_onsite; i++) {
      Basis.LRRL.Inv[i] = -1;
   }
   
   Basis.LLLR.LL.clear();
   Basis.LLLR.LR.clear();
   Basis.LLLR.Count_Enviro.clear();
   
   Basis.LLLR.QN_LLLR.resize(num_of_qn);
   for (int i = 0; i < num_of_qn; i++) {
      Basis.LLLR.QN_LLLR[i].clear();
   }
   
   Basis.LLRR.LL.clear();
   Basis.LLRR.RR.clear();
   
   Basis.LLRL.LL.clear();
   Basis.LLRL.RL.clear();

   Basis.LRRR.LR.clear();
   Basis.LRRR.RR.clear();
   
   Basis.LRRL.LR.clear();
   Basis.LRRL.RL.clear();

   Basis.RRRL.RR.clear();
   Basis.RRRL.RL.clear();
   
   int count_lllr = 0;
   int count_llrr = 0;
   int count_llrl = 0;
   int count_lrrr = 0;
   int count_lrrl = 0;
   int count_rrrl = 0;
   for (int i = 0; i < dim_LLLRRRRL; i++) {
      int LL   = Basis.LLLRRRRL.LL[i];
      int LR   = Basis.LLLRRRRL.LR[i];
      int RR   = Basis.LLLRRRRL.RR[i];
      int RL   = Basis.LLLRRRRL.RL[i];
      int LLLR = LL*dim_onsite + LR;
      int LLRR = LL*dim_RR     + RR;
      int LLRL = LL*dim_onsite + RL;
      int LRRR = LR*dim_RR     + RR;
      int LRRL = LR*dim_onsite + RL;
      int RRRL = RR*dim_onsite + RL;
      
      if (Basis.LLLR.Inv[LLLR] == -1) {
         
         for (int i = 0; i < num_of_qn; i++) {
            int qn = Basis_System.QN_LL_LL_Stored[LL_site][i][LL] + Basis_System.QN_LL_LL_Stored[0][i][LR];
            if (Basis_System.qn_LL_LL_stored_parity_start <= i && i < Basis_System.qn_LL_LL_stored_parity_end) {
               qn = qn%2;
            }
            Basis.LLLR.QN_LLLR[i].push_back(qn);
         }
         
         Basis.LLLR.LL.push_back(LL);
         Basis.LLLR.LR.push_back(LR);
         Basis.LLLR.Count_Enviro.push_back(1);
         Basis.LLLR.Inv[LLLR] = count_lllr;
         
         count_lllr++;
      }
      else {
         Basis.LLLR.Count_Enviro[Basis.LLLR.Inv[LLLR]]++;
      }
      
      if (Basis.LLRR.Inv[LLRR] == -1) {
         Basis.LLRR.LL.push_back(LL);
         Basis.LLRR.RR.push_back(RR);
         Basis.LLRR.Inv[LLRR] = count_llrr;
         count_llrr++;
      }
      
      if (Basis.LLRL.Inv[LLRL] == -1) {
         Basis.LLRL.LL.push_back(LL);
         Basis.LLRL.RL.push_back(RL);
         Basis.LLRL.Inv[LLRL] = count_llrl;
         count_llrl++;
      }
      
      if (Basis.LRRR.Inv[LRRR] == -1) {
         Basis.LRRR.LR.push_back(LR);
         Basis.LRRR.RR.push_back(RR);
         Basis.LRRR.Inv[LRRR] = count_lrrr;
         count_lrrr++;
      }
      
      if (Basis.LRRL.Inv[LRRL] == -1) {
         Basis.LRRL.LR.push_back(LR);
         Basis.LRRL.RL.push_back(RL);
         Basis.LRRL.Inv[LRRL] = count_lrrl;
         count_lrrl++;
      }
      
      if (Basis.RRRL.Inv[RRRL] == -1) {
         Basis.RRRL.RR.push_back(RR);
         Basis.RRRL.RL.push_back(RL);
         Basis.RRRL.Inv[RRRL] = count_rrrl;
         count_rrrl++;
      }
      
   }
   
   if (count_lllr != (int)Basis.LLLR.LL.size()) {
      std::cout << "Error in DMRG_Construct_Superblock at 4" << std::endl;
      std::exit(0);
   }
   if (count_rrrl != (int)Basis.RRRL.RR.size()) {
      std::cout << "Error in DMRG_Construct_Superblock at 5" << std::endl;
      std::exit(0);
   }
   if (count_lrrl != (int)Basis.LRRL.LR.size()) {
      std::cout << "Error in DMRG_Construct_Superblock at 6" << std::endl;
      std::exit(0);
   }
   if (count_llrr != (int)Basis.LLRR.LL.size()) {
      std::cout << "Error in DMRG_Construct_Superblock at 7" << std::endl;
      std::exit(0);
   }
   if (count_llrl != (int)Basis.LLRL.LL.size()) {
      std::cout << "Error in DMRG_Construct_Superblock at 8" << std::endl;
      std::exit(0);
   }
   if (count_lrrr != (int)Basis.LRRR.LR.size()) {
      std::cout << "Error in DMRG_Construct_Superblock at 9" << std::endl;
      std::exit(0);
   }
   
   Time.make_basis = omp_get_wtime() - start;

}
