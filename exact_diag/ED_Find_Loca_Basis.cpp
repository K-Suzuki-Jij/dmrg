
int ED_Find_Local_Basis(long basis, int site, int dim_onsite) {
     
     for (int dummy_i = 0; dummy_i < site; dummy_i++) {
        basis = basis/dim_onsite;
     }
     basis = basis%dim_onsite;
     
     return (int)basis;
  }
