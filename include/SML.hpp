#ifndef SML_hpp
#define SML_hpp

#include <vector>
#include <string>

struct Diag_Param {
   
   int diag_num;
   
   //Lanczos or LOBPCG
   std::string Diag_Method;
   std::string Calc_Vec;
   std::string Lanczos_Initial_Guess;
   double diag_acc;
   int diag_min_step;
   int diag_max_step;
   
   //Conjugate gradient Method or Minimal Residual
   double cg_acc;
   int cg_max_step;

   //Inverse iteration
   std::string ii_Method;
   double ii_acc, ii_diag_add;
   int ii_max_step;
   
   //Matrix Type
   std::string Mat_Type;
   
};

struct CRS {
   
   int row_dim;
   int col_dim;
   
   std::vector<int>    Col;
   std::vector<long>   Row;
   std::vector<double> Val;
   
};

struct CCS {
   
   int row_dim;
   int col_dim;
   
   std::vector<long>   Col;
   std::vector<int>    Row;
   std::vector<double> Val;
   
};

template <typename T> long Binary_Search(const std::vector<T> &Array, long min, long max, T target_val);
long   Binomial_Coefficient(long n, long k);
void   Change_Matrix_Basis(const CRS &M, const CCS &T, CRS &Out, int p_threads);
void   Check_Symmetric_Matrix(const CRS &M, double zero_precision);
template <class SPM> void Clear_Matrix(SPM &M);
void   Conjugate_Gradient(const CRS &M, const std::vector<double> &Vec_In, std::vector<double> &Vec_Out, const Diag_Param &Diag_Param, int p_threads);
void   Copy_Vector(const std::vector<double> &V1, std::vector<double> &Copyed_V, int p_threads);
template <typename T> int Delta_Function(T i, T j);
void   Diag_Add_Matrix(CRS &M, double add, int p_threads);
void   Find_Nth_Permutation(std::vector<int> &Array, long target_number);
void   Free_CRS(CRS &M);
double Inner_Product(const std::vector<double> &V1, const std::vector<double> &V2, int p_threads);
void   Inverse_Iteration(CRS &M, std::vector<double> &Eigen_Vec, double eigen_val, const Diag_Param &Diag_Param, int p_threads);
double L2_Norm(const std::vector<double> &V, int p_threads);
void   Lanczos_Ex1(const CRS &M, const std::vector<double> &Eigen_Vec, double &out_val, const Diag_Param &Diag_Param, int p_threads);
void   Lanczos(const CRS &M, std::vector<double> &Out_Vec, double &out_val, const Diag_Param &Diag_Param, int p_threads);
void   Lapack_Dstev(int dim, const std::vector<double> &Diag, const std::vector<double> &Off_Diag, std::vector<double> &Vec, double &val);
void   Lapack_Dsyev(const CRS &M, std::vector<double> &Vec, double &val, int state);
void   Lapack_Dsyev(const std::vector<std::vector<double>> &M, std::vector<std::vector<double>> &Vec, std::vector<double> &Val);
void   Matrix_Constant_Multiplication(CRS &M, double coeef, int p_threads);
void   Matrix_Matrix_Product(const CRS &M1, const CRS &M2, CRS &Out);
void   Matrix_Matrix_Sum(const CRS &M1, const CRS &M2, CRS &Out);
void   Matrix_Transpose(const CRS &M, CRS &M_Out);
void   Matrix_Transpose(const CCS &M, CCS &M_Out);
void   Matrix_Vector_Product(const CRS &M, const std::vector<double> &V, std::vector<double> &V_Out,  std::vector<std::vector<double>> &Temp_V, std::string Mat_Type, int p_threads);
long   Multinomial_Coefficient(std::vector<int> &Array);
void   Normalize(std::vector<double> &V, int p_threads);
void   Output_Step_Number(int step_num, double time, std::string file_name);
void   Partition_Integer(int div_num, int max_num, std::vector<std::vector<int>> &Out_Array);
void   Print_CCS(const CCS &M, std::string Name);
void   Print_CRS(const CRS &M, std::string Name);
template <typename T> void Print_Vector(std::vector<T> &Vec, std::string Name);
template <typename T> void Print_Vector(std::vector<std::vector<T>> &Vec, std::string Name);
template <typename T_Target, typename T_1> void Quick_Sort_Vector(std::vector<T_Target> &Array_Target, std::vector<T_1> &Array1, long left, long right);
template <typename T_Target, typename T_1, typename T_2> void Quick_Sort_Vector(std::vector<T_Target> &Array_Target, std::vector<T_1> &Array1, std::vector<T_2> &Array2, long left, long right);
double Residual_Error_Eigenpair(const CRS &M, const std::vector<double> &Eigen_Vec, double eigen_val, std::string Mat_Type, int p_threads);
void   Sort_Col_CRS(CRS &M, int p_threads);
void   Sort_Col_CCS(CCS &M, int p_threads);
double Vector_Matrix_Vector_Product(std::vector<double> &V1, CRS &M, std::vector<double> &V2, std::string Mat_Type, int p_threads);

#endif /* SML_hpp */
