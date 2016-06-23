#include <iostream>
#include <TMatrixD.h>
#include <TVectorD.h>

using namespace std;

void vertex(){
  cout << "ciao" << endl;

  vector<TMatrixD> measures;
  vector<TMatrixD> cov;

  TMatrixD r_1(3,1);
  r_1[0][0] = 1; r_1[1][0] = 2; r_1[2][0] = 1;
  TMatrixD cov_1(3,3);
  cov_1[0][0] = 0.1; cov_1[1][1] = 0.1; cov_1[2][2] = 0.01;
  

  TMatrixD r_2(3,1);
  r_2[0][0] = 4; r_2[1][0] = 2; r_2[2][0] = -1;
  TMatrixD cov_2(3,3);
  cov_2[0][0] = 0.1; cov_2[1][1] = 0.1; cov_2[2][2] = 0.01;

  measures.push_back(r_1); cov.push_back(cov_1);
  measures.push_back(r_2); cov.push_back(cov_2);
  
  // Derived quantitities
  TMatrixD C_n(2,2);
  TMatrixD x_v(2,1);
 
  for(int k = 0; k < measures.size(); ++k){
    TMatrixD G_k = cov[k];
    G_k.Invert();


    TMatrixD ce_k(3,1);
    ce_k[1][0] = measures[k][1][0];

    TMatrixD A_k(3,2);
    A_k[0][0] = 1;
    A_k[0][1] = -measures[k][2][0];
    
    TMatrixD A_k_T = A_k;
    A_k_T.Transpose(A_k_T);

    TMatrixD B_k(3,1);
    B_k[0][0] =  measures[k][1][0];
    B_k[2][0] =  1;

    TMatrixD B_k_T = B_k;
    B_k_T.Transpose(B_k_T);
    

    TMatrixD W_k = (B_k_T*G_k*B_k).Invert();

    TMatrixD G_k_B = G_k - G_k*B_k*W_k*B_k_T*G_k;

    C_n += A_k_T*G_k_B*A_k;
    x_v += A_k_T*G_k_B*(measures[k]-ce_k);
  }

  C_n.Invert();
  x_v = C_n * x_v;

  x_v.Print();
 

}
