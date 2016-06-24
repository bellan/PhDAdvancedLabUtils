#include <iostream>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <TF1.h>
#include <TEllipse.h>

using namespace std;

void vertex(){
  cout << "ciao" << endl;

  //TMatrixD x0(2,1);
  
  TMatrixD r_1(3,1);
  r_1[0][0] = 1; r_1[1][0] = 2; r_1[2][0] = 1;
  TMatrixD cov_1(3,3);
  cov_1[0][0] = 0.1; cov_1[1][1] = 0.1; cov_1[2][2] = 0.01;
  

  TMatrixD r_2(3,1);
  r_2[0][0] = 4; r_2[1][0] = 2; r_2[2][0] = -1;
  TMatrixD cov_2(3,3);
  cov_2[0][0] = 0.1; cov_2[1][1] = 0.1; cov_2[2][2] = 0.01;


  vector<TMatrixD> measures;
  vector<TMatrixD> cov;

  measures.push_back(r_1); cov.push_back(cov_1);
  measures.push_back(r_2); cov.push_back(cov_2);
  
 
  // -------------------------------------------------------------------------
  


  // Derived quantitities
  TMatrixD C_n(2,2);
  TMatrixD x_v(2,1);
 
  vector<TMatrixD> A; vector<TMatrixD> W; vector<TMatrixD> B_T; vector<TMatrixD> G; vector<TMatrixD> ce;
  vector<TMatrixD> B;

  vector<TF1*> measureGraphs;
  for(uint k = 0; k < measures.size(); ++k){
    measureGraphs.push_back(new TF1("a","pol1",-10,10));
    measureGraphs[k]->SetParameter(0,measures[k][1][0]-measures[k][2][0]*measures[k][0][0]);
    measureGraphs[k]->SetParameter(1,measures[k][2][0]);
    measureGraphs[k]->Draw(k == 0 ? "" : "same");

    TMatrixD G_k = cov[k];
    G_k.Invert();
    G.push_back(G_k);


    TMatrixD ce_k(3,1);
    ce_k[1][0] = measures[k][1][0];
    ce.push_back(ce_k);

    TMatrixD A_k(3,2);
    A_k[0][0] = 1;
    A_k[0][1] = -measures[k][2][0];
    A.push_back(A_k);

    TMatrixD A_k_T = A_k;
    A_k_T.Transpose(A_k_T);

    TMatrixD B_k(3,1);
    B_k[0][0] =  measures[k][1][0];
    B_k[2][0] =  1;
    B.push_back(B_k);

    TMatrixD B_k_T = B_k;
    B_k_T.Transpose(B_k_T);
    B_T.push_back(B_k_T);

    TMatrixD W_k = (B_k_T*G_k*B_k).Invert();
    W.push_back(W_k);
    
    TMatrixD G_k_B = G_k - G_k*B_k*W_k*B_k_T*G_k;

    C_n += A_k_T*G_k_B*A_k;
    x_v += A_k_T*G_k_B*(measures[k]-ce_k);
  }

  C_n.Invert();
  x_v = C_n * x_v;

  x_v.Print();
  C_n.Print();

  // TODO: add correlatons
  TEllipse *vpos = new TEllipse(x_v[0][0],x_v[1][0],C_n[0][0],C_n[1][1]);
  vpos->Draw();

  // Compute the new direction and the full covariance matrix
  vector<TMatrixD> newdirections;
  TMatrixD chi2(1,1);
  for(uint k = 0; k < measures.size(); ++k){

    newdirections.push_back(W[k]*B_T[k]*G[k]*(measures[k]-ce[k]-A[k]*x_v));
    newdirections[k].Print();

    // compute full chi2
    TMatrixD p_n_k = ce[k] + A[k]*x_v + B[k]*newdirections[k];
    TMatrixD r_n_k = measures[k] - p_n_k;
    TMatrixD r_n_k_T = r_n_k;
    r_n_k_T.Transpose(r_n_k_T);

    chi2 += r_n_k_T*G[k]*r_n_k;

  }
  chi2.Print();
  

}
