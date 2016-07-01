#include <iostream>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <TF1.h>
#include <TEllipse.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TString.h>
#include <sstream>

using namespace std;

class KFOperations{
public:
  KFOperations(){}

  void filter(const TMatrixD& measure, const TMatrixD& covariance){
    
    G.ResizeTo(covariance);
    G = covariance;
    G.Invert();

    ce.ResizeTo(3,1);
    ce[1][0] = measure[1][0];

    A.ResizeTo(3,2);
    A[0][0] = 1;
    A[0][1] = -measure[2][0];

    A_T.ResizeTo(A);
    A_T = A;
    A_T.Transpose(A_T);

    B.ResizeTo(3,1);
    B[0][0] =  measure[1][0];
    B[2][0] =  1;

    B_T.ResizeTo(B);
    B_T = B;
    B_T.Transpose(B_T);

    TMatrixD Wtemp = (B_T*G*B).Invert();
    W.ResizeTo(Wtemp);
    W = Wtemp;
    
    TMatrixD G_Btemp = G - G*B*W*B_T*G;
    G_B.ResizeTo(G_Btemp);
    G_B = G_Btemp;

    C_vtx_.ResizeTo(2,2);
    C_vtx_ = A_T*G_B*A;
    x_vtx_.ResizeTo(2,1);
    x_vtx_ = A_T*G_B*(measure-ce);


  }


  void smooth(const TMatrixD &measure, const TMatrixD & x_v){
    TMatrixD newDir = W * B_T * G * (measure-ce-A * x_v);
    // compute full chi2
    TMatrixD p_n = ce + A*x_v + B* newDir;
    TMatrixD r_n = measure - p_n;
    TMatrixD r_n_T = r_n;
    r_n_T.Transpose(r_n_T);

    chi2_ = (r_n_T * G * r_n)[0][0];
    
  } 

  void inverseFilter(const TMatrixD &measure, const TMatrixD &x_v, const TMatrixD &C_n){
    TMatrixD C_n_1(TMatrixD::kInverted,C_n);
    TMatrixD C_n_star_1 = C_n_1 - A_T*G_B*A;
    TMatrixD C_n_star(TMatrixD::kInverted,C_n_star_1);
    TMatrixD x_v_star = C_n_star*(C_n_1*x_v - A_T*G_B*(measure-ce));
    TMatrixD r = x_v-x_v_star;
    TMatrixD r_T(TMatrixD::kTransposed,r);
    // cout << "-- chi2 " << chi2_ << endl;
    // cout << (r_T*C_n_star_1*r)[0][0] << endl;

    chi2_S_ = (r_T*C_n_star_1*r)[0][0] + chi2_ ;
    // cout << "-- chi2 " << chi2_S_ << endl;
  }


  double chi2() const {return chi2_;}
  double chi2_S() const {return chi2_S_;}

  TMatrixD C_vtx() const {return C_vtx_;}
  TMatrixD x_vtx() const {return x_vtx_;} 

  private:
  // Matrices for KF
  TMatrixD G;// = measures[k].covariance();
  TMatrixD ce;//(3,1);
  TMatrixD A;//(3,2);
  TMatrixD A_T;// = A_k;
  TMatrixD B;//(3,1);
  TMatrixD B_T;// = B_k;
  TMatrixD W;
  TMatrixD G_B;
  TMatrixD C_vtx_;
  TMatrixD x_vtx_; 

  double chi2_;
  double chi2_S_;

};


class Measurement{
public:
  // Be careful, x_k[2][0] is dx/dy, so when translating into a straight line in the (x,y) plane m = 1/x_k[2][0]
  Measurement(const double& x, const double& y, const double& dxdy, 
	      const double& sx2, const double& sy2, const double& sdxdy2, 
	      const double& cxy, const double& cxdxdy, const double& cydxdy){
  

    measure_.ResizeTo(3,1);
    covariance_.ResizeTo(3,3);

    measure_[0][0] = x; measure_[1][0] = y; measure_[2][0] = dxdy;
    covariance_[0][0] = sx2; covariance_[1][1] = sy2;    covariance_[2][2] = sdxdy2;
    covariance_[0][1] = cxy; covariance_[0][2] = cxdxdy; covariance_[1][2] = cydxdy;

    func_ = new TF1("a","pol1",-100,100);
    func_->SetParameter(0,measure()[1][0]-measure()[0][0]/measure()[2][0]);
    func_->SetParameter(1,1./measure()[2][0]);

    // do the KF k-step
    KFOperations_ = new KFOperations();

  }

  static Measurement convert(const double& q, const double m, const double& sq, const double &sm, const double &y_ref, const double &sy_ref){
    
    return Measurement((y_ref-q)/m, y_ref, 1./m, 
		       sqrt(pow(sy_ref/m,2)+pow(sq/m,2)+pow(sm*(y_ref-q)/(m*m),2)), sy_ref, sm, 
		       0., 0., 0.);
    
    
    
  }

  //  Measurement(const Measurement &meas){
  //  (*this).measure_      = meas.measure();
  //  (*this).covariance_   = meas.covariance();
  //  //(*this).func_         = (TF1*)meas.function()->Clone();
  //  (*this).KFOperations_ = meas.KF();
  //}

  
  TMatrixD measure() const {return measure_;}
  TMatrixD covariance() const {return covariance_;}

  TF1 *function()  {return func_;}
  KFOperations *KF() const {return KFOperations_;}
  void filter()                                               {KF()->filter(measure(),covariance());}
  void smooth(const TMatrixD& x_v)                            {KF()->smooth(measure(),x_v);}
  void inverseFilter(const TMatrixD &x_v, const TMatrixD &C_n){KF()->inverseFilter(measure(), x_v, C_n);}


private:
  TMatrixD measure_;
  TMatrixD covariance_;
  TF1 *func_;
  KFOperations *KFOperations_;
};


pair<double, int> searchVertex(vector<Measurement>& measures, TCanvas *canvas){
  
  TMatrixD C_n(2,2);
  TMatrixD x_v(2,1);
  

  canvas->cd();

  for(uint k = 0; k < measures.size(); ++k){
    measures[k].function()->Draw(k == 0 ? "" : "same");
    //measures[k].function()->Draw("same");
    measures[k].function()->GetYaxis()->SetRangeUser(-100,200);
    //measures[k].function()->SetLineColor(col);

    measures[k].filter();
    
    C_n += measures[k].KF()->C_vtx();
    x_v += measures[k].KF()->x_vtx();
  }
  
  C_n.Invert();
  x_v = C_n * x_v;
  
  cout << "Vertex position and its uncertainties is" << endl; 
  x_v.Print();
  C_n.Print();
  cout << endl;
  
  // TODO: add correlatons
  //TVectorD eigenValues;
  //TMatrixD eigenVectors = C_n.EigenVectors(eigenValues); 
  //eigenValues.Print();
  //eigenVectors.Print();
  
  TEllipse *vpos = new TEllipse(x_v[0][0],x_v[1][0],C_n[0][0],C_n[1][1],0,360);
  vpos->Draw();
  
  // Compute the new direction and the full covariance matrix
  double chi2 = 0.;
  int ndof = 0;
  
  for(uint k = 0; k < measures.size(); ++k){
    
    measures[k].smooth(x_v);
    
    chi2 += measures[k].KF()->chi2();
    
  }
  
  ndof = 2*measures.size()-2;
  cout << "Chi2 of the vertex: " << chi2 << " with n.d.f.: " << ndof << " Probability: " << TMath::Prob(chi2, ndof) << endl;
  return make_pair(chi2, ndof);


  
  // // Search for sub-vertices
  // vector<Measurement> ensamble1;
  // vector<Measurement> ensamble2;

  // for(uint k = 0; k < measures.size(); ++k){
    
  //   measures[k].inverseFilter(x_v,C_n);
  //   cout << "Measurement " << k << ": chi2_S = " << measures[k].KF()->chi2_S() << ", reduced: " << measures[k].KF()->chi2_S()/(ndof-2) 
  // 	 << " Prob: " << TMath::Prob(measures[k].KF()->chi2_S(), ndof-2) <<endl;

  //   if(TMath::Prob(measures[k].KF()->chi2_S(), ndof-2) < 0.05)
  //     ensamble2.push_back(measures[k]);
  //   else
  //     ensamble1.push_back(measures[k]);
  // }
  // return make_pair(ensamble1, ensamble2);
  
}

void vertex(){
  cout << "ciao" << endl;

  vector<Measurement> measures0;
  measures0.push_back(Measurement::convert(50.353, 0.00001, 0.224, 0.1, 10, 0.1));
  measures0.push_back(Measurement::convert(94.781, 0.00001, 0.307, 0.1, 10, 0.1));

  vector<Measurement> measures45;
  measures45.push_back(Measurement::convert(49.332, -1     , 0.259, 0.1, 10, 0.1));
  measures45.push_back(Measurement::convert(93.909, -1     , 0.229, 0.1, 10, 0.1)); 

  vector<Measurement> measures120;
  measures120.push_back(Measurement::convert(60.367,  1.732 , 0.229, 0.1, 10, 0.1)); 

  vector<Measurement> measures240;
  measures240.push_back(Measurement::convert(57.190, -1.732 , 0.280, 0.1, 10, 0.1)); 
  measures240.push_back(Measurement::convert(95.361, -1.732 , 0.182, 0.1, 10, 0.1)); 

  vector<vector<Measurement> > measSets;

  for(vector<Measurement>::const_iterator it0 = measures0.begin(); it0 != measures0.end(); ++it0)
    for(vector<Measurement>::const_iterator it45 = measures45.begin(); it45 != measures45.end(); ++it45)
      for(vector<Measurement>::const_iterator it120 = measures120.begin(); it120 != measures120.end(); ++it120)
	for(vector<Measurement>::const_iterator it240 = measures240.begin(); it240 != measures240.end(); ++it240){
	  vector<Measurement> a;
	  a.push_back(*it0); a.push_back(*it45); a.push_back(*it120); a.push_back(*it240);
	  measSets.push_back(a);
	}

  cout << "Number of starting ensambles " << measSets.size() << endl;

  //vector<Measurement> measures;
  //measures.push_back(Measurement(1,2,0.5,0.1,0.1,0.01,0.,0.,0.));
  //measures.push_back(Measurement(4,2,-0.2,0.1,0.1,0.01,0.,0.,0.));
  //measures.push_back(Measurement(4,2,1,0.1,0.1,0.01,0.,0.,0.));
  //measures.push_back(Measurement(1,2,1,0.1,0.1,0.01,0.,0.,0.));

  //measures.push_back(Measurement::convert(50.353, 0.00001, 0.224, 0.1, 10, 0.1));
  //measures.push_back(Measurement::convert(94.781, 0.00001, 0.307, 0.1, 10, 0.1));
  //measures.push_back(Measurement::convert(49.332, -1     , 0.259, 0.1, 10, 0.1));
  //measures.push_back(Measurement::convert(93.909, -1     , 0.229, 0.1, 10, 0.1)); 
  //measures.push_back(Measurement::convert(60.367,  1.732 , 0.229, 0.1, 10, 0.1)); 
  //measures.push_back(Measurement::convert(57.190, -1.732 , 0.280, 0.1, 10, 0.1)); 
  //measures.push_back(Measurement::convert(95.361, -1.732 , 0.182, 0.1, 10, 0.1)); 
 
  // -------------------------------------------------------------------------
  //cout<<measures.size()<<endl;

  vector<TCanvas*> canvasses;

  vector<int> goodEnsambles;

  for(uint e = 0; e < measSets.size(); ++e){
    cout << "\n\n\n\n Ensamble #" << e << endl;
    string res;
    ostringstream convert;
    convert << e;
    res = convert.str();

    canvasses.push_back(new TCanvas(TString("c_"+res)  , TString("measurements "+res)  , 200,10,600,400));
    pair<double, int> chi2_ndof = searchVertex(measSets[e],canvasses[e]);
    if(TMath::Prob(chi2_ndof.first, chi2_ndof.second) > 0.05)
      goodEnsambles.push_back(e);
  }

  for(uint v = 0; v < goodEnsambles.size(); ++v){
    cout << goodEnsambles[v] << endl;
  }
    

  //if(split.second.size() == 0) return;

  //return;
  
  //cout << "Search for vertex #1" << endl;
  // pair<vector<Measurement>, vector<Measurement> > split1 = searchVertex(split.first,3); // green

  //  cout << "Search for vertex #2" << endl;
  //pair<vector<Measurement>, vector<Measurement> > split2 = searchVertex(split.second,4); // blue
  

}
