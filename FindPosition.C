#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TString.h>
#include <TSystem.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TMatrixD.h"

using namespace std;

double* rateVsSpace(int);
void   FindPosition();
void FillHistoPlane(int, double, double);
void ClusterFinder();

const int nangles = 4;
const int nruns = 22;
const int nfnc  = 2*nangles;

double angleD[nangles] = {0.,45.,120.,240.};
double angleR[nangles];
int suff[nruns] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150};

vector<Double_t> pos;
vector<Double_t> errpos;
vector<Double_t> weight;
vector<Double_t> angcoef;


double nbinsX = 282;
double nbinsY = 242;
double xmin = -70.5;
double xmax = 70.5;
double ymin = -0.5;
double ymax = 120.5;
TH2F* histoPlane = new TH2F("histoPlane","source positions 2d image; x (mm); y (mm)",nbinsX,xmin,xmax,nbinsY,ymin,ymax);

//_______________________________________________________________________________________

void FindPosition() {
    
    for(int i = 0; i < nangles; i++) angleR[i] = TMath::Pi()*angleD[i]/180.;
    TF1 *f0pol1[nangles];
    TF1 *f1pol1[nangles];
    Double_t *ymeas;
   // Double_t m1[nangles];
    Double_t q1[nangles];
    Double_t q2[nangles];
    vector<Double_t> m;
    vector<Double_t> q;
    vector<Double_t> errq;
    
    pos.clear();
    errpos.clear();
    weight.clear();
    angcoef.clear();
    
    for(int i = 0; i < nangles; i++) {
        ymeas = rateVsSpace(i);
        f0pol1[i] = 0x0;
        f1pol1[i] = 0x0;
        if(ymeas[0]!=-999) {
            m.push_back(-TMath::Tan(angleR[i]));
            q.push_back(ymeas[0]);
            errq.push_back(ymeas[1]);
        }
        if(ymeas[2]!=-999) {
            m.push_back(-TMath::Tan(angleR[i]));
            q.push_back(ymeas[2]);
            errq.push_back(ymeas[3]);
        }
    }
    
    double errm = 0.099668;
    TCanvas *c4 = new TCanvas("c4","max counts positions",200,10,700,500);
    //    c4->Divide(2,2);
    //    c4->cd(1);
    TH1F *h = new TH1F("h","max counts positions; x (mm); y (mm)",400,-100,100);
    h->GetYaxis()->SetRangeUser(ymeas[0]-100.,ymeas[0]+100.);
    h->DrawCopy();
    
    TF1 *invf[nfnc];
    
    for(int ii = 0; ii < (int)q.size(); ii++) {
        TF1 *f = new TF1(Form("f%d",ii),"pol1",-100,100);
        invf[ii] = new TF1(Form("invf%d",ii),"pol1",-100,100);
        f->SetParameter(0,q[ii]);
        f->SetParameter(1,m[ii]);
        f->SetLineColor(ii+1);
        if(ii==0) f->Draw();
        else
   /*     if(ii==1 || ii==3 || ii==4|| ii==7)*/f->Draw("same");
        if(m[ii]==0) {
            m.at(ii)=0.001;
            errm = 0.001;
            q.at(ii) = 0.001;
            errq.at(ii) = 0.001;
        }
        else errm = 0.099668/(m[ii]*m[ii]);
        invf[ii]->SetParameter(1,1/m[ii]);
        invf[ii]->SetParameter(0,q[ii]/m[ii]);
        invf[ii]->SetParError(0,q[ii]/m[ii]*TMath::Sqrt((1/m[ii])*(1/m[ii])*(errq[ii]/q[ii])*(errq[ii]/q[ii])+q[ii]*q[ii]*(errm/m[ii])*(errm/m[ii])));
        printf("-------> %f +- %f\n",1/m[ii],q[ii]/m[ii]*TMath::Sqrt((1/m[ii])*(1/m[ii])*(errq[ii]/q[ii])*(errq[ii]/q[ii])+q[ii]*q[ii]*(errm/m[ii])*(errm/m[ii])));
        invf[ii]->SetParError(1,errm);
        invf[ii]->SetLineColor(ii+1);
        invf[ii]->SetLineStyle(kDashed);
//        if(ii==2 || ii==4)invf[ii]->Draw("same");
        //else f->Draw("same");
        f->Print();
        invf[ii]->Print();
    }
    
    ///////////////////// Kalman filter implementation
    
    
    TMatrixD Ak(2,2);
    TMatrixD Bk(2,1);
    TMatrixD Akt(2,2);
    TMatrixD Bkt(1,2);
    TMatrixD Vk(2,2);
    TMatrixD Gk(2,2);
    TMatrixD GkB(2,2);
    TMatrixD Wk(1,1);
    TMatrixD Ck(2,2);//?
    TMatrixD Cke(2,1);//?
    TMatrixD tmpCk(2,2);//??
    TMatrixD vPosEstimate(2,1);//??
    TMatrixD tmpvPosEstimate(2,1);//??
    TMatrixD qEstimate(1,1);//??
    TMatrixD Pk(2,1);//??
    TMatrixD Ptildek(2,1);//??
    TMatrixD m0(2,2);
    TMatrixD m1(1,2);
    TMatrixD m2(1,2);
    TMatrixD chi2(1,1);
    TMatrixD tmpChi2(1,1);
    
    
    tmpCk[0][0] = 10000000000;
    tmpCk[0][1] = 10000000000;
    tmpCk[1][0] = 10000000000;
    tmpCk[1][1] = 10000000000;
    
    tmpvPosEstimate[0][0] = 20;
    tmpvPosEstimate[1][0] = 20;

    Double_t xref = 10.;
    chi2[0][0]    = 0.;
    tmpChi2[0][0] = 0.;
    //Set the elements
    for(int ii = 2; ii < 4/*(int)q.size()*/; ii++) {

        if(ii==3 || ii==7 || ii==7) continue;
        Double_t yref = invf[ii]->GetX(xref);
  

        //Pk
        Pk[0][0] = yref;
        Pk[1][0] = 1/m[ii];
        
        //Ak matrix
        Ak[0][0] = 1;
        Ak[0][1] = -1/m[ii];
        Ak[1][0] = 0;
        Ak[1][1] = 0;
        Akt.Transpose(Ak);
        
        //Bk matrix
        Bk[0][0] = yref;
        Bk[1][0] = 1;
        Bkt.Transpose(Bk);

        printf("*** Ak ***\n");
        Ak.Print();
        printf("*** Bk ***\n");
        Bk.Print();
        printf("*** Bkt ***\n");
        Bkt.Print();
        
        //Vk matrix
        Double_t qVar = /*errq[ii];*/(m[ii]/q[ii]*TMath::Sqrt((1/m[ii])*(1/m[ii])*(errq[ii]/q[ii])*(errq[ii]/q[ii])+q[ii]*q[ii]*(errm/m[ii])*(errm/m[ii])));
        Double_t mVar = TMath::ATan(1/10.)/(m[ii]);//*m[ii]);
        Double_t xVar = TMath::Sqrt(qVar*qVar+mVar*mVar);
        Vk[0][0] = xVar*xVar;
        Vk[0][1] = xVar*mVar;
        Vk[1][0] = mVar*xVar;
        Vk[1][1] = mVar*mVar;

        printf("*** Vk ***\n");

        Vk.Print();
        //Gk
        Gk = Vk;
        Gk.SetTol(1.e-50);
        Double_t *determ_ptr;
        Gk.Invert(determ_ptr);
        printf("*** Gk ***\n");
        Gk.Print();

        //Wk
        Wk = Bkt*Gk*Bk;
        printf("*** Wk ***\n");
        determ_ptr=0;
        Wk.SetTol(1.e-50);
        Wk.Invert(determ_ptr);
        Wk.Print();

        //GkB
        GkB = Gk*Bk*Wk*Bkt*Gk;
        GkB = Gk - GkB;
        printf("*** Gkb ***\n");
        GkB.Print();
        
        //Ck
        Ck = Akt*GkB*Ak;
        Ck.SetTol(1.e-50);
        tmpCk.SetTol(1.e-50);
        Ck += tmpCk.Invert(determ_ptr);
        Ck.Invert(determ_ptr);
        printf("*** Ck ***\n");
        Ck.Print();

        //Cke
        Cke[0][0] = yref;
        Cke[1][0] = 1/m[ii];
        Cke = Cke - Bk*(1/m[ii]);
        
        tmpCk.SetTol(1.e-50);
        vPosEstimate = tmpCk.Invert(determ_ptr)*tmpvPosEstimate;
        m0 = Akt*GkB;
        vPosEstimate = vPosEstimate + m0*(Pk-Cke);
        vPosEstimate = Ck*vPosEstimate;
        
        qEstimate = Wk*Bkt*Gk*(Pk - Cke - Ak*vPosEstimate);
        
        Ptildek = Cke + Ak*vPosEstimate + Bk*qEstimate;
        
        chi2 = tmpChi2 + m1.Transpose(vPosEstimate-tmpvPosEstimate)*tmpCk.Invert(determ_ptr)*(vPosEstimate-tmpvPosEstimate) + m2.Transpose(Pk - Ptildek)*Gk*(Pk - Ptildek);

        printf("\n*******************************\n");
        printf("STEP %d:\n",ii);
        printf("1) posizione del vertice:\n ");
        vPosEstimate.Print();
        printf("2) direzione:\n ");
        qEstimate.Print();
        printf("3) Chi2: %f, Chi2/ndf = %f\n ",chi2[0][0],chi2[0][0]/(ii+1-2));

        printf("\n*******************************\n\n");
        tmpChi2 = chi2;
        tmpCk = Ck;
        tmpvPosEstimate = vPosEstimate;
    }
    
    
    ////////////////////////// Now let's eliminate one by one track
    
    
    
    
//    //Find intersection points
//    vector<Double_t> xint;
//    vector<Double_t> yint;
//    vector<Double_t> errxint;
//    vector<Double_t> erryint;
//    for(int ii = 0; ii < (int)q.size()-1; ii++) {
//        for(int jj = ii+1; jj < (int)q.size(); jj++) {
//            if(m[jj]!=m[ii]) {
//                xint.push_back((q[jj]-q[ii])/(m[ii]-m[jj]));
//                errxint.push_back(TMath::Sqrt(errq[jj]*errq[jj]+errq[ii]*errq[ii])/(m[ii]-m[jj]));
//                yint.push_back((m[ii]*q[jj]-m[jj]*q[ii])/(m[ii]-m[jj]));
//                erryint.push_back(TMath::Sqrt(m[ii]*m[ii]*errq[jj]*errq[jj]+m[jj]*m[jj]*errq[ii]*errq[ii])/(m[ii]-m[jj]));
//            }  else continue;
//        }
//    }
//    const int n = xint.size();
//    Double_t Xint[n], Yint[n], errXint[n], errYint[n], w[n];
//    for(int ii = 0; ii < (int)xint.size(); ii++) {
//        
//        Xint[ii] = xint[ii];
//        Yint[ii] = yint[ii];
//        errXint[ii] = errxint[ii];
//        errYint[ii] = erryint[ii];
//        w[ii] = 1.;
//    }
//    
//    TGraphErrors *gr2 = new TGraphErrors(n,Xint,Yint,errXint,errYint);
//    gr2->SetMarkerStyle(20);
//    gr2->SetMarkerSize(0.4);
//    gr2->Draw("sameP");
//    
//    q1[0] = q[1];
//    q1[1] = q[3];
//    q1[2] = q[4];
//    q1[3] = q[7];
//    m1[0] = m[1];
//    m1[1] = m[3];
//    m1[2] = m[4];
//    m1[3] = m[7];
    
//    q1[0] = q[1];
//    q1[1] = q[2];
//    q1[2] = q[5];
    
 //   ClusterFinder();
    
    //    c4->cd(2);
    //    TGraph *gr3 = new TGraph(n,Xint,w);
    //    gr3->SetMarkerStyle(20);
    //    gr2->SetMarkerSize(0.4);
    //    gr3->Draw();
    //
    //    c4->cd(3);
    //    TGraph *gr4 = new TGraph(n,Yint,w);
    //    gr4->SetMarkerStyle(20);
    //    gr2->SetMarkerSize(0.4);
    //    gr4->Draw();
    
    
//        Double_t aa = 0., bb = 0., cc = 0., dd = 0., ee = 0., ff = 0., gg = 0.;
//        for(int i = 0; i < nangles; i++) {
//            aa += q1[i]/(1+m1[i]*m1[i]);
//            bb += m1[i]/(1+m1[i]*m1[i]);
//            cc += q1[i]*m1[i]/(1+m1[i]*m1[i]);
//            dd += 1/(1+m1[i]*m1[i]);
//            ee += m1[i]*m1[i]/(1+m1[i]*m1[i]);
//        }
//        Double_t xs1[1];
//        xs1[0] = (aa*bb - cc*dd)/(ee*dd - bb*bb);
//    
//        Double_t ys1[1];
//        ys1[0] = (xs1[0]*bb+aa)/dd;
//        printf("xp = %f, yp = %f\n",xs1[0],ys1[0]);
//        TGraph *gr = new TGraph(1,xs1,ys1);
//        gr->SetMarkerStyle(20);
//        gr->Draw("sameP");
    //
    //    Double_t aa2 = 0., bb2 = 0., cc2 = 0., dd2 = 0., ee2 = 0., ff2 = 0., gg2 = 0.;
    //    for(int i = 0; i < nangles; i++) {
    //        aa2 += q2[i]/(1+m[i]*m[i]);
    //        bb2 += m[i]/(1+m[i]*m[i]);
    //        cc2 += q2[i]*m[i]/(1+m[i]*m[i]);
    //        dd2 += 1/(1+m[i]*m[i]);
    //        ee2 += m[i]*m[i]/(1+m[i]*m[i]);
    //    }
    //    Double_t xs2[1];
    //    xs2[0] = (aa2*bb2 - cc2*dd2)/(ee2*dd2 -bb2*bb2);
    //    Double_t ys2[1];
    //    ys2[0] = (xs2[0]*bb2+aa2)/dd2;
    //    printf("xp2 = %f, yp2 = %f\n",xs2[0],ys2[0]);
    //    TGraph *gr2 = new TGraph(1,xs2,ys2);
    //    gr2->SetMarkerStyle(20);
    //    gr2->Draw("sameP");
    //
    //    printf("-------> %f\n",f0pol1[1]->Eval());
    
    
    
    
        TCanvas *c3 = new TCanvas("c3","source position 2d",200,10,700,700);
        histoPlane->Draw("colz");
}

//_______________________________________________________________________________________
Double_t* rateVsSpace(int iangle) {
    
    TString dir = "../20160616/";
    TString filename;
    
    Double_t X[nruns];
    Double_t rate[nruns];
    Double_t errX[nruns];
    Double_t errRate[nruns];
    
    TFile *f;
    TTree *t1;
    
    TH1F *hRawEnergyCh0;
    TH1F *hRawEnergyCh1;
    TH1F *hDeltat;
    
    hRawEnergyCh0 = 0x0;
    hRawEnergyCh1 = 0x0;
    hDeltat = 0x0;
    int maxEntries = 1;
    
    for(int irun = 0; irun < nruns; irun++) {
        
        double tmpmintime0;
        double tmpmintime1;
        double tmpmaxtime0;
        double tmpmaxtime1;
        
        if(hRawEnergyCh0) delete hRawEnergyCh0;
        if(hRawEnergyCh1) delete hRawEnergyCh1;
        if(hDeltat) delete hDeltat;
        hRawEnergyCh0 = new TH1F("hRawEnergyCh0","; keV; entries",600,0,3000);
        hRawEnergyCh1 = new TH1F("hRawEnergyCh1","; keV; entries",600,0,3000);
        hDeltat = new TH1F("hDeltat","#Deltat (ch1 - ch0); #Deltat",20,-95,105);
        
        f = TFile::Open(Form("%sAnalysisResults_Calibrated_run_20160616_2Na22_%d_%.0f.root",dir.Data(),suff[irun],angleD[iangle]));
        if(f) {
            t1 = (TTree*)f->Get("t");
            if(t1) {
                Double_t en0, en1, time0, time1;
                t1->SetBranchAddress("time0",&time0);
                t1->SetBranchAddress("en0",&en0);
                t1->SetBranchAddress("time1",&time1);
                t1->SetBranchAddress("en1",&en1);
                
                Long64_t nentries = t1->GetEntries();
                Double_t dt;
                
                t1->GetEntry(0);
                tmpmintime0 = time0;
                tmpmintime1 = time1;
                tmpmaxtime0 = time0;
                tmpmaxtime1 = time1;
                
                for (Long64_t iev = 0; iev < nentries; iev++) {
                    
                    t1->GetEntry(iev);
                    if(time0<tmpmintime0) tmpmintime0 = time0;
                    if(time1<tmpmintime1) tmpmintime1 = time1;
                    if(time0>tmpmaxtime0) tmpmaxtime0 = time0;
                    if(time1>tmpmaxtime1) tmpmaxtime1 = time1;
                    
                    dt = 10*(time1 - time0);
                    if((en0 > 460 && en0 < 560) && (en1 > 460 && en1 < 550) && TMath::Abs(dt)<20) {
                        hDeltat->Fill(dt);
                        hRawEnergyCh0->Fill(en0);
                        hRawEnergyCh1->Fill(en1);
                    }
                }
            }
        } else {
            printf("file %s not there\n",Form("%sAnalysisResults_Calibrated_run_20160616_2Na22_%d_%.0f.root",dir.Data(),suff[irun],angleD[iangle]));
            rate[irun] = 0.;
            errRate[irun] = 0.;
            continue;
        }
        double time = TMath::Power(10,-8)*(tmpmaxtime0+tmpmaxtime0-tmpmintime0-tmpmintime1)/2.;
        double sig = (double)hRawEnergyCh0->GetEntries();
        rate[irun] = sig/time;
        errX[irun] = 0.1;
        errRate[irun] = TMath::Sqrt(sig)/time;
        X[irun] = suff[irun]/10.;
        
        if(rate[irun]>0) {
            double m = -TMath::Tan(angleR[iangle]);
            pos.push_back(X[irun]);
            errpos.push_back(errX[irun]);
            angcoef.push_back(m);
            weight.push_back(rate[irun]);
        }
        FillHistoPlane(iangle,X[irun],rate[irun]);
        
        //                        TCanvas *c = new TCanvas(Form("c%d",irun),"",800,500);
        //                        c->Divide(2);
        //                        c->cd(1);
        //                        hDeltat->DrawCopy();
        //                        c->cd(2);
        //                        hRawEnergyCh0->SetLineColor(kBlue);
        //                        hRawEnergyCh1->SetLineColor(kRed);
        //                        hRawEnergyCh0->DrawCopy();
        //                        hRawEnergyCh1->DrawCopy("same");
        
    }
    
    
    vector<double> localMaxima;
    int counter=0;
    for (int i = 1; i < nruns-1; i++)
        if (rate[i] > rate[i-1] && rate[i] > rate[i+1] && rate[i]>0.5) {
            localMaxima.push_back(X[i]);
            printf("-------------------------> %f\n",localMaxima[counter]);
            counter++;
        }
    int ind1largest = -1;
    int ind2largest = -1;
    if(localMaxima.size()>2) {
        printf("be careful, more than 2 local maxima have been found!\n");
        int largest = -1;
        int second_largest = -1;
        
        for(int j = 0; j < (int)localMaxima.size(); j++) { // find the largest
            if (localMaxima[j] >= largest) {
                largest = localMaxima[j];
                ind1largest = j;
            }
        }
        for (int k = 0; k < (int)localMaxima.size(); k++) { // find the second largest
            if (k != ind1largest) {// skip over the largest one
                if (localMaxima[k] >= second_largest) {
                    second_largest = localMaxima[k];
                    ind2largest = k;
                }
            }
        }
    } else {
        ind1largest = 0;
        ind2largest = 1;
    }
    
    
    printf("*****************************************************\n");
    printf("*****************************************************\n");
    printf("%f, %f\n",localMaxima[ind1largest],localMaxima[ind2largest]);
    printf("*****************************************************\n");
    printf("*****************************************************\n");
    TGraphErrors *gr = new TGraphErrors(nruns,X,rate,errX,errRate);
    TCanvas *c1 = new TCanvas(Form("c111%d",iangle),"Rate vs Source Position",200,10,700,500);
    gr->SetTitle(Form("Rate vs Source Position, angle %.0f#circ",angleD[iangle]));
    gr->GetXaxis()->SetTitle("position (mm)");
    gr->GetYaxis()->SetTitle("rate (s^{-1})");
    gr->Draw("AP");
    
    Double_t mean[4];
    for(int m = 0; m < 4; m++) mean[m] = -999.;
    TF1 *f1gaus = new TF1(Form("f1gaus%d",iangle),"gaus",localMaxima[ind1largest]-15,localMaxima[ind1largest]+15);
    f1gaus->SetParameter(1,localMaxima[ind1largest]);
    f1gaus->SetParameter(2,5.);
    //gStyle->SetOptFit(kTRUE);
    gr->Fit(f1gaus,"R");
    mean[0] = f1gaus->GetParameter(1);
    mean[1] = f1gaus->GetParError(1);
    f1gaus->SetLineColor(kBlue);
    f1gaus->DrawCopy("same");
    
    if(localMaxima.size()>1 and localMaxima[ind2largest]>10) {
        TF1 *f2gaus = new TF1(Form("f2gaus%d",iangle),"gaus",localMaxima[ind2largest]-15,localMaxima[ind2largest]+15);
        f2gaus->SetParameter(1,localMaxima[ind2largest]);
        f2gaus->SetParameter(2,5.);
        gr->Fit(f2gaus,"R");
        mean[2] = f2gaus->GetParameter(1);
        mean[3] = f2gaus->GetParError(1);
        f2gaus->SetLineColor(kBlue);
        f2gaus->DrawCopy("same");
    }
    c1->SaveAs(Form("RateVsSpace_%.0f.eps",angleD[iangle]));
    //    f1->SetParameter(1,localMaxima[ind2largest]);
    //    f1->SetParameter(2,5.);
    //    gr->Fit(fgaus);
    //    double mean = fgaus->GetParameter(1);
    
    printf("VALUES from fit: %.3f +- %.3f, %.3f +- %.3f\n",mean[0],mean[1],mean[2],mean[3]);
    return mean;
}

//_______________________________________________________________________________________

void FillHistoPlane(int iangle, double q, double rate) {
    
    double wbinX = histoPlane->GetXaxis()->GetBinWidth(1);
    double wbinY = histoPlane->GetYaxis()->GetBinWidth(1);
    
    double m =  -TMath::Tan(angleR[iangle]);
    TF1 *fpol1 = new TF1(Form("fpol1%d",iangle),"pol1",-100,100);
    fpol1->SetParameter(0,q);
    fpol1->SetParameter(1,m);
    
    for(int ibinX = 1; ibinX <=nbinsX; ibinX++) {
        double x = xmin +(ibinX-1)*wbinX;
        histoPlane->Fill(x,fpol1->Eval(x),rate);
        for(int ib = 1; ib <=10; ib++) {
            histoPlane->Fill(x,fpol1->Eval(x)+ib*wbinY,rate);
            histoPlane->Fill(x,fpol1->Eval(x)-ib*wbinY,rate);
        }
    }
}

//_______________________________________________________________________________________

void ClusterFinder() {
    
    
    TF1 *fline = 0x0;
    TCanvas *c5 = new TCanvas("c5","source position 2d",200,10,700,700);
    TH2F* histo2d = new TH2F("histo2d","histo 2d; x (mm); y (mm)",nbinsX,xmin,xmax,nbinsY,ymin,ymax);
    histo2d->DrawCopy();
    
    for(int i = 0; i < (int)pos.size(); i++) {
        if(fline) delete fline;
        fline = new TF1(Form("fline%d",i),"pol1",-100,100);
        fline->SetParameter(0,pos[i]);
        fline->SetParameter(1,angcoef[i]);
        
        fline->DrawCopy("same");
    }
    
    
    
    
    
}













