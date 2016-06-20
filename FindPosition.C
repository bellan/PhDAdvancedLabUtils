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
#include "TF1.h"
#include "TTree.h"
#include "TGraphErrors.h"

using namespace std;

double* rateVsSpace(int);
void   FindPosition();

const int nangles = 3;
const int nruns = 22;

double angleD[nangles] = {0.,120.,240.};
double angleR[nangles];
int suff[nruns] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150};

//_______________________________________________________________________________________

void FindPosition() {
    
    for(int i = 0; i < nangles; i++) angleR[i] = TMath::Pi()*angleD[i]/180.;
    TF1 *f0pol1[nangles];
    TF1 *f1pol1[nangles];
    Double_t *ymeas;
    Double_t m[nangles];
    Double_t q1[nangles];
    Double_t q2[nangles];
    for(int i = 0; i < nangles; i++) {
        ymeas = rateVsSpace(i);
        f0pol1[i] = 0x0;
        f1pol1[i] = 0x0;
        printf("RISULTATI ----------------> %f, %f\n",ymeas[0],ymeas[1]);
        m[i] =  -TMath::Tan(angleR[i]);
        if(ymeas[0]!=-999) {
            f0pol1[i] = new TF1(Form("f0pol1%d",i),"pol1",-100,100);
            f0pol1[i]->SetParameter(0,ymeas[0]);
            f0pol1[i]->SetParameter(1,m[i]);
            q1[i] = ymeas[0];
        }
        if(ymeas[1]!=-999) {
            f1pol1[i] = new TF1(Form("f1pol1%d",i),"pol1",-100,100);
            f1pol1[i]->SetParameter(0,ymeas[1]);
            f1pol1[i]->SetParameter(1,m[i]);
            q2[i] = ymeas[1];
        }
    }
    TCanvas *c2 = new TCanvas("c2","prova",200,10,700,500);
    TH1F *h = new TH1F("h","prova",400,-100,100);
    h->GetYaxis()->SetRangeUser(ymeas[0]-100.,ymeas[0]+100.);
    h->DrawCopy();
    for(int i = 0; i < nangles; i++) {
        if(f0pol1[i]) {
            f0pol1[i]->SetLineColor(i+1);
            f0pol1[i]->Draw("same");
            f0pol1[i]->Print();
        }
        if(f1pol1[i]) {
            f1pol1[i]->SetLineColor(i+6);
            f1pol1[i]->Draw("same");
            f1pol1[i]->Print();
        }
    }
    
    Double_t aa = 0., bb = 0., cc = 0., dd = 0., ee = 0., ff = 0., gg = 0.;
    for(int i = 0; i < nangles; i++) {
        aa += q1[i]/(1+m[i]*m[i]);
        bb += m[i]/(1+m[i]*m[i]);
        cc += q1[i]*m[i]/(1+m[i]*m[i]);
        dd += 1/(1+m[i]*m[i]);
        ee += m[i]*m[i]/(1+m[i]*m[i]);
    }
    Double_t xs1[1];
    xs1[0] = (aa*bb - cc*dd)/(ee*dd - bb*bb);

    Double_t ys1[1];
    ys1[0] = (xs1[0]*bb+aa)/dd;
    printf("xp = %f, yp = %f\n",xs1[0],ys1[0]);
    TGraph *gr = new TGraph(1,xs1,ys1);
    gr->SetMarkerStyle(20);
    gr->Draw("sameP");

    Double_t aa2 = 0., bb2 = 0., cc2 = 0., dd2 = 0., ee2 = 0., ff2 = 0., gg2 = 0.;
    for(int i = 0; i < nangles; i++) {
        aa2 += q2[i]/(1+m[i]*m[i]);
        bb2 += m[i]/(1+m[i]*m[i]);
        cc2 += q2[i]*m[i]/(1+m[i]*m[i]);
        dd2 += 1/(1+m[i]*m[i]);
        ee2 += m[i]*m[i]/(1+m[i]*m[i]);
    }
    Double_t xs2[1];
    xs2[0] = (aa2*bb2 - cc2*dd2)/(ee2*dd2 -bb2*bb2);
    Double_t ys2[1];
    ys2[0] = (xs2[0]*bb2+aa2)/dd2;
    printf("xp2 = %f, yp2 = %f\n",xs2[0],ys2[0]);
    TGraph *gr2 = new TGraph(1,xs2,ys2);
    gr2->SetMarkerStyle(20);
    gr2->Draw("sameP");
    
//    printf("-------> %f\n",f0pol1[1]->Eval());

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
    
    int maxEntries = 1;
    
    for(int irun = 0; irun < nruns; irun++) {
        
        double tmpmintime0;
        double tmpmintime1;
        double tmpmaxtime0;
        double tmpmaxtime1;
        
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
            continue;
        }
        double time = TMath::Power(10,-8)*(tmpmaxtime0+tmpmaxtime0-tmpmintime0-tmpmintime1)/2.;
        double sig = (double)hRawEnergyCh0->GetEntries();
        rate[irun] = sig/time;
        errX[irun] = 0.1;
        errRate[irun] = TMath::Sqrt(sig)/time;
        X[irun] = suff[irun]/10.;
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
        if (rate[i] > rate[i-1] && rate[i] > rate[i+1] && rate[i]>0.4) {
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
    TCanvas *c1 = new TCanvas(Form("c1%d",iangle),"Rate vs Source Position",200,10,700,500);
    gr->SetTitle("Rate vs Source Position");
    gr->GetXaxis()->SetTitle("position (mm)");
    gr->GetYaxis()->SetTitle("rate (s^{-1})");
    gr->Draw("AP");
    
    Double_t mean[2];
    mean[0]=-999.;
    mean[1]=-999.;
    TF1 *f1gaus = new TF1("f1gaus","gaus",localMaxima[ind1largest]-15,localMaxima[ind1largest]+15);
    f1gaus->SetParameter(1,localMaxima[ind1largest]);
    f1gaus->SetParameter(2,5.);
    gr->Fit(f1gaus,"R");
    mean[0] = f1gaus->GetParameter(1);
    
    if(localMaxima.size()>1) {
        TF1 *f2gaus = new TF1("f2gaus","gaus",localMaxima[ind2largest]-15,localMaxima[ind2largest]+15);
        f2gaus->SetParameter(1,localMaxima[ind2largest]);
        f2gaus->SetParameter(2,5.);
        gr->Fit(f2gaus,"R");
        mean[1] = f2gaus->GetParameter(1);
    }
    //    f1->SetParameter(1,localMaxima[ind2largest]);
    //    f1->SetParameter(2,5.);
    //    gr->Fit(fgaus);
    //    double mean = fgaus->GetParameter(1);
    
    return mean;
}











