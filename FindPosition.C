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

double rateVsSpace(double angle) {
    
    TString dir = "20160613_2cm/";
    TString filename;
    
    const int nruns = 6;
    int suff[nruns] = {650,675,700,725,750,800};
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
        
        hRawEnergyCh0 = new TH1F("hRawEnergyCh0","; keV; entries",600,0,3000);
        hRawEnergyCh1 = new TH1F("hRawEnergyCh1","; keV; entries",600,0,3000);
        hDeltat = new TH1F("hDeltat","#Deltat (ch1 - ch0); #Deltat",20,-95,105);
        
        f = new TFile(Form("%sAnalysisResults_Calibrated_run_20160613_%d_%.0f_2cm.root",dir.Data(),suff[irun],angle));
        
        double tmpmintime0;
        double tmpmintime1;
        double tmpmaxtime0;
        double tmpmaxtime1;
        
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
                    if((en0 > 460 && en0 < 580) && (en1 > 460 && en1 < 550) && TMath::Abs(dt)<20) {
                        hDeltat->Fill(dt);
                        hRawEnergyCh0->Fill(en0);
                        hRawEnergyCh1->Fill(en1);
                    }
                }
            }
        }
        
        double time = TMath::Power(10,-8)*(tmpmaxtime0+tmpmaxtime0-tmpmintime0-tmpmintime1)/2.;
        double sig = (double)hRawEnergyCh0->GetEntries();
        rate[irun] = sig/time;
        errX[irun] = 0.1;
        errRate[irun] = TMath::Sqrt(sig)/time;
        X[irun] = suff[irun]/10.;
//                TCanvas *c = new TCanvas(Form("c%d",irun),"",800,500);
//                c->Divide(2);
//                c->cd(1);
//                hDeltat->DrawCopy();
//                c->cd(2);
//                hRawEnergyCh0->SetLineColor(kBlue);
//                hRawEnergyCh1->SetLineColor(kRed);
//                hRawEnergyCh0->DrawCopy();
//                hRawEnergyCh1->DrawCopy("same");
        
    }
    
    TGraphErrors *gr = new TGraphErrors(nruns,X,rate,errX,errRate);
    TCanvas *c1 = new TCanvas("c1","Rate vs Source Position",200,10,700,500);
    gr->SetTitle("Rate vs Source Position");
    gr->GetXaxis()->SetTitle("position (mm)");
    gr->GetYaxis()->SetTitle("rate (s^{-1})");
    gr->Draw("AP");
    
    TF1 *fgaus = new TF1("fgaus","gaus");
    gr->Fit(fgaus);
    double mean = fgaus->GetParameter(1);
    
    return mean;
}

//_______________________________________________________________________________________

void FindPosition() {
    
    double angle[3] = {0.,120.,240.};
    
    for(int i = 0; i < 3; i++) angle[i] = TMath::Pi()*angle[i]/180.;
    TF1 *fpol1[3];
    double ymeas;
    for(int i = 0; i < 3; i++) {
        
        ymeas = rateVsSpace(0.);
        double m =  TMath::Tan(angle[i]);
//        if((angle[i]>90 && angle[i]<180)) TMath::Tan(angle[i]);
        
        fpol1[i] = new TF1(Form("fpol1%d",i),"pol1",-5,5);
        fpol1[i]->SetParameter(0,ymeas);
        fpol1[i]->SetParameter(1,m);
        
    }
    
    TCanvas *c2 = new TCanvas("c2","prova",200,10,700,500);
    TH1F *h = new TH1F("h","prova",60,-10,10);
    h->GetYaxis()->SetRangeUser(ymeas-20.,ymeas+20.);
    h->DrawCopy();
    for(int i = 0; i < 3; i++) {
        fpol1[i]->SetLineColor(i+1);
        fpol1[i]->Draw("same");
        fpol1[i]->Print();
    }
}











