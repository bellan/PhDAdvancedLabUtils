#include <iostream>
#include <sstream>
#include <fstream>
#include "TString.h"
#include "TLegend.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"

using namespace std;

//______________________________________________________________

void analyzeSpectrum(/*TString dirname ="../Na22_20160629_tripla/", TString filename = "run_20160630_Na22_teflon_tripla", Int_t nChannels = 3*/ TString dirname ="../20160629_Ba/", TString filename = "run_20160629_Ba", Int_t nChannels = 2) {
    
    TFile *f = new TFile(Form("%sAnalysisResults_Calibrated_%s.root",dirname.Data(),filename.Data()));
    if(f) {
        TTree *t1 = (TTree*)f->Get("t");
        if(t1) {
            Double_t en0, en1, en2, time0, time1, time2;
            t1->SetBranchAddress("time0",&time0);
            t1->SetBranchAddress("en0",&en0);
            t1->SetBranchAddress("time1",&time1);
            t1->SetBranchAddress("en1",&en1);
            if(nChannels>2) {
                t1->SetBranchAddress("time2",&time2);
                t1->SetBranchAddress("en2",&en2);
            }
            
            //Ba133
            TH1F *hRawEnergyCh0 = new TH1F("hRawEnergyCh0","; keV; entries",600,0,3000);
            TH1F *hRawEnergyCh1 = new TH1F("hRawEnergyCh1","; keV; entries",600,0,3000);
            TH1F *hDeltat = new TH1F("hDeltat","#Deltat (ch1 - ch0); #Deltat",200,-995,1005);
            hRawEnergyCh0->SetLineColor(kRed);
            hRawEnergyCh1->SetLineColor(kBlue);
            hRawEnergyCh0->SetStats(0);
            hRawEnergyCh1->SetStats(0);
            
            printf("Analyze events\n");
            Long64_t nentries = t1->GetEntries();
            Double_t dt;
            
            for (Long64_t iev = 0; iev < nentries; iev++) {
                
                t1->GetEntry(iev);
                dt = 10*(time1 - time0);
//                hDeltat->Fill(dt);
                
                if((en0 > 320 && en0 < 380) && (en1 > 72 && en1 < 92)) {
                    hDeltat->Fill(dt);
                    hRawEnergyCh0->Fill(en0);
                    hRawEnergyCh1->Fill(en1);
                }
                
            }
            
            hDeltat->SetTitle("#Deltat");
            
            TLegend *legend = new TLegend(0.6,0.75,0.9,0.9);
            legend->AddEntry(hRawEnergyCh0,"Canale 0");
            legend->AddEntry(hRawEnergyCh1,"Canale 1");
            legend->SetFillColor(0);
            TCanvas *c = new TCanvas("c","",1200,500);
            c->Divide(2);
            c->cd(1);
            gPad->SetLogy();
            hRawEnergyCh1->Draw();
            hRawEnergyCh0->Draw("same");
            legend->Draw("same");
            c->cd(2);
            gPad->SetLogy();
            hDeltat->Draw();
            
            
            
            //Na22
            //            TH1F *hRawEnergyCh0 = new TH1F("hRawEnergyCh0","; keV; entries",600,0,3000);
            //            TH1F *hRawEnergyCh1 = new TH1F("hRawEnergyCh1","; keV; entries",600,0,3000);
            //            TH1F *hRawEnergyCh2 = new TH1F("hRawEnergyCh2","; keV; entries",600,0,3000);
            //            TH1F *hDeltat01 = new TH1F("hDeltat01","#Deltat (ch1 - ch0); #Deltat",200,-995,1005);//20,-95,105);
            //            TH1F *hDeltat02 = new TH1F("hDeltat02","#Deltat (ch2 - ch0); #Deltat",200,-995,1005);//20,-95,105);
            //            TH1F *hDeltat3g = new TH1F("hDeltat3g","#Deltat (3#gamma); #Deltat",20,-95,105);
            //            TH1F *hDeltat2g = new TH1F("hDeltat2g","#Deltat (2#gamma); #Deltat",20,-95,105);
            //            hRawEnergyCh0->SetLineColor(kRed);
            //            hRawEnergyCh1->SetLineColor(kBlue);
            //            hRawEnergyCh2->SetLineColor(kGreen+1);
            //            hRawEnergyCh0->SetStats(0);
            //            hRawEnergyCh1->SetStats(0);
            //            hRawEnergyCh2->SetStats(0);
            //
            //            printf("Analyze events\n");
            //            Long64_t nentries = t1->GetEntries();
            //            Double_t dt01,dt02,dt21;
            //            Double_t dt3, dt2;
            //
            //            for (Long64_t iev = 0; iev < nentries; iev++) {
            //
            //                t1->GetEntry(iev);
            //                dt01 = 10*(time1 - time0);
            //                dt02 = 10*(time2 - time0);
            //                dt21 = 10*(time2 - time1);
            ////                hDeltat01->Fill(dt01);
            ////                hDeltat02->Fill(dt02);
            //                hRawEnergyCh0->Fill(en0);
            //                hRawEnergyCh1->Fill(en1);
            //                hRawEnergyCh2->Fill(en2);
            //
            //                //3 gamma from ortho-positronium <511 keV, triggered by one gamma 1.27 + 2 gammas with E < 511 keV
            //                if(en0 < 640 && en1 < 560  && (en2 > 1180 && en2 < 1350)) {
            //                    dt3 = 10*(time1-time2);
            //                    hDeltat3g->Fill(dt3);
            //                }
            //                else if(en0 < 640 && en2 < 560  && (en1 > 1180 && en1 < 1350)) {
            //                    dt3 = 10*(time2-time1);
            //                    hDeltat3g->Fill(dt3);
            //                }
            //                else if(en1 < 560 && en2 < 560  && (en0 > 1160 && en0 < 1350)) {
            //                    dt3 = 10*(time1-time0);
            //                    hDeltat3g->Fill(dt3);
            //                }
            //
            //                //2 gamma 511 back-to-back from para-positronium: triggered by 1.27 in coincidence with one around 511
            //                if((en0 > 450 && en0 < 640) && (en2 > 1180 && en2 < 1350)) {
            //                    dt2 = 10*(time0-time2);
            //                    hDeltat2g->Fill(dt2);
            //                }
            //                else if((en1 > 450 && en1 < 560)  && (en2 > 1180 && en2 < 1350)) {
            //                    dt2 = 10*(time1-time2);
            //                    hDeltat2g->Fill(dt2);
            //                }
            //                else if((en0 > 450 && en0 < 640)  && (en1 > 1180 && en1 < 1350)) {
            //                    dt2 = 10*(time0-time1);
            //                    hDeltat2g->Fill(dt2);
            //                }
            //                else if((en2 > 450 && en2 < 560)  && (en1 > 1180 && en1 < 1350)) {
            //                    dt2 = 10*(time2-time1);
            //                    hDeltat2g->Fill(dt2);
            //                }
            //                else if((en1 > 450 && en1 < 560)  && (en0 > 1160 && en0 < 1350)) {
            //                    dt2 = 10*(time1-time0);
            //                    hDeltat2g->Fill(dt2);
            //                }
            //                else if((en2 > 450 && en2 < 560)  && (en0 > 1160 && en0 < 1350)) {
            //                    dt2 = 10*(time2-time0);
            //                    hDeltat2g->Fill(dt2);
            //                }
            //
            //
            //        }
            
            
            //            hDeltat2g->SetTitle("#Deltat");
            //            hDeltat3g->SetTitle("#Deltat");
            //
            //            TLegend *legend = new TLegend(0.6,0.75,0.9,0.9);
            //            legend->AddEntry(hRawEnergyCh0,"Canale 0");
            //            legend->AddEntry(hRawEnergyCh1,"Canale 1");
            //            legend->AddEntry(hRawEnergyCh2,"Canale 2");
            //            legend->SetFillColor(0);
            //            TCanvas *c = new TCanvas("c","",1200,500);
            //            c->Divide(2);
            //            c->cd(1);
            //            gPad->SetLogy();
            //            hRawEnergyCh1->Draw();
            //            hRawEnergyCh0->Draw("same");
            //            hRawEnergyCh2->Draw("same");
            //            legend->Draw("same");
            //            c->cd(2);
            //            gPad->SetLogy();
            //            //            hDeltat01->SetLineColor(kRed);
            //            //            hDeltat02->SetLineColor(kGreen+1);
            //            hDeltat3g->SetLineColor(kRed);
            //            hDeltat2g->SetLineColor(kBlue);
            //            hDeltat3g->Scale(1/hDeltat3g->Integral());
            //            hDeltat2g->Scale(1/hDeltat2g->Integral());
            //            hDeltat2g->Draw();
            //            hDeltat3g->Draw("sames");
            //            TLegend *legend2 = new TLegend(0.6,0.75,0.9,0.9);
            //            legend2->AddEntry(hDeltat3g,"3#gamma");
            //            legend2->AddEntry(hDeltat2g,"2#gamma");
            //            legend2->SetFillColor(0);
            //            legend2->Draw("same");
        } else printf("no TTree found!\n");
    } else printf("no data file found!\n");
}










