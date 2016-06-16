#include <iostream>
#include <sstream>
#include <fstream>
#include "TString.h"
#include "TLegend.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TTree.h"

using namespace std;

vector<Double_t> TimeCalo0,TimeCalo1,EnCalo0,EnCalo1;
TString filenameData = "";

bool ReadEvents() {
    
    TimeCalo0.clear();
    TimeCalo1.clear();
    EnCalo0.clear();
    EnCalo1.clear();
    
    ifstream inFileData(filenameData.Data());
    if(inFileData) {
        string line;
        while(getline(inFileData, line)) {
            if (line[0] == '#') continue;
            std::istringstream iss(line);
            Double_t t0,e0,t1,e1;
            iss>>t0>>e0>>t1>>e1;
            if(e0 == 16383 || e1 == 16383) continue;
            TimeCalo0.push_back(t0);
            EnCalo0.push_back(e0);
            TimeCalo1.push_back(t1);
            EnCalo1.push_back(e1);
        }
        printf("Unpacked data. Events = %lu\n",TimeCalo0.size());
    } else {
        printf("No file data found!\n");
        return false;
    }
    return true;
    
}

//______________________________________________________________

void doCalib(TString filename = "run1_20160531_Na22_Teflon.txt") {
    
    int numbit = 14;
    int maxRawEnergy = pow(2,numbit)-1;
    
    filenameData = Form("processed_%s",filename.Data());
    bool raws = ReadEvents();
    Double_t ePeak[3] = {511,1274,1785};
    Double_t lowlimsCh0[3], lowlimsCh1[3], uplimsCh0[3], uplimsCh1[3];
    
    if(filename == "run3_260516.txt" || filename == "run4Rot_260516.txt") {
        lowlimsCh0[0] = 2600;
        lowlimsCh0[1] = 6500;
        lowlimsCh0[2] = 9300;
        uplimsCh0[0] = 3000;
        uplimsCh0[1] = 7100;
        uplimsCh0[2] = 9700;
        
        lowlimsCh1[0] = 2600;
        lowlimsCh1[1] = 6400;
        lowlimsCh1[2] = 9100;
        uplimsCh1[0] = 2850;
        uplimsCh1[1] = 6800;
        uplimsCh1[2] = 10200;
    }
    if(filename == "run0_20160601_Na22_Al.txt") {
        lowlimsCh0[0] = 2800;
        lowlimsCh0[1] = 6800;
        lowlimsCh0[2] = 9300;
        uplimsCh0[0] = 3100;
        uplimsCh0[1] = 7300;
        uplimsCh0[2] = 9900;
        
        lowlimsCh1[0] = 2800;
        lowlimsCh1[1] = 6900;
        lowlimsCh1[2] = 9100;
        uplimsCh1[0] = 3200;
        uplimsCh1[1] = 7500;
        uplimsCh1[2] = 10400;
    }
    if(filename == "run1_20160531_Na22_Teflon.txt") {
        lowlimsCh0[0] = 2800;
        lowlimsCh0[1] = 7100;
        lowlimsCh0[2] = 9300;
        uplimsCh0[0] = 3300;
        uplimsCh0[1] = 7600;
        uplimsCh0[2] = 9900;
        
        lowlimsCh1[0] = 2800;
        lowlimsCh1[1] = 6900;
        lowlimsCh1[2] = 9100;
        uplimsCh1[0] = 3050;
        uplimsCh1[1] = 7300;
        uplimsCh1[2] = 10400;
    }
    if(filename == "run0_20160531_Na22_Teflon.txt") {
        lowlimsCh0[0] = 2800;
        lowlimsCh0[1] = 7100;
        lowlimsCh0[2] = 9300;
        uplimsCh0[0] = 3200;
        uplimsCh0[1] = 7700;
        uplimsCh0[2] = 9900;
        
        lowlimsCh1[0] = 2750;
        lowlimsCh1[1] = 6850;
        lowlimsCh1[2] = 9100;
        uplimsCh1[0] = 3100;
        uplimsCh1[1] = 7300;
        uplimsCh1[2] = 10400;
    }
    if(filename == "run_20160613_650_0.txt") {
        lowlimsCh0[0] = 2900;
        lowlimsCh0[1] = 7100;
        lowlimsCh0[2] = 9300;
        uplimsCh0[0] = 3300;
        uplimsCh0[1] = 7700;
        uplimsCh0[2] = 9900;
        
        lowlimsCh1[0] = 2750;
        lowlimsCh1[1] = 6850;
        lowlimsCh1[2] = 9100;
        uplimsCh1[0] = 3000;
        uplimsCh1[1] = 7300;
        uplimsCh1[2] = 10400;
    }
    vector<Double_t> gPeaks0, gPeaks1;
    
    TH1F *hRawEnergyCh0, *hRawEnergyCh1;
    if(raws) {
        
        hRawEnergyCh0 = new TH1F("hRawEnergyCh0","; channels; entries",512,0,maxRawEnergy);
        hRawEnergyCh1 = new TH1F("hRawEnergyCh1","; channels; entries",512,0,maxRawEnergy);
        hRawEnergyCh0->SetLineColor(kRed);
        hRawEnergyCh1->SetLineColor(kBlue);
        hRawEnergyCh0->SetStats(0);
        hRawEnergyCh1->SetStats(0);
        
        for(int i = 0; i < (int)TimeCalo1.size(); i++) {
            hRawEnergyCh0->Fill(EnCalo0[i]);
            hRawEnergyCh1->Fill(EnCalo1[i]);
        }
        
        //    }
        //        new TCanvas;
        //        hRawEnergyCh0->Draw();
        //        hRawEnergyCh1->Draw("same");
        
        TF1 *fgaus = new TF1("fgaus","gaus", 0,16000);
        //channel 1
        hRawEnergyCh1->Fit("fgaus","","",lowlimsCh1[0],uplimsCh1[0]);
        Double_t m = fgaus->GetParameter(1);
        gPeaks1.push_back(m);
        hRawEnergyCh1->Fit("fgaus","","",lowlimsCh1[1],uplimsCh1[1]);
        m = fgaus->GetParameter(1);
        gPeaks1.push_back(m);
        hRawEnergyCh1->Fit("fgaus","","",lowlimsCh1[2],uplimsCh1[2]);
        m = fgaus->GetParameter(1);
        gPeaks1.push_back(m);
        //channel 0
        //        new TCanvas;
        //        hRawEnergyCh0->SetLineColor(kBlue);
        hRawEnergyCh0->Fit("fgaus","","",lowlimsCh0[0],uplimsCh0[0]);
        m = fgaus->GetParameter(1);
        gPeaks0.push_back(m);
        hRawEnergyCh0->Fit("fgaus","","",lowlimsCh0[1],uplimsCh0[1]);
        m = fgaus->GetParameter(1);
        gPeaks0.push_back(m);
        hRawEnergyCh0->Fit("fgaus","","",lowlimsCh0[2],uplimsCh0[2]);
        m = fgaus->GetParameter(1);
        gPeaks0.push_back(m);
        
        //        TCanvas *cg = new TCanvas("cg","",700,600);
        //        gPad->SetLogy();
        //        hRawEnergyCh1->Draw();
    }
    
    TH1F *h0Calib = new TH1F("h0Calib","Energy vs peak; peak mean; energy (keV)",maxRawEnergy,0,maxRawEnergy);
    TH1F *h1Calib = new TH1F("h1Calib","Energy vs peak; peak mean; energy (keV)",maxRawEnergy,0,maxRawEnergy);
    int ipmax = (int)gPeaks0.size();
    if(filename == "run1_20160531_Na22_Teflon.txt" || filename == "run0_20160531_Na22_Teflon.txt" || filename == "run0_20160601_Na22_Al.txt") ipmax = 2;
    for(int ip = 0; ip < ipmax; ip++) {
        h0Calib->Fill(gPeaks0[ip],ePeak[ip]);
        h1Calib->Fill(gPeaks1[ip],ePeak[ip]);
    }
    TF1 *f0Calib = new TF1("f0Calib","pol1", 0,16400);
    h0Calib->Fit("f0Calib");
    TF1 *f1Calib = new TF1("f1Calib","pol1", 0,16400);
    h1Calib->Fit("f1Calib");
    
    Double_t e0Min = f0Calib->Eval(0);
    Double_t e0Max = f0Calib->Eval(16000);
    Double_t e1Min = f1Calib->Eval(0);
    Double_t e1Max = f1Calib->Eval(16000);
    
    TH1F *hCalibRawEnergyCh0 = new TH1F("hCalibRawEnergyCh0", "Raw energy ch0 - calibrated; keV",1000,e0Min,e0Max);
    TH1F *hCalibRawEnergyCh1 = new TH1F("hCalibRawEnergyCh1", "Raw energy ch1 - calibrated; keV",1000,e1Min,e1Max);
    
    
    if(filename == "run3_260516.txt") filename = "run3_260516";
    if(filename == "run4Rot_260516.txt") filename = "run4Rot_260516";
    if(filename == "run0_20160601_Na22_Al.txt") filename = "run0_20160601_Na22_Al";
    if(filename == "run0_20160531_Na22_Teflon.txt") filename = "run0_20160531_Na22_Teflon";
    if(filename == "run1_20160531_Na22_Teflon.txt") filename = "run1_20160531_Na22_Teflon";
    
    TFile f(Form("AnalysisResults_Calib_%s.root",filename.Data()),"recreate");
    TTree *tree = new TTree("t","");
    Double_t en0, en1, time0, time1;
    tree->Branch("time0",&time0,"time0/D");
    tree->Branch("en0",&en0,"en0/D");
    tree->Branch("time1",&time1,"time1/D");
    tree->Branch("en1",&en1,"en1/D");
    
    
    for(int i = 0; i < (int)TimeCalo1.size(); i++) {
        en0 = f0Calib->Eval(EnCalo0[i]);
        en1 = f1Calib->Eval(EnCalo1[i]);
        time0 = TimeCalo0[i];
        time1 = TimeCalo1[i];
        if(en0 == 16383 || en1 == 16383) continue;
        
        tree->Fill();
        hCalibRawEnergyCh0->Fill(en0);
        hCalibRawEnergyCh1->Fill(en1);
    }
    
    tree->Write();
    hCalibRawEnergyCh0->Write();
    hCalibRawEnergyCh1->Write();
    
    TCanvas *cc = new TCanvas("cc","",700,600);
    gPad->SetLogy();
    hCalibRawEnergyCh1->Draw();
    hCalibRawEnergyCh0->SetLineColor(kRed);
    hCalibRawEnergyCh0->Draw("same");
    //hRawEnergyCh1->Draw("same");
    return;
    
}

