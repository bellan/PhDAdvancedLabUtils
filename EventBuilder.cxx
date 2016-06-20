#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "EventBuilder.h"

using namespace std;

ClassImp(EventBuilder);

//_______________________________________________

EventBuilder::EventBuilder():
fTimeCalo0(0),
fTimeCalo1(0),
fEnCalo0(0),
fEnCalo1(0),
fEvTime0(0),
fEvTime1(0),
fEvEn0(0),
fEvEn1(0),
fCalibFileName(0),
fCalibFuncChan0(0x0),
fCalibFuncChan1(0x0) {
    
    fFileName = "";
    fCalibFileName = "";
    fdir = "";
    
}

//_______________________________________________

EventBuilder::~EventBuilder() {
    
    //Default destructor
    
    delete fCalibFuncChan0;
    delete fCalibFuncChan1;
}
//_______________________________________________
void EventBuilder::SetFitRangeForCalibration(Double_t lowlimsCh0[], Double_t lowlimsCh1[], Double_t uplimsCh0[], Double_t uplimsCh1[]) {
    
    for(int i = 0; i < 2; i++) {
        fLowLimsCh0[i] = lowlimsCh0[i];
        fLowLimsCh1[i] = lowlimsCh1[i];
        fUpLimsCh0[i]  = uplimsCh0[i];
        fUpLimsCh1[i]  = uplimsCh1[i];
    }
}

//_______________________________________________
bool EventBuilder::ReadRaws() {
    
    fTimeCalo0.clear();
    fTimeCalo1.clear();
    fEnCalo0.clear();
    fEnCalo1.clear();
    
    fEvTime0.clear();
    fEvTime1.clear();
    fEvEn0.clear();
    fEvEn1.clear();
    
    ifstream inFileData(Form("%s%s",fdir.Data(),fFileName.Data()));
    if(inFileData) {
        string line;
        while(getline(inFileData, line)) {
            if (line[0] == '#') continue;
            std::istringstream iss(line);
            double ch,t,e;
            if (line[0] == '0') {
                iss>>ch>>t>>e;
                fTimeCalo0.push_back(t);
                fEnCalo0.push_back(e);
            }
            if (line[0] == '2') {
                iss>>ch>>t>>e;
                fTimeCalo1.push_back(t);
                fEnCalo1.push_back(e);
            }
        }
        printf("Unpacked calo data. Events for channel 0 = %lu, for channel 2 = %lu\n",fTimeCalo0.size(),fTimeCalo1.size());
    } else {
        printf("No file data %s found!\n",Form("%s%s",fdir.Data(),fFileName.Data()));
        return false;
    }
    return true;
    
}

//_______________________________________________
void EventBuilder::BuildEvents() {
    
    int numbit = 14;
    int maxRawEnergy = pow(2,numbit)-1;
    
    //    TH1F *hDeltat = new TH1F("hDeltat","#Deltat Ch1-Ch0",20,-95,105);
    //    TH1F *hRawEnergyCh0 = new TH1F("hRawEnergyCh0","Raw Energy channel 0; channels; entries",maxRawEnergy,0,(float)maxRawEnergy);
    //    TH1F *hRawEnergyCh1 = new TH1F("hRawEnergyCh1","Raw Energy channel 1; channels; entries",maxRawEnergy,0,(float)maxRawEnergy);
    
    //    ofstream myfile;
    //    myfile.open (Form("%sprocessed_%s.txt",fdir.Data(),fFileName.Data()));
    
    Long_t startIndex = 0;
    for(Long_t iev0 = 0; iev0 < (Long_t)fTimeCalo0.size(); iev0++) {
        
        if((iev0%5000)==0) printf("event %ld\n",iev0);
        
        double time0 = fTimeCalo0[iev0];
        if(iev0 > 0)
            if(time0 < fTimeCalo0[iev0-1] || time0 > fTimeCalo0[iev0+1]) continue;
        double tmpdt = TMath::Power(10,12);
        double time1, dt;
        Long_t indexCh0 = iev0;
        Long_t indexCh1;
        for(Long_t iev1 = startIndex; iev1 < (Long_t)fTimeCalo1.size(); iev1++) {
            time1 = fTimeCalo1[iev1];
            if(iev1 > 0)
                if(time1 < fTimeCalo1[iev1-1] || time1 > fTimeCalo1[iev1+1]) continue;
            dt = TMath::Abs(time1 - time0);
            if(dt > tmpdt) break;
            else {
                tmpdt = dt;
                indexCh1 = iev1;
            }
        }
        startIndex = indexCh1 -10;
        if(startIndex < 0) startIndex = 0;
        if(fEnCalo1[indexCh1] == 16383 || fEnCalo0[indexCh0] == 16383) continue;
        
        fEvTime0.push_back(fTimeCalo0[indexCh0]);
        fEvTime1.push_back(fTimeCalo1[indexCh1]);
        fEvEn0.push_back(fEnCalo0[indexCh0]);
        fEvEn1.push_back(fEnCalo1[indexCh1]);
        
        
        
        
        //    myfile <<Form("%.0f",fTimeCalo0[indexCh0])<<" "<<Form("%.0f",fEnCalo0[indexCh0])<<" "<<Form("%.0f",fTimeCalo1[indexCh1])<<" "<<Form("%.0f",fEnCalo1[indexCh1])<<endl;
        //        hDeltat->Fill(10.*(fTimeCalo1[indexCh1]-fTimeCalo0[indexCh0]));
        //        hRawEnergyCh0->Fill(fEnCalo0[indexCh0]);
        //        hRawEnergyCh1->Fill(fEnCalo1[indexCh1]);
        if(indexCh1 == (Long_t)(fTimeCalo1.size()-2))
            break;
    }
    // myfile.close();
    
    //    TCanvas *c = new TCanvas("c","",1000,600);
    //    c->Divide(2);
    //    c->cd(1);
    //    gPad->SetLogy();
    //    hDeltat->DrawCopy();
    //    c->cd(2);
    //    hRawEnergyCh0->SetLineColor(kBlue);
    //    hRawEnergyCh1->SetLineColor(kRed);
    //    hRawEnergyCh0->DrawCopy();
    //    hRawEnergyCh1->DrawCopy("same");
    //
    //    TLegend *legend = new TLegend(0.45,0.7,0.65,0.85);
    //    legend->AddEntry(hRawEnergyCh0,"Canale 0");
    //    legend->AddEntry(hRawEnergyCh1,"Canale 1");
    //    legend->SetFillColor(0);
    //    legend->Draw("same");
    
    
    
}

//_______________________________________________
void EventBuilder::SetCalibrationCurve() {
    
    
    int numbit = 14;
    int maxRawEnergy = pow(2,numbit)-1;
    
    Double_t ePeak[3] = {511,1274,1785};
    vector<Double_t> gPeaks0, gPeaks1;
    
    TH1F *hRawEnergyCh0, *hRawEnergyCh1;
    
    hRawEnergyCh0 = new TH1F("hRawEnergyCh0","; channels; entries",512,0,maxRawEnergy);
    hRawEnergyCh1 = new TH1F("hRawEnergyCh1","; channels; entries",512,0,maxRawEnergy);
    hRawEnergyCh0->SetLineColor(kRed);
    hRawEnergyCh1->SetLineColor(kBlue);
    hRawEnergyCh0->SetStats(0);
    hRawEnergyCh1->SetStats(0);
    
    for(int i = 0; i < (int)fEvTime1.size(); i++) {
        hRawEnergyCh0->Fill(fEvEn0[i]);
        hRawEnergyCh1->Fill(fEvEn1[i]);
    }
    
    TF1 *fgaus = new TF1("fgaus","gaus", 0,16000);
    
    //channel 1
    hRawEnergyCh1->Fit("fgaus","","",fLowLimsCh1[0],fUpLimsCh1[0]);
    Double_t m = fgaus->GetParameter(1);
    gPeaks1.push_back(m);
    hRawEnergyCh1->Fit("fgaus","","",fLowLimsCh1[1],fUpLimsCh1[1]);
    m = fgaus->GetParameter(1);
    gPeaks1.push_back(m);
    
    //channel 0
    hRawEnergyCh0->Fit("fgaus","","",fLowLimsCh0[0],fUpLimsCh0[0]);
    m = fgaus->GetParameter(1);
    gPeaks0.push_back(m);
    hRawEnergyCh0->Fit("fgaus","","",fLowLimsCh0[1],fUpLimsCh0[1]);
    m = fgaus->GetParameter(1);
    gPeaks0.push_back(m);
    
    TH1F *h0Calib = new TH1F("h0Calib","Energy vs peak; peak mean; energy (keV)",maxRawEnergy,0,maxRawEnergy);
    TH1F *h1Calib = new TH1F("h1Calib","Energy vs peak; peak mean; energy (keV)",maxRawEnergy,0,maxRawEnergy);
    int ipmax = 2;
    for(int ip = 0; ip < ipmax; ip++) {
        h0Calib->Fill(gPeaks0[ip],ePeak[ip]);
        h1Calib->Fill(gPeaks1[ip],ePeak[ip]);
    }
    fCalibFuncChan0 = new TF1("fCalibFuncChan0","pol1", 0,16400);
    fCalibFuncChan1 = new TF1("fCalibFuncChan1","pol1", 0,16400);
    
    h0Calib->Fit("fCalibFuncChan0");
    h1Calib->Fit("fCalibFuncChan1");
    
    
}

//_______________________________________________
void EventBuilder::DoCalibration() {
    
    int numbit = 14;
    int maxRawEnergy = pow(2,numbit)-1;
    
    Double_t e0Min = fCalibFuncChan0->Eval(0);
    Double_t e0Max = fCalibFuncChan0->Eval(16000);
    Double_t e1Min = fCalibFuncChan1->Eval(0);
    Double_t e1Max = fCalibFuncChan1->Eval(16000);
    
    //    TH1F *hCalibRawEnergyCh0 = new TH1F("hCalibRawEnergyCh0", "Raw energy ch0 - calibrated; keV",1000,e0Min,e0Max);
    //    TH1F *hCalibRawEnergyCh1 = new TH1F("hCalibRawEnergyCh1", "Raw energy ch1 - calibrated; keV",1000,e1Min,e1Max);
    
    string newname = (string)fFileName;
    newname = newname.substr(0, newname.find("."));
    
    TFile f(Form("%sAnalysisResults_Calibrated_%s.root",fdir.Data(),newname.data()),"recreate");
    TTree *tree = new TTree("t","");
    Double_t en0, en1, time0, time1;
    tree->Branch("time0",&time0,"time0/D");
    tree->Branch("en0",&en0,"en0/D");
    tree->Branch("time1",&time1,"time1/D");
    tree->Branch("en1",&en1,"en1/D");
    
    
    for(int i = 0; i < (int)fEvTime1.size(); i++) {
        en0 = fCalibFuncChan0->Eval(fEvEn0[i]);
        en1 = fCalibFuncChan1->Eval(fEvEn1[i]);
        time0 = fEvTime0[i];
        time1 = fEvTime1[i];
        tree->Fill();
        //        hCalibRawEnergyCh0->Fill(en0);
        //        hCalibRawEnergyCh1->Fill(en1);
    }
    
    tree->Write();
    //    hCalibRawEnergyCh0->Write();
    //    hCalibRawEnergyCh1->Write();
    
    //    TCanvas *cc = new TCanvas("cc","",700,600);
    //    gPad->SetLogy();
    //    hCalibRawEnergyCh1->Draw();
    //    hCalibRawEnergyCh0->SetLineColor(kRed);
    //    hCalibRawEnergyCh0->Draw("same");
    //hRawEnergyCh1->Draw("same");
    
}















