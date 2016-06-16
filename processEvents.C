#include <iostream>
#include <sstream>
#include <fstream>
#include "TString.h"
#include "TLegend.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;

vector<double> TimeCalo0,TimeCalo1,EnCalo0,EnCalo1;
TString filenameData = "";

bool ReadRaws() {
    
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
            double ch,t,e;
            if (line[0] == '0') {
                iss>>ch>>t>>e;
                TimeCalo0.push_back(t);
                EnCalo0.push_back(e);
            }
            if (line[0] == '2') {
                iss>>ch>>t>>e;
                TimeCalo1.push_back(t);
                EnCalo1.push_back(e);
            }
        }
        printf("Unpacked calo data. Events for channel 0 = %lu, for channel 2 = %lu\n",TimeCalo0.size(),TimeCalo1.size());
    } else {
        printf("No file data found!\n");
        return false;
    }
    return true;
}
//______________________________________________________________

void BuildEvents() {
    
    int numbit = 14;
    int maxRawEnergy = pow(2,numbit)-1;
    
    TH1F *hDeltat = new TH1F("hDeltat","#Deltat Ch1-Ch0",20,-95,105);
    TH1F *hRawEnergyCh0 = new TH1F("hRawEnergyCh0","Raw Energy channel 0; channels; entries",maxRawEnergy,0,(float)maxRawEnergy);
    TH1F *hRawEnergyCh1 = new TH1F("hRawEnergyCh1","Raw Energy channel 1; channels; entries",maxRawEnergy,0,(float)maxRawEnergy);

    ofstream myfile;
    myfile.open (Form("processed_%s",filenameData.Data()));
    
    Long_t startIndex = 0;
    for(Long_t iev0 = 0; iev0 < (Long_t)TimeCalo0.size(); iev0++) {
        
        if((iev0%5000)==0) printf("event %ld\n",iev0);
 
        double time0 = TimeCalo0[iev0];
        if(iev0 > 0)
            if(time0 < TimeCalo0[iev0-1] || time0 > TimeCalo0[iev0+1]) continue;
        double tmpdt = TMath::Power(10,12);
        double time1, dt;
        Long_t indexCh0 = iev0;
        Long_t indexCh1;
        for(Long_t iev1 = startIndex; iev1 < (Long_t)TimeCalo1.size(); iev1++) {
            time1 = TimeCalo1[iev1];
            if(iev1 > 0)
                if(time1 < TimeCalo1[iev1-1] || time1 > TimeCalo1[iev1+1]) continue;
            dt = TMath::Abs(time1 - time0);
            if(dt > tmpdt) break;
            else {
                tmpdt = dt;
                indexCh1 = iev1;
            }
        }
        startIndex = indexCh1 -10;
        if(startIndex < 0) startIndex = 0;
        if(EnCalo1[indexCh1] == 16383 || EnCalo0[indexCh0] == 16383) continue;
        myfile <<Form("%.0f",TimeCalo0[indexCh0])<<" "<<Form("%.0f",EnCalo0[indexCh0])<<" "<<Form("%.0f",TimeCalo1[indexCh1])<<" "<<Form("%.0f",EnCalo1[indexCh1])<<endl;
        hDeltat->Fill(10.*(TimeCalo1[indexCh1]-TimeCalo0[indexCh0]));
        if(TMath::Abs(10.*(TimeCalo1[indexCh1]-TimeCalo0[indexCh0]))==470) printf("%.0f, %.0f, %.0f, %.0f\n",TimeCalo0[indexCh0],EnCalo0[indexCh0],TimeCalo1[indexCh1],EnCalo1[indexCh1]);
        hRawEnergyCh0->Fill(EnCalo0[indexCh0]);
        hRawEnergyCh1->Fill(EnCalo1[indexCh1]);
        if(indexCh1 == (Long_t)(TimeCalo1.size()-2))
            break;
    }
    myfile.close();
    
    TCanvas *c = new TCanvas("c","",1000,600);
    c->Divide(2);
    c->cd(1);
    gPad->SetLogy();
    hDeltat->DrawCopy();
    c->cd(2);
    hRawEnergyCh0->SetLineColor(kBlue);
    hRawEnergyCh1->SetLineColor(kRed);
    hRawEnergyCh0->DrawCopy();
    hRawEnergyCh1->DrawCopy("same");
    
    TLegend *legend = new TLegend(0.45,0.7,0.65,0.85);
    legend->AddEntry(hRawEnergyCh0,"Canale 0");
    legend->AddEntry(hRawEnergyCh1,"Canale 1");
    legend->SetFillColor(0);
    legend->Draw("same");
}

//______________________________________________________________

void processEvents(TString filename = "run3_260516.txt") {
    
    filenameData = filename;
    bool raws = ReadRaws();
    if(raws) {
        BuildEvents();
    }
    return;
}
