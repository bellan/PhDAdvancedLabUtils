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

void mergeCalibFiles() {
    
    TString filename0 = "20160613_2cm/run_20160613_CalibCh0_2cm.txt";
    TString filename1 = "20160613_2cm/run_20160613_CalibCh2_2cm.txt";
    TString fout = "20160613_2cm/run_20160613_calibration_2cm.txt";
    
    TimeCalo0.clear();
    TimeCalo1.clear();
    EnCalo0.clear();
    EnCalo1.clear();
    
    ifstream inFileData0(filename0.Data());
    if(inFileData0) {
        string line;
        while(getline(inFileData0, line)) {
            if (line[0] == '#') continue;
            std::istringstream iss(line);
            double ch,t,e;
            if (line[0] == '0') {
                iss>>ch>>t>>e;
                TimeCalo0.push_back(t);
                EnCalo0.push_back(e);
            }
          }
        printf("Unpacked calo data. Events for channel 0 = %lu\n",TimeCalo0.size());
    } else {
        printf("No file data found!\n");
        return;
    }

    ifstream inFileData1(filename1.Data());
    if(inFileData1) {
        string line;
        while(getline(inFileData1, line)) {
            if (line[0] == '#') continue;
            std::istringstream iss(line);
            double ch,t,e;
            if (line[0] == '2') {
                iss>>ch>>t>>e;
                TimeCalo1.push_back(t);
                EnCalo1.push_back(e);
            }
        }
        printf("Unpacked calo data. Events for channel 2 = %lu\n",TimeCalo1.size());
    } else {
        printf("No file data found!\n");
        return;
    }

  
    ofstream myfile;
    myfile.open(fout.Data());
    for(Long_t iev0 = 0; iev0 < (Long_t)TimeCalo0.size(); iev0++) {
        myfile<<"0 "<<Form("%.0f",TimeCalo0[iev0])<<" "<<Form("%.0f",EnCalo0[iev0])<<endl;
    }
    for(Long_t iev1 = 0; iev1 < (Long_t)TimeCalo1.size(); iev1++) {
        myfile<<"2 "<<Form("%.0f",TimeCalo1[iev1])<<" "<<Form("%.0f",EnCalo1[iev1])<<endl;
    }
    myfile.close();

    return;

    
    
    
}
