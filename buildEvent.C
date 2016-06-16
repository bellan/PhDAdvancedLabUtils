#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TString.h>
#include <TSystem.h>
#include "EventBuilder.h"

using namespace std;


void buildEvent() {
    
    TString dir = "20160613_2cm/";
    TString filename;
    TString calibFilename = "run_20160613_calibration_2cm";
    
    const int nruns = 6;
    int suff[nruns] = {650,675,700,725,750,800};
    
    EventBuilder *event = new EventBuilder();
    
    //set the calibration
    event->SetDirectory(dir.Data());
    event->SetRawsFile(calibFilename.Data());
    event->SetCalibrationFile(calibFilename.Data());
    bool raws = event->ReadRaws();
    if(!raws) {
        printf("raw file not found!\n");
        return;
    }
    printf("Set the calibration file %s.txt \n",calibFilename.Data());
    event->BuildEvents();
    event->SetCalibrationCurve();
    event->DoCalibration();
   
    
    //process all the runs
    for(int irun = 0; irun < nruns; irun++) {
        filename = Form("run_20160613_%d_0_2cm",suff[irun]);
        printf("Building file %s%s.txt \n",dir.Data(),filename.Data());
        event->SetRawsFile(filename.Data());
        raws = event->ReadRaws();
        if(!raws) {
            printf("raw file %s not found!\n",filename.Data());
            return;
        }
        event->BuildEvents();
        event->DoCalibration();
    }
    
}

