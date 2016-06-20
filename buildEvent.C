#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TString.h>
#include <TSystem.h>
#include "EventBuilder.h"

using namespace std;
Double_t lowlimsCh0[2], lowlimsCh1[2], uplimsCh0[2], uplimsCh1[2];


//////////////////////////////////////////////////////////////////////////////////////////////
// Here we build the events.                                                                //
// buildEvent() needs as input a list "runs.txt" of the files of run inside the proper dir. //
// Set the limits for the fits used to do calibration in SetLimitsForCalibFit().            //
//////////////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________________
void SetLimitsForCalibFit(string run) {
    
    if(run == "run_20160613_calibration_2cm.txt") {
        lowlimsCh0[0] = 2850;
        lowlimsCh0[1] = 7100;
        uplimsCh0[0]  = 3300;
        uplimsCh0[1]  = 7700;
        
        lowlimsCh1[0] = 2700;
        lowlimsCh1[1] = 6700;
        uplimsCh1[0]  = 3000;
        uplimsCh1[1]  = 7200;
    }
    else if(run == "run_20160616_2Na22_calibration.txt") {
        lowlimsCh0[0] = 2600;
        lowlimsCh0[1] = 6400;
        uplimsCh0[0]  = 2850;
        uplimsCh0[1]  = 6900;
        
        lowlimsCh1[0] = 2750;
        lowlimsCh1[1] = 6800;
        uplimsCh1[0]  = 3000;
        uplimsCh1[1]  = 7250;
    }
    else {
        lowlimsCh0[0] = 2600;
        lowlimsCh0[1] = 6400;
        uplimsCh0[0]  = 2850;
        uplimsCh0[1]  = 6900;
        
        lowlimsCh1[0] = 2750;
        lowlimsCh1[1] = 6800;
        uplimsCh1[0]  = 3000;
        uplimsCh1[1]  = 7250;
    }
}

//_____________________________________________________________________________________
void buildEvent(string dir = "../20160616/") {
    
    EventBuilder *event;
    bool raws;
    string run;
    ifstream myfile (Form("%sruns.txt",dir.data()));

    //set the calibration
    if (myfile.is_open()) {
        while ( getline (myfile,run) ) {
            std::string search;
            search = "calibration";
            size_t pos;
            pos = run.find(search);
            if (pos != std::string::npos) {
                event = new EventBuilder();
                SetLimitsForCalibFit(run);
                event->SetDirectory(dir.data());
                event->SetRawsFile(run.data());
                event->SetCalibrationFile(run.data());
                event->SetFitRangeForCalibration(lowlimsCh0,lowlimsCh1,uplimsCh0,uplimsCh1);
                raws = event->ReadRaws();
                if(!raws) {
                    printf("raw file not found!\n");
                    return;
                }
                printf("Set the calibration file %s \n",run.data());
                event->BuildEvents();
                event->SetCalibrationCurve();
                event->DoCalibration();
            }
        }
        myfile.close();
    } else {
        printf("Unable to open file\n");
        return;
    }
    
    
    //process all the runs
    ifstream file (Form("%sruns.txt",dir.data()));
    if (file.is_open()) {
        while ( getline (file,run) ) {
            std::string search;
            search = "alib";
            size_t pos;
            pos = run.find(search);
            if (pos == std::string::npos) {
                printf("Building file %s%s \n",dir.data(),run.data());
                event->SetRawsFile(run.data());
                raws = event->ReadRaws();
                if(!raws) {
                    printf("raw file %s not found!\n",run.data());
                    return;
                }
                event->BuildEvents();
                event->DoCalibration();
            }
        }
        file.close();

    }
    
}


