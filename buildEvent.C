#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TString.h>
#include <TSystem.h>
#include "EventBuilder.h"

using namespace std;
Double_t lowlimsCh0[2], lowlimsCh1[2],  lowlimsCh2[2], uplimsCh0[2], uplimsCh1[2], uplimsCh2[2];
Int_t fNChannels;
string ch0ID, ch1ID, ch2ID;
TString source;

//////////////////////////////////////////////////////////////////////////////////////////////
// Here we build the events.                                                                //
// buildEvent() needs as input a list "runs.txt" of the files of run inside the proper dir. //
// Set the limits for the fits used to do calibration in SetLimitsForCalibFit().            //
//////////////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________________
void SetLimitsForCalibFit(string run) {
    
    if(run == "run_20160613_calibrazione_4.txt") {
        lowlimsCh0[0] = 2850;
        lowlimsCh0[1] = 7100;
        uplimsCh0[0]  = 3300;
        uplimsCh0[1]  = 7700;
        
        lowlimsCh1[0] = 2700;
        lowlimsCh1[1] = 6700;
        uplimsCh1[0]  = 3000;
        uplimsCh1[1]  = 7200;
        
        lowlimsCh2[0] = 0;
        lowlimsCh2[1] = 0;
        uplimsCh2[0]  = 0;
        uplimsCh2[1]  = 0;
        
        source = "Na22";
        fNChannels = 2;
        ch0ID = "0";
        ch1ID = "2";
    }
    if(run == "run_20160613_calibration_2cm.txt") {
        lowlimsCh0[0] = 2850;
        lowlimsCh0[1] = 7100;
        uplimsCh0[0]  = 3300;
        uplimsCh0[1]  = 7700;
        
        lowlimsCh1[0] = 2700;
        lowlimsCh1[1] = 6700;
        uplimsCh1[0]  = 3000;
        uplimsCh1[1]  = 7200;
        
        lowlimsCh2[0] = 0;
        lowlimsCh2[1] = 0;
        uplimsCh2[0]  = 0;
        uplimsCh2[1]  = 0;
        
        source = "Na22";
        fNChannels = 2;
        ch0ID = "0";
        ch1ID = "2";
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
        
        lowlimsCh2[0] = 0;
        lowlimsCh2[1] = 0;
        uplimsCh2[0]  = 0;
        uplimsCh2[1]  = 0;
        
        source = "Na22";
        fNChannels = 2;
        
    }
    else if(run == "run_20160630_Na22_teflon_tripla.txt") {
        lowlimsCh0[0] = 2500;
        lowlimsCh0[1] = 6300;
        uplimsCh0[0]  = 2900;
        uplimsCh0[1]  = 6800;
        
        lowlimsCh1[0] = 2850;
        lowlimsCh1[1] = 7000;
        uplimsCh1[0]  = 3200;
        uplimsCh1[1]  = 7500;
        
        lowlimsCh2[0] = 2700;
        lowlimsCh2[1] = 6700;
        uplimsCh2[0]  = 3000;
        uplimsCh2[1]  = 7200;
        
        fNChannels = 3;
        ch0ID = "0";
        ch1ID = "1";
        ch2ID = "2";
        source = "Na22";
        
    }
    else if(run == "run_20160629_Ba.txt") {
        lowlimsCh0[0] = 400;
        lowlimsCh0[1] = 1720;
        uplimsCh0[0]  = 525;
        uplimsCh0[1]  = 2100;
        
        lowlimsCh1[0] = 450;
        lowlimsCh1[1] = 1880;
        uplimsCh1[0]  = 560;
        uplimsCh1[1]  = 2120;
        
        lowlimsCh2[0] = 0;
        lowlimsCh2[1] = 0;
        uplimsCh2[0]  = 0;
        uplimsCh2[1]  = 0;
        
        fNChannels = 2;
        ch0ID = "0";
        ch1ID = "2";
        source = "Ba133";
        
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
        
        lowlimsCh2[0] = 0;
        lowlimsCh2[1] = 0;
        uplimsCh2[0]  = 0;
        uplimsCh2[1]  = 0;
        
        fNChannels = 2;
    }
}

//_____________________________________________________________________________________
//void buildEvent(string dir = "../Na22_20160629_tripla/"/*"../20160616/"*/) {
//void buildEvent(string dir = "../20160629_Ba/") {
void buildEvent(string dir = "../20160613_2cm/") {
    
    EventBuilder *event;
    bool raws;
    string run;
    ifstream myfile (Form("%sruns.txt",dir.data()));
    string lastrun;
    
    //set the calibration
    bool calibDone = false;
    if (myfile.is_open()) {
        event = new EventBuilder();
        event->SetDirectory(dir.data());
        while ( getline (myfile,run) ) {
            std::string search;
            search = "calibration";
            size_t pos;
            pos = run.find(search);
            if (pos != std::string::npos) {
                SetLimitsForCalibFit(run);
                event->SetRawsFile(run.data());
                event->SetCalibrationFile(run.data());
                event->SetFitRangeForCalibration(lowlimsCh0,uplimsCh0,lowlimsCh1,uplimsCh1,lowlimsCh2,uplimsCh2);
                event->SetSource(source);
                event->SetNChannels(fNChannels);
                event->SetChanID(ch0ID,ch1ID,ch2ID);
                raws = event->ReadRaws();
                if(!raws) {
                    printf("raw file not found!\n");
                    return;
                }
                printf("Set the calibration file %s \n",run.data());
                event->BuildEvents();
                event->SetCalibrationCurve();
                event->DoCalibration();
                calibDone = true;
            }
            lastrun = run;
        }
        if(!calibDone) {
            SetLimitsForCalibFit(lastrun);
            event->SetRawsFile(lastrun.data());
            event->SetCalibrationFile(lastrun.data());
            event->SetFitRangeForCalibration(lowlimsCh0,uplimsCh0,lowlimsCh1,uplimsCh1,lowlimsCh2,uplimsCh2);
            event->SetSource(source);
            event->SetNChannels(fNChannels);
            event->SetChanID(ch0ID,ch1ID,ch2ID);
            raws = event->ReadRaws();
            if(!raws) {
                printf("raw file not found!\n");
                return;
            }
            printf("Set the calibration file %s \n",lastrun.data());
            event->BuildEvents();
            event->SetCalibrationCurve();
            event->DoCalibration();
            calibDone = true;
            
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
                    event->CheckTimeAlignment();
                }
            }
            file.close();
    
        }
    
}


