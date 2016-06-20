#ifndef EVENTBUILDER_H
#define EVENTBUILDER_H

#include <TObject.h>
#include <iostream>
#include "TString.h"
#include "TF1.h"

using namespace std;

class EventBuilder : public TObject {
public:
    
    EventBuilder();
    virtual ~EventBuilder();

    
    void SetDirectory(TString fname)       {fdir = fname;}
    void SetRawsFile(TString fname)        {fFileName = fname;}
    void SetCalibrationFile(TString fname) {fCalibFileName = fname;}
    void SetFitRangeForCalibration(Double_t lowlimsCh0[], Double_t lowlimsCh1[], Double_t uplimsCh0[], Double_t uplimsCh1[]);

    TString GetDirectory()       {return fdir;}
    TString GetRawsFile()        {return fFileName;}
    TString GetCalibrationFile() {return fCalibFileName;}

    bool ReadRaws();
    void BuildEvents();
    void SetCalibrationCurve();
    void DoCalibration();
    
 private:

    vector<double> fTimeCalo0; //unmatched
    vector<double> fTimeCalo1;
    vector<double> fEnCalo0;
    vector<double> fEnCalo1;
    
    vector<double> fEvTime0;  //processed
    vector<double> fEvTime1;
    vector<double> fEvEn0;
    vector<double> fEvEn1;

    TF1 *fCalibFuncChan0;
    TF1 *fCalibFuncChan1;

    Double_t fLowLimsCh0[2];
    Double_t fLowLimsCh1[2];
    Double_t fUpLimsCh0[2];
    Double_t fUpLimsCh1[2];

    
    TString fFileName;
    TString fCalibFileName;
    TString fdir;
    
    ClassDef(EventBuilder,1);

};

#endif
