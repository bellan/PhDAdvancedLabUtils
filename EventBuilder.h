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
    
    TString fFileName;
    TString fCalibFileName;
    TString fdir;
    
    ClassDef(EventBuilder,1);

};

#endif
