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

    void SetSource(TString source)                      {fSource = source;}
    void SetNChannels(Int_t nChann)                     {fNChannels = nChann;}
    void SetDirectory(TString fname)                    {fdir = fname;}
    void SetRawsFile(TString fname)                     {fFileName = fname;}
    void SetCalibrationFile(TString fname)              {fCalibFileName = fname;}
    void SetChanID(string id0, string id1, string id2)  {fCh0Id = id0; fCh1Id = id1; fCh2Id = id2;}
    void SetFitRangeForCalibration(Double_t lowlimsCh0[], Double_t uplimsCh0[], Double_t lowlimsCh1[], Double_t uplimsCh1[], Double_t lowlimsCh2[], Double_t uplimsCh2[]);
    void CheckTimeAlignment();
    
    Int_t   GetNChannels()       {return fNChannels;}
    TString GetSource()          {return fSource;}
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
    vector<double> fTimeCalo2;
    vector<double> fEnCalo0;
    vector<double> fEnCalo1;
    vector<double> fEnCalo2;
    
    vector<double> fEvTime0;  //processed
    vector<double> fEvTime1;
    vector<double> fEvTime2;
    vector<double> fEvEn0;
    vector<double> fEvEn1;
    vector<double> fEvEn2;

    TF1 *fCalibFuncChan0;
    TF1 *fCalibFuncChan1;
    TF1 *fCalibFuncChan2;

    Double_t fLowLimsCh0[2];
    Double_t fLowLimsCh1[2];
    Double_t fLowLimsCh2[2];
    Double_t fUpLimsCh0[2];
    Double_t fUpLimsCh1[2];
    Double_t fUpLimsCh2[2];

    Int_t fNChannels;
    
    TString fSource;
    TString fFileName;
    TString fCalibFileName;
    TString fdir;
    
    string fCh0Id;
    string fCh1Id;
    string fCh2Id;
    
    ClassDef(EventBuilder,1);

};

#endif
