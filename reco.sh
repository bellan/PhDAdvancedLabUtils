#!/bin/bash

#root -l <<EOF
#.x mergeCalibFiles.C () 
#.q
#EOF

#find . -name "run_*.txt"  -exec basename \{} .po \; > runs.txt

root -l <<EOF
.x CompileMyClass.C ("../20160616/") 
.q
EOF


#find . -name "AnalysisResults_Calibrated*_0.root"  -exec basename \{} .po \; > calibRuns0.txt
#find . -name "AnalysisResults_Calibrated*_120.root"  -exec basename \{} .po \; > calibRuns120.txt
#find . -name "AnalysisResults_Calibrated*_240.root"  -exec basename \{} .po \; > calibRuns240.txt
#
#


