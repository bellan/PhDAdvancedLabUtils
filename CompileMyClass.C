void CompileMyClass() {
    
    gSystem->CompileMacro("EventBuilder.cxx","kfg");
    gSystem->CompileMacro("buildEvent.C","kfg");
    buildEvent();
    
}
