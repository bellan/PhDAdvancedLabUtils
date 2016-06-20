void CompileMyClass(string dir) {
    
    gSystem->CompileMacro("EventBuilder.cxx","kfg");
    gSystem->CompileMacro("buildEvent.C","kfg");
    buildEvent(dir);
    
}
