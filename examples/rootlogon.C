/*
   rootlogon.C
   This script is automatically invoked whenever ROOT is started.
   Add session-level configurations as needed.
*/

void rootlogon()
{
  Printf("Starting ROOT version %s.", gROOT->GetVersion());
  Printf("Running %s/rootlogon.C on %s.",
         gSystem->Getenv("PWD"),
         gSystem->HostName());

  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));

  // Include paths
  // =============
  // There are 2 separate paths: one for ACLiC, and one for CINT or CLING.
  // 1. ACLiC include path
  const char *pwd = gSystem->WorkingDirectory();
  const char *inc = gSystem->DirName(pwd);
  gSystem->AddIncludePath(Form("-I%s ", inc));

  // 2. Interpreter include path
#ifdef __CINT__
  gROOT->ProcessLine(Form(".include %s", inc));
#endif
#ifdef __CLING__
  gROOT->ProcessLine(Form(".I %s", inc));
#endif

  gROOT->LoadMacro("../UnfoldingUtils.C+");
  gROOT->LoadMacro("../TestProblems.C+");

}