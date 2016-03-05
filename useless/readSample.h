#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <TFile.h>
#include <TSystem.h>
#include <TSystemDirectory.h>

void readSample(std::string inputFile, std::vector<string>& infiles){

  if( inputFile.find(".root") != std::string::npos )  
    infiles.push_back(inputFile.data());

  else{

    string inputTextFile = "inputdir.txt";
    gSystem->Exec(Form("rm -rf %s",inputTextFile.data()));
    gSystem->Exec(Form("ls -R %s | grep -a \"%s\" >> %s", inputFile.data(),"data7", inputTextFile.data()));
    TSystemDirectory *base = new TSystemDirectory("root","root"); 
    int nfile = 0;
    string tempdir;
    ifstream fin;
    fin.open(inputTextFile.data());
    fin >> tempdir;
    TString realDirName = gSystem->GetFromPipe(Form("a=%s; echo ${a%%:*}",tempdir.data()));
    while( !fin.eof() ){
      if( realDirName.Contains("fail") ){
	fin >> tempdir;
	realDirName = gSystem->GetFromPipe(Form("a=%s; echo ${a%%:*}",tempdir.data()));      
	continue;
      }
      base->SetDirectory(realDirName.Data());
      TList *listOfFiles = base->GetListOfFiles();
      TIter fileIt(listOfFiles);
      TFile *fileH = new TFile();
      while( (fileH = (TFile*)fileIt()) ) {
	std::string fileN = fileH->GetName();
	std::string baseString = "root";
	if( fileH->IsFolder() ) continue;
	if( fileN.find(baseString) == std::string::npos ) continue;
	nfile++;
	string tempfile = Form("%s/%s",realDirName.Data(),fileN.data());
	infiles.push_back(tempfile);
      }
      fin >> tempdir;
      realDirName = gSystem->GetFromPipe(Form("a=%s; echo ${a%%:*}",tempdir.data()));
    }

  }

}
