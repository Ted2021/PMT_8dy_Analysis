#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "iostream"
#include "vector"

void Average_Make(TString filesrc, TString treename, TString wavep, TString waven, std::vector<float> &par, std::vector<float> &par2)
{
    TFile* f  = TFile::Open(filesrc);
    TTree* tr = (TTree*)f->Get(treename);
    //float time_e[1024] = {0};
    float data_time[1024];
    float wfp[1024] = {0};
    float wfn[1024] = {0};
    tr->SetBranchAddress("time", data_time);
    tr->SetBranchAddress(wavep,wfp);
    tr->SetBranchAddress(waven,wfn);
    int nEve = tr->GetEntries();

    float ave_wave[1024] = {0};
    
    for(int i=0;i<nEve;i++)
    {
        tr->GetEntry(i);
        for(int l=0;l<1024;l++)
        {
            if(l%1024==392)
	        {
	          wfp[l] = (wfp[l-1] + wfp[l+1])*0.5; 
	          wfn[l] = (wfn[l-1] + wfn[l+1])*0.5; 
	        }
            ave_wave[l] += (wfp[l]-wfn[l])/nEve;
        }
    }
    f->Close();
    //パラメータのReturn
    for(int i=0;i<1024;i++)
    {
        par.push_back(ave_wave[i]);
        par2.push_back(data_time[i]);
    }
}
