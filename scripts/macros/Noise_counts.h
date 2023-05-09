/*設計思想
目的: ノイズ波形の抽出
アルゴリズム: 
1.電圧の閾値を超えたサンプリング点を抽出
2.閾値を超えたサンプリング店が3点以上の場合はノイズ波形として、電荷量を求める

*/

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "iostream"
#include "vector"


int NoiseCount(TString filesrc, TString treename, float thres, std::vector<float> &ave_wave)
{
    TFile* f  = TFile::Open(filesrc);
    TTree* tr = (TTree*)f->Get(treename);
    float data_time[1024];
    float wfp[1024] = {0};
    float wfn[1024] = {0};
    tr->SetBranchAddress("time", data_time);
    tr->SetBranchAddress("wform1",wfp);
    tr->SetBranchAddress("wform0",wfn);
    int nEve = tr->GetEntries();
    float crest_value;
    int noise_count = 0;

    for(int i=0;i<nEve;i++)
    {
        tr->GetEntry(i);
        for(int l=10;l<1024;l++)
        {
            if(l%1024==392)
	        {
	          wfp[l] = (wfp[l-1] + wfp[l+1])*0.5; 
	          wfn[l] = (wfn[l-1] + wfn[l+1])*0.5; 
	        }
            crest_value = (wfp[l]-wfn[l]) - ave_wave[l];

            if(crest_value > thres)
            {
                //std::cout << l << std::endl;
                noise_count += 1;
            }
        }
    }
    f->Close();
    
    return noise_count;

}



int NoiseAnalysis(TString filesrc, TString treename, TString wf_p, TString wf_n, float thres, std::vector<float> &ave_wave, std::vector<int> &event, std::vector<int> &cell, std::vector<float> &chrg)
{
    TFile* f  = TFile::Open(filesrc);
    TTree* tr = (TTree*)f->Get(treename);
    float data_time[1024];
    float wfp[1024] = {0};
    float wfn[1024] = {0};
    tr->SetBranchAddress("time", data_time);
    tr->SetBranchAddress(wf_p,wfp);
    tr->SetBranchAddress(wf_n,wfn);
    int nEve = tr->GetEntries();
    int state = 0;
    int start_cell;
    int end_cell; 
    float crest_value;
    float charge = 0;
    int noise_count = 0;

    for(int i=0;i<nEve;i++)
    {
        state = 0;
        charge = 0;
        start_cell = 0;
        end_cell = 0;
        tr->GetEntry(i);
        for(int l=10;l<1019;l++)
        {
            if(l%1024==392)
	        {
	          wfp[l] = (wfp[l-1] + wfp[l+1])*0.5; 
	          wfn[l] = (wfn[l-1] + wfn[l+1])*0.5; 
	        }
            crest_value = (wfp[l]-wfn[l]) - ave_wave[l];

            if(state == 0)
            {
                if(crest_value > thres)
                {
                    start_cell = l;
                    state = 1;
                
                }/*else
                {
                    continue;
                }
                */
                
            }
            //else if(state == 1)
            else
            {
                if(crest_value <= thres)
                {
                    end_cell = l-1;
                    if(end_cell - start_cell >= 2)
                    {
                        for(int k=start_cell-5;k<end_cell+5;k++)
                        {
                            charge += ((wfp[k]-wfn[k]) - ave_wave[k])*(data_time[k] - data_time[k-1]);
                        }
                        event.push_back(i);
                        cell.push_back(start_cell);
                        chrg.push_back(charge);
                        //std::cout << start_cell << std::endl;
                        //std::cout << end_cell << std::endl;
                        noise_count += 1;
                    }
                    state = 0;
                    start_cell = 0;
                    end_cell = 0;
                    charge = 0;
                }
                /*
                else{
                    continue;
                }
                */

            }
        }
    }
    f->Close();

    return noise_count;
}