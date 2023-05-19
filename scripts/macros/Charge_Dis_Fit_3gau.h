#include "vector"
#include "TH1F.h"
#include "TF1.h"
void Fit_Charge_Dist(TH1F* hst,std::vector<float> &par,std::vector<float> &par_err, double c_x, double c_y, double c_wid, double fit_min, double fit_max)
{
    //初期値の設定
    Double_t par_ped[3]={1e4,0,3};
    Double_t par_single[3]={c_y,c_x,c_wid};
    Double_t par_double[2]={14,32};
    //Double_t par_double[2]={4,8};
    par={};
    par_err={};
    //まずはSingle GausianでのPedestalのFitting
    TF1* f_ped = new TF1("gaus","gaus",-20,20);
    //Fittingパラメーターの設定
    f_ped->SetParameters(par_ped);
    //Fitの実行
    hst->Fit(f_ped,"","",-20,10);
    //1 photoのFitting
    TF1* f_ped_single = new TF1("gaus2","gaus(0)+gaus(3)",-10,100);
    //パラメーターの値域設定
    for(int i=0;i<3;i++)
    {
        f_ped_single->FixParameter(i,f_ped->GetParameter(i));
        f_ped_single->SetParameter(i+3,par_single[i]);
    }
    hst->Fit(f_ped_single,"","",fit_min,fit_max);
    //For 2p.e. fitting
    //Define the fitting function
    TF1* f_inter = new TF1("gaus2","gaus(0)+gaus(3)+[6]*TMath::Exp(-0.5*((x-2*[4])/[7])*((x-2*[4])/[7]))",-20,200);
    //パラメーターの値域設定
    for(int i=0;i<8;i++)
    {
        if (i<6)
        {
            f_inter->FixParameter(i,f_ped_single->GetParameter(i));
        }
        else
        {
            f_inter->SetParameter(i,par_double[i-6]);
        }
    }
    //Set parameter range
    Double_t p_lim = f_ped_single->GetParameter(5);
    f_inter->SetParLimits(7,10.0,50.0);
    //Fitting
    hst->Fit(f_inter,"","",fit_min,fit_max*2);

    //For 2p.e. fitting
    //Define the fitting function
    TF1* f_ped_double = new TF1("gaus3","gaus(0)+gaus(3)+[6]*TMath::Exp(-0.5*((x-2*[4])/[7])*((x-2*[4])/[7]))",-10,200);
    //パラメーターの値域設定
    for(int i=0;i<8;i++)
    {
        if (i<3)
        {
            f_ped_double->FixParameter(i,f_inter->GetParameter(i));
        }
        else
        {
            f_ped_double->SetParameter(i,f_inter->GetParameter(i));
        }
    }
    //set parameter range
    Double_t p_lim2 = f_inter->GetParameter(5);
    f_ped_double->SetParLimits(7,10.0,50.0);
    //Fitting
    hst->Fit(f_ped_double,"","",fit_min,fit_max*2);



    //パラメータのReturn
    for(int i=0;i<8;i++)
    {
        par.push_back(f_ped_double->GetParameter(i));
	    par_err.push_back(f_ped_double->GetParError(i));
    }
    
}