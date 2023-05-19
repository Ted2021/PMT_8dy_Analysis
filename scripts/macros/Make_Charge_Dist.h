#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "iostream"

//Singleの場合
TH1F* Make_Charge_Dist(TString filesrc,TString filedrk,TString treename = "Treesource_0", TString treedark = "Treedark_0", int t_min=320, int t_max=360)
{
  TH1F* hst = new TH1F("hst","hst",330, -20.5,149.5);//Histgramの作成 ヒストグラムの範囲は要検討
  //ファイルを開く
  TFile* f  = TFile::Open(filesrc);
  TFile* f_drk  = TFile::Open(filedrk);
  //開いたファイルからTree名を取得
  TTree* tr = (TTree*)f->Get(treename);
  TTree* tr_drk = (TTree*)f_drk->Get(treedark);
  //配列を用意
  float time_e[1024] = {0};
  float time_d[1024] = {0};
  float wfp[1024] = {0};//Positiveの波形のデータ格納場所
  float wfp_d[1024] = {0};//Positiveの波形(Dark)のデータ格納場所
  float wfn[1024] = {0};//Negativeの波形のデータ格納場所
  float wfn_d[1024] = {0};//Negativeの波形(Dark)のデータ格納場所
  int sc = 0;
  int sc_d = 0;
  //用意した配列とtree内イベントの対応付け
  tr->SetBranchAddress("time",time_e);
  tr->SetBranchAddress("wform1",wfp);
  tr->SetBranchAddress("wform0",wfn);
  tr->SetBranchAddress("stopcell",&sc);
  tr_drk->SetBranchAddress("time",time_d);
  tr_drk->SetBranchAddress("wform1",wfp_d);
  tr_drk->SetBranchAddress("wform0",wfn_d);
  tr_drk->SetBranchAddress("stopcell",&sc_d);
  
  int int_wnd[2] = {t_min,t_max};//64-74 ns
  int ped_region[2] = {50,150};//pedestalの除去に使用する領域(Dark, Event共通)
  float chrg = 0;

  //Pedestalの平均chargeの初期化
  //1 event毎に計算する
  int nEve = tr->GetEntries(); //tr内の全イベント数を格納
  int nEve_ped = tr_drk->GetEntries(); //tr_drk内の全イベント数を格納
  float ped = 0;
  float ped_drk = 0;
  float ave_dark[1024] = {0};
  
  //Darkの平均波形の作成
  for(int l=0;l<nEve_ped;l++)
    {
      ped_drk = 0;
      tr_drk->GetEntry(l); //l番目のイベントに対応したデータを繰り返す
      //Darkの平均的なPedestalを作成する
      //t_max以下の場合に繰り返す
      for(int m=0;m<int_wnd[1];m++)
        {
          //50以上150未満のときに繰り返す
          if (m>=ped_region[0]&&m<ped_region[1])
          {
            ped_drk += (wfp_d[m]-wfn_d[m])/(ped_region[1]-ped_region[0]);
          }
        }
      //Darkの平均波形を作成する
      for(int m=0;m<int_wnd[1];m++)
        {
          if((sc_d+m)%1024==392)
	        {
	          wfp_d[m] = (wfp_d[m-1] + wfp_d[m+1])*0.5; 
	          wfn_d[m] = (wfn_d[m-1] + wfn_d[m+1])*0.5; 
	        }
          ave_dark[m] += (wfp_d[m]-wfn_d[m] - ped_drk)/nEve_ped;
        }
      //std::cout<<l<<endl;
    }

  //Pedestalの計算

  for(int i=0;i<nEve;i++)//イベント回数で回す
    {
      tr->GetEntry(i);
      //値の初期化
      chrg = 0;  
      ped = 0;
      //Pedestalの作成
      for(int j=0;j<int_wnd[1];j++)
        {
          if (j>=ped_region[0]&&j<ped_region[1])
          {
            ped += (wfp[j]-wfn[j])/(ped_region[1]-ped_region[0]);
          }
        }
      for(int j=0;j<int_wnd[1];j++)//積分時間
      {
        if((sc_d+j)%1024==392)
        {
          wfp_d[j] = (wfp_d[j-1] + wfp_d[j+1])*0.5; 
          wfn_d[j] = (wfn_d[j-1] + wfn_d[j+1])*0.5;
        }
        if((sc+j)%1024==392)
        {
          wfp_d[j] = (wfp_d[j-1] + wfp_d[j+1])*0.5; 
	        wfn_d[j] = (wfn_d[j-1] + wfn_d[j+1])*0.5;
        }
        else if(j>int_wnd[0])
        {
          chrg += ((wfp[j]-wfn[j]- ave_dark[j] - ped)*(time_e[j]-time_e[j-1]));
        }
      }
      hst->Fill(chrg);
    }
  //hst->Draw();
  f->Close();
  f_drk->Close();
  return hst;
}

//Multiの場合
TH1F* Make_Charge_Dist_multi(TString filesrc,TString filedrk, int t_min=320, int t_max=360)
{
  //TH1F* hst = new TH1F("hst","hst", 1000, 100, 2000.5);//Histgramの作成
  TH1F* hst = new TH1F("hst","hst", 1000, 1000, 1500);//Histgramの作成
  TString treename = "Treesource_0";
  TString treedark = "Treedark_0";
  TFile* f  = TFile::Open(filesrc);
  TFile* f_drk  = TFile::Open(filedrk);
  TTree* tr = (TTree*)f->Get(treename);
  TTree* tr_drk = (TTree*)f_drk->Get(treedark);
  float time_e[1024] = {0};
  float time_d[1024] = {0};
  float wfp[1024] = {0};//Positiveの波形のデータ格納場所
  float wfp_d[1024] = {0};//Positiveの波形(Dark)のデータ格納場所
  float wfn[1024] = {0};//Negativeの波形のデータ格納場所
  float wfn_d[1024] = {0};//Negativeの波形(Dark)のデータ格納場所
  int sc = 0;
  int sc_d = 0;
  tr->SetBranchAddress("time",time_e);
  tr->SetBranchAddress("wform3",wfp);
  tr->SetBranchAddress("wform2",wfn);
  tr->SetBranchAddress("stopcell",&sc);
  tr_drk->SetBranchAddress("time",time_d);
  tr_drk->SetBranchAddress("wform3",wfp_d);
  tr_drk->SetBranchAddress("wform2",wfn_d);
  tr_drk->SetBranchAddress("stopcell",&sc_d);
  
  int int_wnd[2] = {t_min,t_max};//64-74 ns
  int ped_region[2] = {50,150};//pedestalの除去に使用する領域(Dark, Event共通)
  float chrg = 0;

  //Pedestalの平均chargeの初期化
  //1 event毎に計算する
  int nEve = tr->GetEntries();
  int nEve_ped = tr_drk->GetEntries();
  float ped = 0;
  float ped_drk = 0;
  float ave_dark[1024] = {0};
  
  //Darkの平均波形の作成
  for(int l=0;l<nEve_ped;l++)
    {
      ped_drk = 0;
      tr_drk->GetEntry(l);
      //Darkの平均的なPedestalを作成する
      for(int m=0;m<int_wnd[1];m++)
        {
          if (m>=ped_region[0]&&m<ped_region[1])
          {
            ped_drk += (wfp_d[m]-wfn_d[m])/(ped_region[1]-ped_region[0]);
          }
        }
      //Darkの平均波形を作成する
      for(int m=0;m<int_wnd[1];m++)
        {
          if((sc_d+m)%1024==392)
	        {
	          wfp_d[m] = (wfp_d[m-1] + wfp_d[m+1])*0.5; 
	          wfn_d[m] = (wfn_d[m-1] + wfn_d[m+1])*0.5; 
	        }
          ave_dark[m] += (wfp_d[m]-wfn_d[m] - ped_drk)/nEve_ped;
        }
      //std::cout<<l<<endl;
    }

  //Pedestalの計算

  for(int i=0;i<nEve;i++)//イベント回数で回す
    {
      tr->GetEntry(i);
      //値の初期化
      chrg = 0;  
      ped = 0;
      //Pedestalの作成
      for(int j=0;j<int_wnd[1];j++)
        {
          if (j>=ped_region[0]&&j<ped_region[1])
          {
            ped += (wfp[j]-wfn[j])/(ped_region[1]-ped_region[0]);
          }
        }
      for(int j=0;j<int_wnd[1];j++)//積分時間
      {
        if((sc_d+j)%1024==392)
        {
          wfp_d[j] = (wfp_d[j-1] + wfp_d[j+1])*0.5; 
          wfn_d[j] = (wfn_d[j-1] + wfn_d[j+1])*0.5;
        }
        if((sc+j)%1024==392)
        {
          wfp_d[j] = (wfp_d[j-1] + wfp_d[j+1])*0.5; 
	        wfn_d[j] = (wfn_d[j-1] + wfn_d[j+1])*0.5;
        }
        else if(j>int_wnd[0])
        {
          chrg += ((wfp[j]-wfn[j]- ave_dark[j] - ped)*(time_e[j]-time_e[j-1]));
        }
      }
      hst->Fill(chrg);
    }
  //hst->Draw();
  f->Close();
  f_drk->Close();
  return hst;
}

TH1F* Make_Charge_Dist_nodark(TString filesrc,int int_start=325,int int_end=375,TString treename="Treetest_0")
{
  std::cout<<filesrc<<std::endl;
  TH1F* hst = new TH1F("hst","hst",1000,0,2000);//Histgramの作成
  TFile* f  = TFile::Open(filesrc);
  TTree* tr = (TTree*)f->Get(treename);
  float time_e[1024] = {0};

  float wfp[1024] = {0};//Positiveの波形のデータ格納場所
  float wfn[1024] = {0};//Negativeの波形のデータ格納場所
  int sc = 0;
  tr->SetBranchAddress("time",time_e);
  tr->SetBranchAddress("wform1",wfp);
  tr->SetBranchAddress("wform0",wfn);
  tr->SetBranchAddress("stopcell",&sc);
  //TH1F* hst = new TH1F("hst","hst",400,0,800);//Histgramの作成
  
  //int int_wnd[2] = {325,375};
  //int ped_region[2] = {50,150};
  int ped_region[2] = {int_start-150,int_start-100};//pedestalの除去に使用する領域(Dark, Event共通)
  if(int_start-150<10)
    {
      ped_region[0]=0;
      ped_region[1]=0;
    }

  float chrg = 0;

  //Pedestalの平均chargeの初期化
  //1 event毎に計算する
  int nEve = tr->GetEntries();
  float ped = 0;
  
  for(int i=0;i<nEve;i++)//イベント回数で回す
    {
      tr->GetEntry(i);
      //値の初期化
      chrg = 0;  
      ped = 0;
      for(int j=0;j<int_end;j++)//積分時間
	{
	  if((sc+j)%1024==392)
	  {
	    wfp[j] = (wfp[j-1] + wfp[j+1])*0.5; 
	    wfn[j] = (wfn[j-1] + wfn[j+1])*0.5; 
	  }
	  
	  if (j>=ped_region[0]&&j<ped_region[1])
	    {
	      ped += (wfp[j]-wfn[j])/(ped_region[1]-ped_region[0]);
	    }
	  else if(j>int_start)
	    {
	      chrg += ((wfp[j]-wfn[j] - ped)*(time_e[j]-time_e[j-1]));
	    }
	}
      hst->Fill(chrg);
    }
  //hst->Draw();
  f->Close();
  return hst;
}
