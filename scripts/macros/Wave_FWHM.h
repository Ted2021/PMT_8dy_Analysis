#include "vector"
#include "TH1F.h"
#include "TF1.h"
#include "algorithm"
#include "iostream"

int practice()
{
    std::vector<int> arr = {11, 42, -5, 9, -8, 34};

    auto minmax = std::minmax_element(arr.begin(), arr.end());

    std::cout << "min: " << *(minmax.first) << std::endl;
    std::cout << "max: " << *(minmax.second) << std::endl;

    vector<int> num = { 1, 2, 10, 3, 4, 5 };
  
    int maxIndex = distance(num.begin(), max_element(num.begin(), num.end()));
    
    cout << maxIndex;

    return 0;
}

void hoge(TString filesrc,TString filedrk, std::vector<float> &par, TString treename = "Treesource_0", TString treedark = "Treedark_0", int t_min=320, int t_max=360)
{
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

    std::vector<float> temp;
    std::vector<int> result;
    
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
            else if(j>=int_wnd[0])
            {
            //tempベクトルに積分範囲の波高値を入れる
            temp.push_back(wfp[j]-wfn[j]- ave_dark[j] - ped);
            }
        }

        //tempベクトルの最大値を求める(ピーク電圧)
        int maxElement = *std::max_element(temp.begin(), temp.end());
        //std::cout << "max: " << maxElement << std::endl;

        //ピーク電圧のcell番号を取得
        //int maxIndex = distance(temp.begin(), max_element(temp.begin(), temp.end()));
        //std::cout << "max_index: " << maxIndex+int_wnd[0] << std::endl;

        //ピーク半値
        float half = maxElement/2;


        int bin = 0;
        int state = 0;
        for (float num : temp) {
            if (state == 0){
                if (num >= half) {
                    state = 1;
                    result.push_back(bin);
                    bin +=1;
                }
            }
            else if (state ==1)
            {
                if (num < half){
                    state = 2;
                }
                else{
                    result.push_back(bin);
                    bin += 1;
                }
            }
        }

        int firstElement = result.front(); //ピーク半値を超えた位置
        int lastElement = result.back();   //ピーク半値を下がった位置
        //par.push_back(time_e[int_wnd[0]+lastElement] - time_e[int_wnd[0]+firstElement]);
        par.push_back(bin);

        temp.clear();
        result.clear();
        }
    //hst->Draw();
    f->Close();
    f_drk->Close();

}

//Singleの場合
void Make_Charge(TString filesrc,TString filedrk,std::vector<float> &par,TString treename = "Treesource_0", TString treedark = "Treedark_0", int t_min=320, int t_max=360)
{
  //TH1F* hst = new TH1F("hst","hst",330, -20.5,149.5);//Histgramの作成 ヒストグラムの範囲は要検討
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
      par.push_back(chrg);
      //hst->Fill(chrg);
    }
  //hst->Draw();
  f->Close();
  f_drk->Close();
  //return hst;
}

void Peak_wf(TString filesrc,TString filedrk, std::vector<float> &par, TString treename = "Treesource_0", TString treedark = "Treedark_0", int t_min=320, int t_max=360)
{
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

    std::vector<float> temp;
    
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
            else if(j>=int_wnd[0])
            {
            //tempベクトルに積分範囲の波高値を入れる
            temp.push_back(wfp[j]-wfn[j]- ave_dark[j] - ped);
            }
        }

        //tempベクトルの最大値を求める(ピーク電圧)
        //int maxElement = *std::max_element(temp.begin(), temp.end());
        float max = *max_element(temp.begin(), temp.end());
        par.push_back(max);

        temp.clear();
        }
    //hst->Draw();
    f->Close();
    f_drk->Close();

}