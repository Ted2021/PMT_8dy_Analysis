import ROOT as RT
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Pmtanalysis:
    def __init__(self, file_path):
        self.file_single = file_path
        self.file_multi = file_path
        self.file_dark = file_path
        
        #Single測定から得られるパラメータ
        self.Gain = "Null"
        self.CE = "Null"
        self.QE = "Null"
        self.f_factor = "Null"

        #Single_ChargeのFittingパラメータ
        self.l_para = []
        self.l_para_err = []
    
    def ShowTrees(self):
        f1 = uproot.open(self.file_single)
        return f1.keys()
    
    def GetAverage(self, tree, branch_p, branch_n):
        #ROOTマクロを読み込む
        RT.gROOT.LoadMacro("/Users/kiyomoto/reaserch/8dy_Noise/scripts/macros/average_wave.h")

        #戻り値を定義する
        time = RT.std.vector(float)()
        wave = RT.std.vector(float)()
        RT.Average_Make(self.file_single, tree, branch_p, branch_n, wave, time)

        return time, wave

    def SingleWaveFit(self):
        para = 0
        error = 1

        return para, error

    def MultiWaveFit(self):
        para = 0
        error = 1

        return para, error
    
    def SingleChargeDist(self, tree, treedark, tmin, tmax):
        #ROOTマクロの読み込み
        RT.gROOT.LoadMacro("/Users/kiyomoto/reaserch/8dy_Noise/scripts/macros/Make_Charge_Dist.h")
        self.hst2d = RT.Make_Charge_Dist(self.file_single, self.file_dark, tree, treedark, tmin, tmax)

        return self.hst2d

    def MultiChargeDist(self, treename, tmin, tmax):
        #ROOTマクロの読み込み
        RT.gROOT.LoadMacro("/Users/kiyomoto/reaserch/8dy_Noise/scripts/macros/Make_Charge_Dist.h")
        hst2d = RT.Make_Charge_Dist_nodark(self.file_single, tmin, tmax, treename)
        return hst2d

    def SingleChargeFit(self, hist):
        #ROOTマクロの読み込み
        RT.gROOT.LoadMacro("/Users/kiyomoto/reaserch/8dy_Noise/scripts/macros/Charge_Dis_Fit_3gau.h")
        single_charge = []
        single_charge_err = []
        single_charge.append(hist.GetMean())
        single_charge_err.append(hist.GetMeanError())

        #1 p.e.のFittingする帯域を設定
        #ROOT Histgram to pandas df
        tmp_hist = [[], []]
        for j in range(hist.GetNbinsX()):
            tmp_hist[0].append(hist.GetXaxis().GetBinCenter(j))
            tmp_hist[1].append(hist.GetBinContent(j))
        tmp_df = pd.DataFrame(np.array(tmp_hist).T, columns=['X', 'Y'])
        #Chargeが10以上の部分での1 p.e.のところのX軸の値を抜粋
        c_max_x = np.mean(tmp_df[tmp_df["Y"] == np.max(tmp_df[tmp_df["X"] > 10]["Y"])][tmp_df["X"]>10]["X"])
        c_max_y = np.mean(tmp_df[tmp_df["Y"] == np.max(tmp_df[tmp_df["X"] > 10]["Y"])][tmp_df["X"]>10]["Y"])

        #Fittingする帯域設定
        fit_min = c_max_x/1.5
        fit_max = c_max_x*1.8

        #変数を設定
        para = RT.std.vector(float)()
        para_err = RT.std.vector(float)()
        #self.l_para = []
        #self.l_para_err = []
        #Fittingの実行
        RT.Fit_Charge_Dist(hist, para, para_err, c_max_x, c_max_y, c_max_x/2, fit_min, fit_max)
        self.l_para.extend(list(para))
        self.l_para_err.extend(list(para_err))
        return tmp_df, self.l_para, self.l_para_err
    
    def SingleGain(self, tmp):
        self.Gain = np.mean(tmp[tmp["Y"] == np.max(tmp[tmp["X"] > 10]["Y"])][tmp["X"]>10]["X"])
        c_max_y = np.mean(tmp[tmp["Y"] == np.max(tmp[tmp["X"] > 10]["Y"])][tmp["X"]>10]["Y"])
        return self.Gain, c_max_y

    def Single_CE(self):
        if self.Gain == "Null":
            print("Error! Please Run 'SingleGain' before This.")
        else:
            N_015 = self.Gain * 0.15
            N_03 = self.Gain * 0.3
            N_15 = self.Gain + 1.5
            self.CE = (N_03 - N_015)/(N_15 - N_015)
            return self.CE, N_015, N_03, N_15
    
    def Ffactor(self):
        self.f_factor = 1 + (self.l_para[5]**2/self.l_para[4]**2)
        return self.f_factor
    
    def SingleChargePlot(self, tmp, para = []):
        #Fitting関数の定義
        def gau(x, p0, p1, p2):
            return p0*np.exp(-0.5*((x-p1)/p2)**2)

        def double_gau(x, l_para):
            p0, p1, p2, p3, p4, p5 = l_para[0], l_para[1], l_para[2], l_para[3], l_para[4], l_para[5]
            return gau(x, p0, p1, p2)+gau(x, p3, p4, p5)

        def triple_gau(x, l_para):
            p0, p1, p2, p3, p4, p5, p6, p7, p8 = l_para[0], l_para[1], l_para[2], l_para[3], l_para[4], l_para[5], l_para[6], l_para[7], l_para[8]
            return gau(x, p0, p1, p2)+gau(x, p3, p4, p5)+gau(x, p6, p7, p8)

        def triple_gau_2(x, l_para):
            p0, p1, p2, p3, p4, p5, p6, p7, p8 = l_para[0], l_para[1], l_para[2], l_para[3], l_para[4], l_para[5], l_para[6], 2*l_para[4], l_para[7]
            return gau(x, p0, p1, p2)+gau(x, p3, p4, p5)+gau(x, p6, p7, p8)
        
        #Fittingする領域とその結果を表示
        fig, ax = plt.subplots()
        ax.plot(tmp["X"], tmp["Y"],c="blue", label="Charge Dist.")
        y_best = triple_gau_2(tmp["X"], self.l_para)
        ax.plot(tmp["X"],y_best,c="red",label="Best Fit")
        y_min, y_max = ax.get_ylim()
        ax.axvline(x = self.Gain, c="green")
        ax.axvline(x = self.Gain * 0.15, c="pink")
        ax.axvline(x = self.Gain * 0.3, c="pink")
        ax.axvline(x = self.Gain * 1.5, c="orange")
        ax.set_yscale("log")
        ax.grid(which="both")
        #ax.set_ylim(0.6, y_max)
        ax.set_ylim(0.6, np.max(y_best))
        plt.xlabel("Charge (mV*ns)")
        plt.ylabel("#Events")
        plt.legend()
        plt.show()


    def MultiChargeFit(self):
        hoge = 1
        return hoge