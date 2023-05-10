import ROOT as RT

#マクロの読み込み
RT.gROOT.LoadMacro("/Users/kiyomoto/reaserch/8dy_Noise/scripts/macros/average_wave.h")
RT.gROOT.LoadMacro("/Users/kiyomoto/reaserch/8dy_Noise/scripts/macros/Noise_counts.h")

#平均波形の計算を行う関数
"""
引数:ROOTファイルのパス, Tree名, branch(channel: 0=>DRS4_CH0&1, 1=>DRS4_CH3&2)
戻り値: 平均時間, 平均波形
"""
def CalcAverage(file_path, tree_name, channel):
    if channel == 0:
        wave_p = "wform1"
        wave_n = "wform0"
    elif channel == 1:
        wave_p = "wform3"
        wave_n = "wform2"
    else:
        print("Error! please input correct Channel num: 0 or 1")
    wave = RT.std.vector(float)()
    time = RT.std.vector(float)()
    RT.Average_Make(file_path, tree_name, wave_p, wave_n, wave, time)
    return time, wave

#ノイズを計算する関数
"""
引数:ROOTファイルのパス, Tree名, branch(channel: 0=>DRS4_CH0&1, 1=>DRS4_CH3&2), 平均波形, 閾値
戻り値:ノイズevent, ノイズcell, ノイズ電荷
"""
def NoiseAnalysis(file_path, tree_name, channel, average_wave, thres = 10):
    if channel == 0:
        wave_p = "wform1"
        wave_n = "wform0"
    elif channel == 1:
        wave_p = "wform3"
        wave_n = "wform2"
    else:
        print("Error! please input correct Channel num: 0 or 1")
        
    wave = RT.std.vector(float)(average_wave)
    event = RT.std.vector(int)()
    cell = RT.std.vector(int)()
    charge = RT.std.vector(float)()
    RT.NoiseAnalysis(file_path, tree_name, wave_p, wave_n, thres, wave, event, cell, charge)
    return event, cell, charge