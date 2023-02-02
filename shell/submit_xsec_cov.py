import subprocess

def main():
    filename = './run_xsec_cov.sh'
    for isyst in range(0,21):
        dial = get_dial_from_id(isyst)
        for ikinematics in [0, 2, 3]:
            list_hist = get_hist_list_from_kine(ikinematics)
            for ihist in list_hist:
                subprocess.call([filename + " " + str(ikinematics) + " " + str(dial) + " " + str(ihist)], shell = True)
            
def get_dial_from_id(i):
    if i == 0 :
        dial = "MaCCQE"
    elif i == 1 :
        dial = "QETwk_HighQ2Weight_1"
    elif i == 2 :
        dial = "QETwk_HighQ2Weight_2"
    elif i == 3 :
        dial = "QETwk_HighQ2Weight_3"
    elif i == 4 :
        dial = "SF_OptPotTwkDial_O16"
    elif i == 5 :
        dial = "MECTwkDial_Norm_O16"
    elif i == 6 :
        dial = "MECTwkDial_PDDWeight_O16_NN"
    elif i == 7 :
        dial = "MECTwkDial_PDDWeight_O16_np"
    elif i == 8 :
        dial = "MECTwkDial_PNNN_Shape"
    elif i == 9 :
        dial = "RES_Eb_O_numu"
    elif i == 10 :
        dial = "BgSclRES"
    elif i == 11 :
        dial = "CA5RES"
    elif i == 12 :
        dial = "MaRES"
    elif i == 13 :
        dial = "PionFSI_AbsProb"
    elif i == 14 :
        dial = "PionFSI_CExHighMomProb"
    elif i == 15 :
        dial = "PionFSI_CExLowMomProb"
    elif i == 16 :
        dial = "PionFSI_InelProb"
    elif i == 17 :
        dial = "PionFSI_QEHighMomProb"
    elif i == 18 :
        dial = "PionFSI_QELowMomProb"
    elif i == 19 :
        dial = "TwkDial_FateNucleonFSI"
    elif i == 20 :
        dial = "CC_DIS_MultiPi_Norm_Nu"

    return dial

def get_hist_list_from_kine(i):
    if i == 0 :
        return ["hist_water_total_multi", "hist_water_proton_multi", "hist_water_pion_multi"]
    elif i == 2 :
        return ["hist_muon_mom", "hist_muon_mom_range", "hist_muon_mom_mcs",
                "hist_proton_mom", "hist_proton_mom_range", "hist_proton_mom_mcs",
                "hist_pion_mom", "hist_pion_mom_range", "hist_pion_mom_mcs"]
    elif i == 3 :
        return ["hist_muon_ang_deg", "hist_proton_ang_deg", "hist_pion_ang_deg"]

if __name__ == "__main__":
    main()
