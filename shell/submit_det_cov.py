import subprocess

def main():
    filename = './run_det_cov.sh'
    #for isyst in range(1,2):
    for isyst in [1,6]:
        dial = get_dial_from_id(isyst)
        for ikinematics in [0, 2, 3]:
            list_hist = get_hist_list_from_kine(ikinematics)
            for ihist in list_hist:
                subprocess.call([filename + " " + str(ikinematics) + " " + str(dial) + " " + str(ihist)], shell = True)

def get_dial_from_id(i):
    if i == 0 :
        dial = "BabyMind"
    elif i == 1 :
        dial = "NinjaBabyMindDistance"
    elif i == 2 :
        dial = "HitThreshold"
    elif i == 3 :
        dial = "MPPCNoise"
    elif i == 4 :
        dial = "Alignment"
    elif i == 5 :
        dial = "AngularResolution"
    elif i == 6 : 
        dial = "DetectionEfficiency"
    elif i == 7 :
        dial = "McsScaling"
    elif i == 8 :
        dial = "VphMean"
    elif i == 9 :
        dial = "VphSigma"
    elif i == 10 :
        dial = "MaterialThickness"
    elif i == 11 :
        dial = "PhysicsList"
    elif i == 12 :
        dial = "Nominal"

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
    elif i == 4 : # CC0pi1p
        return ["hist_dpt", "hist_dalphat", "hist_dphit", "hist_dptx", "hist_dpty"]
    elif i == 5 : # CC0pi2p
        return ["hist_open_cos"]

if __name__ == "__main__":
    main()
