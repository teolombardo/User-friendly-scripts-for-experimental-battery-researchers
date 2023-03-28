# Messaging initial
from tkinter import *

root = Tk()

MyLable = Label(root, text="This windows will be closed in 5 seconds, and after that the program will start automatically.")

root.after(5000,lambda:root.destroy())

MyLable.pack()

root.mainloop()

# Author: Teo Lombardo
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import plotly.express as px
import pickle

# Reading the input file
Inputs = pd.read_excel(r'Inputs_ToF.xlsx')

path_data = (Inputs["Parameters' value"][1])
print("Results path =", path_data)

mass_interval_list_name = str(Inputs["Parameters' value"][2])
print("Mass interval list =", mass_interval_list_name)

splitted_peaks = Inputs["Parameters' value"][3].split(",")
for i in range(len(splitted_peaks)):
    important_peaks = [x.replace(" ", "") for x in splitted_peaks]
print("Important ion fragment =", important_peaks)

norm_peak = str(Inputs["Parameters' value"][4])
if norm_peak == "nan":
    norm_peak = ""
print("Normalization peak =", norm_peak)

# Reading the mass interval list
mass_interval_list = pd.read_csv(mass_interval_list_name, sep="	")

path = mass_interval_list_name
with open(path, 'r') as f:
    lines = f.readlines()
    f.close()

# Initializing dictionaries (first run)
Samples_dir = [x[0] for x in os.walk(path_data)]
test_previous = path_data + "\\Box_plots_results"
previous = 0
if test_previous in Samples_dir:
    previous = 1

if previous == 0:
    dict_all = {}
    normalization_dict = {}
    for directory in Samples_dir[1:]:
        if "Box_plots_results" not in directory:
            sample_name = directory.split("\\")[-1]
            print(sample_name)
            dict_all[sample_name] = {}
            normalization_dict[sample_name] = {}
            measurements = os.listdir(directory)
            for i in measurements:
                if ".txt" in i:
                    dict_all[sample_name][i] = {}
                    normalization_dict[sample_name][i] = {}
                    measurement_dict = dict_all[sample_name][i]
                    for line in lines:
                        if line == lines[0]:
                            splitted_line = line.split("	")
                            count_col = 0
                            for j in splitted_line:
                                if j == "Lower Mass (u)":
                                    Lower_mass_n = count_col
                                if j == "Upper Mass (u)":
                                    Upper_mass_n = count_col
                                if j == "Assignment":
                                    Assignment_n = count_col
                                if j == "Description":
                                    Description_n = count_col
                                count_col += 1

                        else:
                            splitted_line = line.split("	")
                            if "total" not in splitted_line:
                                Lower_mass = float(splitted_line[Lower_mass_n])
                                Upper_mass = float(splitted_line[Upper_mass_n])
                                Assignment = splitted_line[Assignment_n]
                                Description = splitted_line[Description_n]
                                measurement_dict[Assignment] = {}
                                measurement_dict[Assignment]["Range masses"] = [Lower_mass, Upper_mass]
                                measurement_dict[Assignment]["m/z"] = []
                                measurement_dict[Assignment]["I"] = []
                                measurement_dict[Assignment]["I_norm"] = []
                                if Assignment == norm_peak:
                                    normalization_dict[sample_name][i][Assignment] = {}
                                    normalization_dict[sample_name][i][Assignment]["Range masses"] = [Lower_mass,
                                                                                                      Upper_mass]
                                    normalization_dict[sample_name][i][Assignment]["m/z"] = []
                                    normalization_dict[sample_name][i][Assignment]["I"] = []

# Filling dictionaries (first run)
if previous == 0:
    for directory in Samples_dir[1:]:
        if "Box_plots_results" not in directory:
            sample_name = directory.split("\\")[-1]
            measurements = os.listdir(directory)
            for m in measurements:
                if ".txt" in m:
                    measurement_dict = dict_all[sample_name][m]
                    Norm_to_fill = normalization_dict[sample_name][m]
                    path_sample = directory + "\\" + m
                    spectrum = pd.read_csv(path_sample, sep="	", skiprows=2)
                    Int_all = 0
                    for i in range(len(spectrum["Intensity"])):
                        Int = float(spectrum["Intensity"][i])
                        Int_all = Int_all + Int
                    for i in range(len(spectrum["m/z"])):
                        mass = float(spectrum["m/z"][i])
                        Int = float(spectrum["Intensity"][i])
                        for j in measurement_dict:
                            peak = measurement_dict[j]
                            if mass >= peak['Range masses'][0] and mass <= peak['Range masses'][1]:
                                peak["m/z"].append(mass)
                                peak["I"].append(Int)
                                peak["I_norm"].append(Int / Int_all)
                                if j == norm_peak:
                                    Norm_to_fill[j]["m/z"].append(mass)
                                    Norm_to_fill[j]["I"].append(Int)

    pkl_filename = "dict_all.pkl"
    with open(pkl_filename, 'wb') as file:
        pickle.dump(dict_all, file)

# Initializing and filling normalization dictionary (not first run)
if previous == 1 and len(norm_peak) > 0:
    print("Restarting")
    pkl_filename = "dict_all.pkl"
    with open(pkl_filename, 'rb') as file:
        dict_all = pickle.load(file)

    normalization_dict = {}
    for directory in Samples_dir[1:]:
        if "Box_plots_results" not in directory:
            sample_name = directory.split("\\")[-1]
            normalization_dict[sample_name] = {}
            measurements = os.listdir(directory)
            for i in measurements:
                if ".txt" in i:
                    normalization_dict[sample_name][i] = {}
                    for line in lines:
                        if line == lines[0]:
                            splitted_line = line.split("	")
                            count_col = 0
                            for j in splitted_line:
                                if j == "Lower Mass (u)":
                                    Lower_mass_n = count_col
                                if j == "Upper Mass (u)":
                                    Upper_mass_n = count_col
                                if j == "Assignment":
                                    Assignment_n = count_col
                                if j == "Description":
                                    Description_n = count_col
                                count_col += 1

                        else:
                            splitted_line = line.split("	")
                            if "total" not in splitted_line:
                                Lower_mass = float(splitted_line[Lower_mass_n])
                                Upper_mass = float(splitted_line[Upper_mass_n])
                                Assignment = splitted_line[Assignment_n]
                                Description = splitted_line[Description_n]
                                if Assignment == norm_peak:
                                    normalization_dict[sample_name][i][Assignment] = {}
                                    normalization_dict[sample_name][i][Assignment]["Range masses"] = [Lower_mass,
                                                                                                      Upper_mass]
                                    normalization_dict[sample_name][i][Assignment]["m/z"] = []
                                    normalization_dict[sample_name][i][Assignment]["I"] = []

    for directory in Samples_dir[1:]:
        if "Box_plots_results" not in directory:
            sample_name = directory.split("\\")[-1]
            measurements = os.listdir(directory)
            for m in measurements:
                if ".txt" in m:
                    Norm_to_fill = normalization_dict[sample_name][m]
                    path_sample = directory + "\\" + m
                    spectrum = pd.read_csv(path_sample, sep="	", skiprows=2)
                    for i in range(len(spectrum["m/z"])):
                        mass = float(spectrum["m/z"][i])
                        Int = float(spectrum["Intensity"][i])
                        if mass >= Norm_to_fill[norm_peak]['Range masses'][0] and mass <= \
                                Norm_to_fill[norm_peak]['Range masses'][1]:
                            Norm_to_fill[norm_peak]["m/z"].append(mass)
                            Norm_to_fill[norm_peak]["I"].append(Int)

# Normalization(s)
dict_area_all = {}
for i in mass_interval_list["Assignment"][1:]:
    dict_area_all[i] = {}
    dict_area_all[i]["sample"] = []
    dict_area_all[i]["Normalized area"] = []
    dict_area_all[i]["Custom normalized area"] = []

for directory in Samples_dir[1:]:
    if "Box_plots_results" not in directory:
        sample_name = directory.split("\\")[-1]
        measurements = os.listdir(directory)
        for m in measurements:
            if ".txt" in m:
                measurement_dict = dict_all[sample_name][m]
                Norm_to_fill = normalization_dict[sample_name][m]
                normalization_area = np.trapz(Norm_to_fill[norm_peak]["I"], Norm_to_fill[norm_peak]["m/z"])
                # print(normalization_area)
                for j in measurement_dict:
                    peak = measurement_dict[j]
                    for i in mass_interval_list["Assignment"][1:]:
                        if j == i:
                            area = np.trapz(peak["I_norm"], peak["m/z"])
                            dict_area_all[j]["Normalized area"].append(area)
                            dict_area_all[j]["sample"].append(sample_name)
                            custom_n_area = np.trapz(peak["I"], peak["m/z"]) / normalization_area
                            dict_area_all[i]["Custom normalized area"].append(custom_n_area)

# Creating results' folder
result_folder_name = "Box_plots_results"
directory_res = path_data + '\\' + result_folder_name
directory_res_tot = directory_res + "\\" + "Normalization by total counts"
directory_res_cust = directory_res + "\\" + "Customized normalization"
directory_res_tot_extra = directory_res_tot + '\\' + "Extra"
directory_res_cust_extra = directory_res_cust + '\\' + "Extra"

try:
    os.makedirs(directory_res)
except:
    pass

try:
    os.makedirs(directory_res_tot)
except:
    pass

if len(norm_peak) > 0:
    try:
        os.makedirs(directory_res_cust)
    except:
        pass

try:
    os.makedirs(directory_res_tot_extra)
except:
    pass

if len(norm_peak) > 0:
    try:
        os.makedirs(directory_res_cust_extra)
    except:
        pass

# Plotting for original and total ions normalization
for peak in dict_area_all:
    df_to_plot = pd.DataFrame.from_dict(dict_area_all[peak])
    fig = px.box(df_to_plot, x="sample", y="Normalized area", points="all", title=peak)
    # fig.update_layout(yaxis_range=[0.55,1.05])

    if peak in important_peaks:
        path_save = directory_res_tot
    else:
        path_save = directory_res_tot_extra

    csv_save = path_save + "\\" + peak + ".csv"
    df_to_plot.to_csv(csv_save, index=False)
    plot_save = path_save + "\\" + peak + ".png"
    fig.write_image(plot_save)

    df_to_plot["Internally normalized area"] = df_to_plot["Normalized area"] / max(df_to_plot["Normalized area"])
    fig = px.box(df_to_plot, x="sample", y="Internally normalized area", points="all", title=peak)
    fig.update_layout(yaxis_range=[0, 1.05])

    if peak in important_peaks:
        path_save = directory_res_tot
    else:
        path_save = directory_res_tot_extra

    csv_save = path_save + "\\Internally_normalized_" + peak + ".csv"
    df_to_plot.to_csv(csv_save, index=False)

    plot_save = path_save + "\\Internally_normalized_" + peak + ".png"
    fig.write_image(plot_save)

# Plotting for customized ions normalization
if len(norm_peak) > 0:
    for peak in dict_area_all:
        df_to_plot = pd.DataFrame.from_dict(dict_area_all[peak])
        fig = px.box(df_to_plot, x="sample", y="Custom normalized area", points="all", title=peak)
        # fig.update_layout(yaxis_range=[0.55,1.05])

        if peak in important_peaks:
            path_save = directory_res_cust
        else:
            path_save = directory_res_cust_extra

        csv_save = path_save + "\\" + peak + ".csv"
        df_to_plot.to_csv(csv_save, index=False)
        plot_save = path_save + "\\" + peak + ".png"
        fig.write_image(plot_save)

# Messaging final

root = Tk()

MyLable = Label(root, text="Box plot analysis done. Enjoy!")
MyLable.pack()
root.mainloop()

##################################################################################### END OF THE PROGRAM #####################################################################################