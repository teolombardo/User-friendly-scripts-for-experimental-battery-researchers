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
Inputs = pd.read_excel(r'Inputs_ToF_Images.xlsx')

path_data = (Inputs["Parameters' value"][1])
print("Results path =", path_data)

splitted_peaks = Inputs["Parameters' value"][2].split(",")
for i in range(len(splitted_peaks)):
    important_peaks = [x.replace(" ", "") for x in splitted_peaks]
print("Important ion fragment =", important_peaks)

norm_peak = str(Inputs["Parameters' value"][3])
if norm_peak == "nan":
    norm_peak = ""
print("Normalization peak =", norm_peak)

# Initializing and filling dictionaries (first run)
Samples_dir = [x[0] for x in os.walk(path_data)]
test_previous = path_data + "\\Box_plots_results"
previous = 0
if test_previous in Samples_dir:
    previous = 1

if previous == 0:
    dict_all = {}
    normalization_dict = {}
    dict_sums = {}
    for directory in Samples_dir[1:]:
        if "Box_plots_results" not in directory:
            sample_name = directory.split("\\")[-1]
            # print(sample_name)
            dict_all[sample_name] = {}
            normalization_dict[sample_name] = {}
            measurements = os.listdir(directory)

            for img_n in measurements:
                measure_name = img_n.split(" (")[0]
                peak_name = img_n.split(") ")[1].split("- ")[1].split(".txt")[0]
                # print(measure_name, peak_name)
                dict_all[sample_name][measure_name] = {}
                normalization_dict[sample_name][measure_name] = {}
            for img_n in measurements:
                measure_name = img_n.split(" (")[0]
                peak_name = img_n.split(") ")[1].split("- ")[1].split(".txt")[0]
                # print(measure_name, peak_name)
                dict_all[sample_name][measure_name][peak_name] = {}
                peak_dict = dict_all[sample_name][measure_name][peak_name]
                if peak_name == norm_peak:
                    normalization_dict[sample_name][measure_name][peak_name] = {}
                    Norm_to_fill = normalization_dict[sample_name][measure_name][peak_name]

                path_img = directory + "\\" + img_n
                if peak_name != "total":
                    dict_sums[peak_name] = {}
                    dict_sums[peak_name]["Sample"] = []
                    dict_sums[peak_name]["Intensity"] = []
                    dict_sums[peak_name]["Normalized intensity"] = []
                    dict_sums[peak_name]["Custom normalized intensity"] = []

                with open(path_img, 'r') as f:
                    lines = f.readlines()
                    f.close()

                count_line = 0
                for line in lines:
                    splitted_line = line.split(" ")
                    # print(splitted_line)
                    if splitted_line[0] != "#":
                        header = count_line + 1
                        break
                    count_line += 1

                max_pixel = 0
                max_I = 0
                for line in lines[header:]:
                    splitted_line = line.split(" ")
                    pixel = float(splitted_line[0])
                    if pixel == 0:
                        max_pixel += 1
                    I = float(splitted_line[2])
                    if I > max_I:
                        max_I = I

                peak_dict["I"] = [[] for x in range(max_pixel)]
                peak_dict["In"] = [[] for x in range(max_pixel)]
                if peak_name == norm_peak:
                    Norm_to_fill["I"] = [[] for x in range(max_pixel)]
                for line in lines[header:]:
                    splitted_line = line.split(" ")
                    X = int(splitted_line[0])
                    Y = int(splitted_line[1])
                    I = float(splitted_line[2])
                    In = I / max_I
                    peak_dict["I"][Y].append(I)
                    peak_dict["In"][Y].append(In)
                    if peak_name == norm_peak:
                        Norm_to_fill["I"][Y].append(I)
                peak_dict["I"] = np.array(peak_dict["I"])
                peak_dict["In"] = np.array(peak_dict["In"])
                if peak_name == norm_peak:
                    yi = 0
                    for y in Norm_to_fill["I"]:
                        xi = 0
                        for x in y:
                            if Norm_to_fill["I"][yi][xi] == 0:
                                Norm_to_fill["I"][yi][xi] = 0.5
                            xi += 1
                        yi += 1
                    Norm_to_fill["I"] = np.array(Norm_to_fill["I"])

    pkl_filename = "dict_all_images.pkl"
    with open(pkl_filename, 'wb') as file:
        pickle.dump(dict_all, file)

# Initializing and filling normalization dictionary (not first run)
if previous == 1 and len(norm_peak) > 0:
    print("Restarting")
    pkl_filename = "dict_all_images.pkl"
    with open(pkl_filename, 'rb') as file:
        dict_all = pickle.load(file)

    normalization_dict = {}
    dict_sums = {}
    for directory in Samples_dir[1:]:
        if "Box_plots_results" not in directory:
            sample_name = directory.split("\\")[-1]
            normalization_dict[sample_name] = {}
            measurements = os.listdir(directory)

            for img_n in measurements:
                measure_name = img_n.split(" (")[0]
                peak_name = img_n.split(") ")[1].split("- ")[1].split(".txt")[0]
                normalization_dict[sample_name][measure_name] = {}
            for img_n in measurements:
                measure_name = img_n.split(" (")[0]
                peak_name = img_n.split(") ")[1].split("- ")[1].split(".txt")[0]
                if peak_name == norm_peak:
                    normalization_dict[sample_name][measure_name][peak_name] = {}
                    Norm_to_fill = normalization_dict[sample_name][measure_name][peak_name]

                path_img = directory + "\\" + img_n
                if peak_name != "total":
                    dict_sums[peak_name] = {}
                    dict_sums[peak_name]["Sample"] = []
                    dict_sums[peak_name]["Intensity"] = []
                    dict_sums[peak_name]["Normalized intensity"] = []
                    dict_sums[peak_name]["Custom normalized intensity"] = []

                if peak_name == norm_peak:
                    with open(path_img, 'r') as f:
                        lines = f.readlines()
                        f.close()

                if peak_name == norm_peak:
                    count_line = 0
                    for line in lines:
                        splitted_line = line.split(" ")
                        if splitted_line[0] != "#":
                            header = count_line + 1
                            break
                        count_line += 1

                if peak_name == norm_peak:
                    max_pixel = 0
                    max_I = 0
                    for line in lines[header:]:
                        splitted_line = line.split(" ")
                        pixel = float(splitted_line[0])
                        if pixel == 0:
                            max_pixel += 1
                        I = float(splitted_line[2])
                        if I > max_I:
                            max_I = I

                if peak_name == norm_peak:
                    Norm_to_fill["I"] = [[] for x in range(max_pixel)]
                    for line in lines[header:]:
                        splitted_line = line.split(" ")
                        X = int(splitted_line[0])
                        Y = int(splitted_line[1])
                        I = float(splitted_line[2])
                        In = I / max_I
                        if peak_name == norm_peak:
                            Norm_to_fill["I"][Y].append(I)
                    if peak_name == norm_peak:
                        yi = 0
                        for y in Norm_to_fill["I"]:
                            xi = 0
                            for x in y:
                                if Norm_to_fill["I"][yi][xi] == 0:
                                    Norm_to_fill["I"][yi][xi] = 0.5
                                xi += 1
                            yi += 1
                        Norm_to_fill["I"] = np.array(Norm_to_fill["I"])

# Sums and normalization(s)
for directory in Samples_dir[1:]:
    if "Box_plots_results" not in directory:
        sample_name = directory.split("\\")[-1]
        # dict_all_norm_t[sample_name] = {}
        measurements = os.listdir(directory)
        for img_n in measurements:
            measure_name = img_n.split(" (")[0]
            peak_name = img_n.split(") ")[1].split("- ")[1].split(".txt")[0]
            # for peak in dict_all[sample_name][measure_name]:
            if peak_name != "total":
                dict_sums[peak_name]["Sample"].append(sample_name)
                dict_sums[peak_name]["Intensity"].append(np.sum(dict_all[sample_name][measure_name][peak_name]["I"]))
                dict_all[sample_name][measure_name][peak_name]["I_norm_t"] = \
                dict_all[sample_name][measure_name][peak_name]["I"] / dict_all[sample_name][measure_name]["total"]["I"]
                dict_all[sample_name][measure_name][peak_name]["Cust_norm"] = \
                dict_all[sample_name][measure_name][peak_name]["I"] / \
                normalization_dict[sample_name][measure_name][norm_peak]["I"]

                if peak_name == norm_peak:
                    yi = 0
                    for y in dict_all[sample_name][measure_name][norm_peak]["Cust_norm"]:
                        xi = 0
                        for x in y:
                            if dict_all[sample_name][measure_name][norm_peak]["Cust_norm"][yi][xi] == 0:
                                dict_all[sample_name][measure_name][norm_peak]["Cust_norm"][yi][xi] = 1
                            xi += 1
                        yi += 1

                dict_sums[peak_name]["Normalized intensity"].append(
                    np.sum(dict_all[sample_name][measure_name][peak_name]["I_norm_t"]))
                dict_sums[peak_name]["Custom normalized intensity"].append(
                    np.sum(dict_all[sample_name][measure_name][peak_name]["Cust_norm"]))

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
for peak in dict_sums:
    df_to_plot = pd.DataFrame.from_dict(dict_sums[peak])
    fig = px.box(df_to_plot, x="Sample", y="Normalized intensity", points="all", title=peak)

    if peak in important_peaks:
        path_save = directory_res_tot
    else:
        path_save = directory_res_tot_extra

    csv_save = path_save + "\\" + peak + ".csv"
    df_to_plot.to_csv(csv_save, index=False)

    plot_save = path_save + "\\" + peak + ".png"
    fig.write_image(plot_save)

    df_to_plot["Normalized intensity"] = df_to_plot["Normalized intensity"] / max(df_to_plot["Normalized intensity"])
    fig = px.box(df_to_plot, x="Sample", y="Normalized intensity", points="all", title=peak)
    fig.update_layout(yaxis_range=[0.55, 1.05])

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
    for peak in dict_sums:
        df_to_plot = pd.DataFrame.from_dict(dict_sums[peak])
        fig = px.box(df_to_plot, x="Sample", y="Custom normalized intensity", points="all", title=peak)

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