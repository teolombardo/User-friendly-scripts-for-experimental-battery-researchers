import pandas as pd
import numpy as np
from glob import glob
import os
import matplotlib.pyplot as plt
import math
import statistics
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter, find_peaks
from termcolor import colored
import csv
import pickle

Inputs = pd.read_excel(r'Inputs_CAP.xlsx')

# Software used and file location

software = str(Inputs["Parameters' value"][1].strip().lower())
print("Software =", software)

folder_path = str(Inputs["Parameters' value"][2])
print("Folder path =", folder_path)

filename = str(Inputs["Parameters' value"][3])
print("Filename = ", filename)

# Rate capability test

C_rates_list = []
try:
    splitted_C_rate = Inputs["Parameters' value"][5].split(",")
    for i in range(len(splitted_C_rate)):
        C_rates_list = [x.replace(" ", "") for x in splitted_C_rate]
    print("C_rate_list =", C_rates_list)
except:
    print("No rate capability test detected")

try:
    cycles_per_formations = int(Inputs["Parameters' value"][6])
except:
    cycles_per_formations = Inputs["Parameters' value"][6]
print("Cycle(s) per formation =", cycles_per_formations)

if cycles_per_formations > 0:
    formation = 1
else:
    formation = 0
print("formation = ", formation)

try:
    cycles_per_Crates = int(Inputs["Parameters' value"][7])
except:
    cycles_per_Crates = Inputs["Parameters' value"][7]
print("Cycle(s) per C_rate =", cycles_per_Crates)

# Cycle life test(s)
C_rates_lct = []
try:
    splitted_C_rate_lct = Inputs["Parameters' value"][9].split(",")
    for i in range(len(splitted_C_rate_lct)):
        C_rates_lct = [x.replace(" ", "") for x in splitted_C_rate_lct]
        print("C_rate_list for life cycle test =", C_rates_lct)
except:
    print("No cycle life test detected")
# print(len(C_rates_lct))

Cycles_lct = []
try:
    splitted_Cycles_lct = Inputs["Parameters' value"][10].split(",")
    for i in range(len(splitted_Cycles_lct)):
        Cycles_lct = [int(x.replace(" ", "")) for x in splitted_Cycles_lct]
    print("Cycles per C-rate for life cycle test =", Cycles_lct)
except:
    pass

try:
    splitted_Cycles_lct = int(Inputs["Parameters' value"][10])
    Cycles_lct = []
    Cycles_lct.append(splitted_Cycles_lct)
    print("Cycles per C-rate for life cycle test =", Cycles_lct)
except:
    print("No cycle life detected")
# print(len(Cycles_lct))

pred_cycles = 0

# Smoothing
try:
    points = int(Inputs["Parameters' value"][12])
    assert points > 0
    print("Points per linear interpolation =", points)
except:
    points = 100
    print("Points per linear interpolation = default - ", points)

try:
    sv_windows_lenght = int(Inputs["Parameters' value"][13])
    assert sv_windows_lenght > 0
    if sv_windows_lenght > points:
        sv_windows_lenght = points
    print("Savgol filter - windows length =", sv_windows_lenght)
except:
    sv_windows_lenght = 100
    if sv_windows_lenght > points:
        sv_windows_lenght = points
    print("Savgol filter - windows length = default - ", sv_windows_lenght)

try:
    sv_pol_ord = int(Inputs["Parameters' value"][14])
    assert sv_pol_ord > 0
    print("Savgol filter - polynomial order =", sv_pol_ord)
except:
    sv_pol_ord = 2
    print("Savgol filter - polynomial order = default - ", sv_pol_ord)

# Normalization
mass_active_mg = float(Inputs["Parameters' value"][16])
mass_active_g = mass_active_mg / 1000
if mass_active_g > 0:
    check_mass = 1
else:
    check_mass = 0
print("Normalization mass =", mass_active_mg, "mg")

# ELab
Elab_exp = str(Inputs["Parameters' value"][18])
print("ELab experiment name = ", Elab_exp)
User_name = str(Inputs["Parameters' value"][19])
print("ELab user name = ", User_name)
token_ELab = str(Inputs["Parameters' value"][20])
print("ELab token = ", token_ELab)

# Creating results' folder
test_folders = 0
result_folder_name = "Results_rate_capability_" + filename
directory_res = folder_path + '\\' + result_folder_name
result_folder_CL = "Results_cycle_life_" + filename
directory_res_CL = folder_path + '\\' + result_folder_CL
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    try:
        os.makedirs(directory_res)
    except:
        test_folders += 1
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    try:
        os.makedirs(directory_res_CL)
    except:
        test_folders += 1

if test_folders == 2:
    path_input_log = directory_res + "\\Inputs.log"
    with open(path_input_log, 'r') as f:
        lines = f.readlines()
        f.close()

    for line in lines:
        try:
            splitted_line = line.split(" = ")
        except:
            pass

        try:
            splitted_line2 = splitted_line[1].split("\n")
        except:
            pass

        try:
            if splitted_line[0] == "Software":
                software = splitted_line2[0]
                # print(software)

        except:
            pass

        try:
            if splitted_line[0] == "C_rate_list":
                C_rates_list = []
                rates = splitted_line2[0].split(", ")
                for rate in rates:
                    if rate != '':
                        C_rates_list.append(rate)
                # print(C_rates_list)
        except:
            pass

        try:
            if splitted_line[0] == "Cycle(s) per formation":
                cycles_per_formations = int(splitted_line2[0])
                if cycles_per_formations > 0:
                    formation = 1
                else:
                    formation = 0
                # print(cycles_per_formations)
        except:
            pass

        try:
            if splitted_line[0] == "Cycle(s) per C_rate":
                cycles_per_Crates = int(splitted_line2[0])
                # print(cycles_per_Crates)
        except:
            pass

        try:
            if splitted_line[0] == "C_rate":
                C_rates_lct = []
                rates = splitted_line2[0].split(", ")
                for rate in rates:
                    if rate != '':
                        C_rates_lct.append(rate)
                # print(C_rates_lct)
        except:
            pass

        try:
            if splitted_line[0] == "Cycles per cycle life test":
                Cycles_lct = []
                cycles = splitted_line2[0].split(", ")
                for cycle in cycles:
                    if cycle != '':
                        Cycles_lct.append(int(cycle))
                # print(Cycles_lct)
        except:
            pass

        try:
            if splitted_line[0] == " Cycle(s) to be predicted":
                pred_cycles = int(splitted_line2[0])
                # print(pred_cycles)
        except:
            pass

        try:
            if splitted_line[0] == "Points per linear interpolation":
                points = int(splitted_line2[0])
                # print(points)
        except:
            pass

        try:
            if splitted_line[0] == "Savgol filter - windows length":
                sv_windows_lenght = int(splitted_line2[0])
                # print(sv_windows_lenght)
        except:
            pass

        try:
            if splitted_line[0] == "Savgol filter - polynomial order":
                sv_pol_ord = int(splitted_line2[0])
                # print(sv_pol_ord)
        except:
            pass

        try:
            if splitted_line[0] == "Normalization mass (if any)":
                mass_active_mg = float(splitted_line2[0].split(" (mg)")[0])
                mass_active_g = mass_active_mg / 1000
                # print(mass_active_mg)
        except:
            pass

if test_folders == 2:
    print("Inputs re-read from the log")

# Maccor
if software == "maccor":

    # Initialization
    if mass_active_g > 0:
        list_variables = ["C_rate", "I", "Q", "Qn", "E"]
        x_to_plot = ["Q", "Qn"]
    else:
        list_variables = ["C_rate", "I", "Q", "E"]
        x_to_plot = ["Q"]

    y_to_plot = ["E"]
    if formation == 1:
        number_of_cycles = cycles_per_formations + cycles_per_Crates * (len(C_rates_list) - 1)
    if formation == 0:
        number_of_cycles = cycles_per_Crates * len(C_rates_list)

    if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
        pass
    else:
        number_of_cycles = 0

    path = folder_path + "\\" + filename
    with open(path, 'r') as f:
        lines = f.readlines()
        f.close()

    count_line = 0
    for line in lines:
        splitted_line = line.split("	")
        if splitted_line[0] == "Rec":
            header = count_line
            Maccor_option = 0
            # print(splitted_line)
            for i in range(len(splitted_line)):
                if splitted_line[i] == "Cycle C":
                    N_full_cycle = i
                if splitted_line[i] == "Cap. [Ah]":
                    Q_colums = i
                if splitted_line[i] == "Current [A]":
                    I_cell_A = i
                if splitted_line[i] == "Voltage [V]":
                    Ecell = i
                if splitted_line[i] == "Md":
                    c_d_direction = i

        if splitted_line[0] == "Rec#":
            header = count_line
            Maccor_option = 1
            # print(splitted_line)
            for i in range(len(splitted_line)):
                if splitted_line[i] == "Cyc#":
                    N_full_cycle = i
                if splitted_line[i] == "Capacity (mAh)":
                    Q_colums = i
                if splitted_line[i] == "Current (mA)":
                    I_cell_A = i
                if splitted_line[i] == "Voltage (V)":
                    Ecell = i
                if splitted_line[i] == "State":
                    c_d_direction = i
        count_line += 1

    with open(path, 'r') as f:
        lines = f.readlines()[header + 2:]
        f.close()

    for line in lines:
        splitted_line = line.split("	")
        full_cycle = int(splitted_line[N_full_cycle])
        direction = str(splitted_line[c_d_direction])
        if direction == "D" or direction == "C" or direction == "D\n" or direction == "C\n":
            initial_count_cycle = full_cycle
            break
    max_cycle_test = 0

    for line in lines:
        splitted_line = line.split("	")
        if Maccor_option == 0:
            full_cycle = int(splitted_line[N_full_cycle])
        if Maccor_option == 1:
            full_cycle = int(splitted_line[N_full_cycle]) - initial_count_cycle + 1
        if full_cycle > max_cycle_test:
            max_cycle_test = full_cycle
    # print(max_cycle_test)
    if max_cycle_test < number_of_cycles:
        finish_check = 0
    if max_cycle_test >= number_of_cycles:
        finish_check = 1

    # Rate capability test
    if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:

        list_n_cycles = []

        if finish_check == 0:
            max_cycle = max_cycle_test
            number_of_cycles = max_cycle - 1
            for i in range(number_of_cycles):
                a = i + 1
                list_n_cycles.append(a)

            if formation == 0:
                # print(max_cycle)
                for i in range(len(C_rates_list)):
                    num_rates = len(C_rates_list) - i
                    num_cycles_i = num_rates * cycles_per_Crates
                    if num_cycles_i >= max_cycle and num_cycles_i < max_cycle + cycles_per_Crates:
                        break
                del C_rates_list[num_rates:]
                cycle_per_last = max_cycle - ((len(C_rates_list) - 1) * cycles_per_Crates) - 1
                # print(cycle_per_last)

            if formation == 1:
                if max_cycle > cycles_per_formations:
                    max_cycle = max_cycle - cycles_per_formations
                    # print(max_cycle)
                    for i in range(len(C_rates_list)):
                        num_rates = len(C_rates_list) - i
                        num_cycles_i = num_rates * cycles_per_Crates
                        if num_cycles_i >= max_cycle and num_cycles_i < max_cycle + cycles_per_Crates:
                            break
                    del C_rates_list[num_rates + 1:]
                    cycle_per_last = max_cycle - ((len(C_rates_list) - 2) * cycles_per_Crates) - 1
                    # print(cycle_per_last)

                elif max_cycle <= cycles_per_formations:
                    cycles_per_formations = max_cycle - 1
                    del C_rates_list[1:]

        if finish_check == 1:
            for i in range(number_of_cycles):
                a = i + 1
                list_n_cycles.append(a)

    # Cycle life test
    num_cycles_dict = {}
    if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
        cycles_lct = 0
        cycles_max = number_of_cycles
        cycles_stop = []
        for i in range(len(C_rates_lct)):
            cycles_lct = Cycles_lct[i]
            cycles_max = cycles_max + cycles_lct
            if max_cycle_test < cycles_max:
                del C_rates_lct[i + 1:]
        for i in range(len(C_rates_lct)):
            ycles_lct = Cycles_lct[i]
            if max_cycle_test < cycles_max:
                cycles_stop.append(max_cycle_test - 1)
                if i == 0:
                    last_num_cycles = cycles_stop[i] - number_of_cycles
                else:
                    last_num_cycles = cycles_stop[i] - cycles_stop[i - 1]
                num_cycles_dict[C_rates_lct[i]] = []
                for n in range(last_num_cycles):
                    a = n + 1
                    num_cycles_dict[C_rates_lct[i]].append(a)
            elif max_cycle_test == cycles_max:
                del C_rates_lct[i + 1:]
                cycles_stop.append(max_cycle_test)
                if i == 0:
                    last_num_cycles = cycles_stop[i] - number_of_cycles
                else:
                    last_num_cycles = cycles_stop[i] - cycles_stop[i - 1]
                num_cycles_dict[C_rates_lct[i]] = []
                for n in range(last_num_cycles):
                    a = n + 1
                    num_cycles_dict[C_rates_lct[i]].append(a)
            else:
                cycles_stop.append(cycles_max)
                num_cycles_dict[C_rates_lct[i]] = []
                for n in range(cycles_lct):
                    a = n + 1
                    num_cycles_dict[C_rates_lct[i]].append(a)

    # Initializing the dictionaries
    counters = {}

    # Rate capability test
    if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
        charges = {}
        if formation == 1:
            print("Formation accounted for charge")
            for i in range(cycles_per_formations):
                for j in list_variables:
                    list_temp_c = "list_charge_" + str(C_rates_list[0]) + "_" + str(j) + "_" + str(i + 1)
                    charges[list_temp_c] = []
                    counter_temp = "count_Charge_" + str(C_rates_list[0])
                    counters[counter_temp] = 1
            for a in range(1, len(C_rates_list)):

                if finish_check == 0 and C_rates_list[a] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for i in range(cycles_to_count):
                    for j in list_variables:
                        list_temp_c = "list_charge_" + str(C_rates_list[a]) + "_" + str(j) + "_" + str(i + 1)
                        charges[list_temp_c] = []
                        counter_temp = "count_Charge_" + str(C_rates_list[a])
                        counters[counter_temp] = 1

        if formation == 0:
            print("No formation")
            for a in range(len(C_rates_list)):

                if finish_check == 0 and C_rates_list[a] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for i in range(cycles_to_count):
                    for j in list_variables:
                        list_temp_c = "list_charge_" + str(C_rates_list[a]) + "_" + str(j) + "_" + str(i + 1)
                        charges[list_temp_c] = []
                        counter_temp = "count_Charge_" + str(C_rates_list[a])
                        counters[counter_temp] = 1

        # print(charges)
        # print(len(charges))

        discharges = {}
        if formation == 1:
            print("Formation accounted for discharge")
            for i in range(cycles_per_formations):
                for j in list_variables:
                    list_temp_d = "list_discharge_" + str(C_rates_list[0]) + "_" + str(j) + "_" + str(i + 1)
                    discharges[list_temp_d] = []
                    counter_temp = "count_Discharge_" + str(C_rates_list[0])
                    counters[counter_temp] = 1
            for a in range(1, len(C_rates_list)):

                if finish_check == 0 and C_rates_list[a] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for i in range(cycles_to_count):
                    for j in list_variables:
                        list_temp_d = "list_discharge_" + str(C_rates_list[a]) + "_" + str(j) + "_" + str(i + 1)
                        discharges[list_temp_d] = []
                        counter_temp = "count_Discharge_" + str(C_rates_list[a])
                        counters[counter_temp] = 1

        if formation == 0:
            # print("No formation")
            for a in range(len(C_rates_list)):

                if finish_check == 0 and C_rates_list[a] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for i in range(cycles_to_count):
                    for j in list_variables:
                        list_temp_d = "list_discharge_" + str(C_rates_list[a]) + "_" + str(j) + "_" + str(i + 1)
                        discharges[list_temp_d] = []
                        counter_temp = "count_Discharge_" + str(C_rates_list[a])
                        counters[counter_temp] = 1

    # Cycle life test(s)
    if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
        cycle_life_charges = {}
        cycle_life_discharges = {}
        for i in range(len(C_rates_lct)):
            for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                name_y = str(C_rates_lct[i]) + "_E_" + str(j)
                cycle_life_charges[name_y] = []
                cycle_life_discharges[name_y] = []
                for k in x_to_plot:
                    name_x = str(C_rates_lct[i]) + "_" + k + "_" + str(j)
                    cycle_life_charges[name_x] = []
                    cycle_life_discharges[name_x] = []

    # Filling the dictionaries with (dis)charge results
    path = folder_path + "\\" + filename
    with open(path, 'r') as f:
        lines = f.readlines()
        f.close()

    with open(path, 'r') as f:
        lines = f.readlines()[header + 2:]
        f.close()

    # Rate capability test
    if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
        if formation == 1:
            for line in lines:
                splitted_line = line.split("	")
                if Maccor_option == 0:
                    full_cycle = int(splitted_line[N_full_cycle])
                    Q_Ah = float(splitted_line[Q_colums])
                    I_mA = float(splitted_line[I_cell_A]) / 1000
                if Maccor_option == 1:
                    full_cycle = int(splitted_line[N_full_cycle]) - initial_count_cycle + 1
                    Q_Ah = float(splitted_line[Q_colums]) / 1000
                    I_mA = float(splitted_line[I_cell_A])
                Q_mAh = Q_Ah * 1000
                if mass_active_g > 0:
                    Q_mAh_g = Q_mAh / mass_active_g
                E_V = float(splitted_line[Ecell])
                direction = str(splitted_line[c_d_direction])
                limit_full_cycles_down = 0
                if direction != "R":
                    for i in range(len(C_rates_list)):
                        C_rate = C_rates_list[i]

                        if C_rates_list[i] == C_rates_list[0]:
                            cycles_to_count = cycles_per_formations

                        elif finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                            cycles_to_count = cycle_per_last

                        else:
                            cycles_to_count = cycles_per_Crates

                        if i == 0:
                            limit_full_cycles_up = limit_full_cycles_down + cycles_to_count
                        if i > 0:
                            limit_full_cycles_down = limit_full_cycles_up
                            limit_full_cycles_up = limit_full_cycles_down + cycles_to_count

                        if full_cycle > limit_full_cycles_down and full_cycle <= limit_full_cycles_up:
                            if Q_Ah > 0:
                                if direction == "C" or direction == "C\n":
                                    if counters["count_Charge_" + str(C_rate)] == counters[
                                        "count_Discharge_" + str(C_rate)] + 1:
                                        counters["count_Discharge_" + str(C_rate)] += 1
                                    for a in range(cycles_to_count):
                                        if (a + 1) == counters["count_Charge_" + str(C_rate)]:
                                            charges["list_charge_" + str(C_rates_list[i]) + "_" + "C_rate" + "_" + str(
                                                a + 1)].append(C_rate)
                                            charges["list_charge_" + str(C_rates_list[i]) + "_" + "I" + "_" + str(
                                                a + 1)].append(I_mA)
                                            charges["list_charge_" + str(C_rates_list[i]) + "_" + "Q" + "_" + str(
                                                a + 1)].append(Q_mAh)
                                            if mass_active_g > 0:
                                                charges["list_charge_" + str(C_rates_list[i]) + "_" + "Qn" + "_" + str(
                                                    a + 1)].append(Q_mAh_g)
                                            charges["list_charge_" + str(C_rates_list[i]) + "_" + "E" + "_" + str(
                                                a + 1)].append(E_V)

                                if direction == "D" or direction == "D\n":
                                    if counters["count_Charge_" + str(C_rate)] == counters[
                                        "count_Discharge_" + str(C_rate)]:
                                        counters["count_Charge_" + str(C_rate)] += 1
                                    for a in range(cycles_to_count):
                                        if (a + 1) == counters["count_Discharge_" + str(C_rate)]:
                                            discharges[
                                                "list_discharge_" + str(C_rates_list[i]) + "_" + "C_rate" + "_" + str(
                                                    a + 1)].append(C_rate)
                                            discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "I" + "_" + str(
                                                a + 1)].append(I_mA)
                                            discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "Q" + "_" + str(
                                                a + 1)].append(Q_mAh)
                                            if mass_active_g > 0:
                                                discharges[
                                                    "list_discharge_" + str(C_rates_list[i]) + "_" + "Qn" + "_" + str(
                                                        a + 1)].append(Q_mAh_g)
                                            discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "E" + "_" + str(
                                                a + 1)].append(E_V)
        count_line2 = 0
        count_line = 0
        if formation == 0:
            for line in lines:
                splitted_line = line.split("	")
                if Maccor_option == 0:
                    full_cycle = int(splitted_line[N_full_cycle])
                    Q_Ah = float(splitted_line[Q_colums])
                    I_mA = float(splitted_line[I_cell_A]) / 1000
                if Maccor_option == 1:
                    full_cycle = int(splitted_line[N_full_cycle]) - initial_count_cycle + 1
                    Q_Ah = float(splitted_line[Q_colums]) / 1000
                    I_mA = float(splitted_line[I_cell_A])
                Q_mAh = Q_Ah * 1000
                if mass_active_g > 0:
                    Q_mAh_g = Q_mAh / mass_active_g
                E_V = float(splitted_line[Ecell])
                direction = str(splitted_line[c_d_direction])
                # print(direction)
                limit_full_cycles_down = 0
                if direction != "R":
                    # count_line+=1
                    for i in range(len(C_rates_list)):
                        C_rate = C_rates_list[i]

                        if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                            cycles_to_count = cycle_per_last

                        else:
                            cycles_to_count = cycles_per_Crates

                        if i == 0:
                            limit_full_cycles_up = limit_full_cycles_down + cycles_to_count
                        if i > 0:
                            limit_full_cycles_down = limit_full_cycles_up
                            limit_full_cycles_up = limit_full_cycles_down + cycles_to_count

                        if full_cycle > limit_full_cycles_down and full_cycle <= limit_full_cycles_up:
                            if Q_Ah > 0:
                                if direction == "C" or direction == "C\n":
                                    # print(Q_Ah)
                                    count_line += 1
                                    if counters["count_Charge_" + str(C_rate)] == counters[
                                        "count_Discharge_" + str(C_rate)] + 1:
                                        counters["count_Discharge_" + str(C_rate)] += 1
                                    for a in range(cycles_to_count):
                                        if (a + 1) == counters["count_Charge_" + str(C_rate)]:
                                            count_line2 += 1
                                            charges["list_charge_" + str(C_rates_list[i]) + "_" + "C_rate" + "_" + str(
                                                a + 1)].append(C_rate)
                                            charges["list_charge_" + str(C_rates_list[i]) + "_" + "I" + "_" + str(
                                                a + 1)].append(I_mA)
                                            charges["list_charge_" + str(C_rates_list[i]) + "_" + "Q" + "_" + str(
                                                a + 1)].append(Q_mAh)
                                            if mass_active_g > 0:
                                                charges["list_charge_" + str(C_rates_list[i]) + "_" + "Qn" + "_" + str(
                                                    a + 1)].append(Q_mAh_g)
                                            charges["list_charge_" + str(C_rates_list[i]) + "_" + "E" + "_" + str(
                                                a + 1)].append(E_V)

                                if direction == "D" or direction == "D\n":
                                    count_line += 1
                                    if counters["count_Charge_" + str(C_rate)] == counters[
                                        "count_Discharge_" + str(C_rate)]:
                                        counters["count_Charge_" + str(C_rate)] += 1
                                    for a in range(cycles_to_count):
                                        if (a + 1) == counters["count_Discharge_" + str(C_rate)]:
                                            count_line2 += 1
                                            discharges[
                                                "list_discharge_" + str(C_rates_list[i]) + "_" + "C_rate" + "_" + str(
                                                    a + 1)].append(C_rate)
                                            discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "I" + "_" + str(
                                                a + 1)].append(I_mA)
                                            discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "Q" + "_" + str(
                                                a + 1)].append(Q_mAh)
                                            if mass_active_g > 0:
                                                discharges[
                                                    "list_discharge_" + str(C_rates_list[i]) + "_" + "Qn" + "_" + str(
                                                        a + 1)].append(Q_mAh_g)
                                            discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "E" + "_" + str(
                                                a + 1)].append(E_V)

    # Cycle life test(s)
    if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
        for line in lines:
            splitted_line = line.split("	")
            if Maccor_option == 0:
                full_cycle = int(splitted_line[N_full_cycle])
                Q_Ah = float(splitted_line[Q_colums])
                I_mA = float(splitted_line[I_cell_A]) / 1000
            if Maccor_option == 1:
                full_cycle = int(splitted_line[N_full_cycle]) - initial_count_cycle + 1
                Q_Ah = float(splitted_line[Q_colums]) / 1000
                I_mA = float(splitted_line[I_cell_A])
            Q_mAh = Q_Ah * 1000
            if mass_active_g > 0:
                Q_mAh_g = Q_mAh / mass_active_g
            E_V = float(splitted_line[Ecell])
            direction = str(splitted_line[c_d_direction])
            if direction != "R":
                for i in range(len(C_rates_lct)):
                    if full_cycle > number_of_cycles and full_cycle <= cycles_stop[i]:
                        num_cycle_life = full_cycle - number_of_cycles - 1
                        if direction == "C" or direction == "C\n":
                            cycle_life_charges[str(C_rates_lct[i]) + "_E_" + str(num_cycle_life)].append(E_V)
                            cycle_life_charges[str(C_rates_lct[i]) + "_Q_" + str(num_cycle_life)].append(Q_mAh)
                            if mass_active_g > 0:
                                cycle_life_charges[str(C_rates_lct[i]) + "_Qn_" + str(num_cycle_life)].append(Q_mAh_g)
                        if direction == "D" or direction == "D\n":
                            cycle_life_discharges[str(C_rates_lct[i]) + "_E_" + str(num_cycle_life)].append(E_V)
                            cycle_life_discharges[str(C_rates_lct[i]) + "_Q_" + str(num_cycle_life)].append(Q_mAh)
                            if mass_active_g > 0:
                                cycle_life_discharges[str(C_rates_lct[i]) + "_Qn_" + str(num_cycle_life)].append(
                                    Q_mAh_g)

if software == "ec-lab":

    # Initialization
    list_variables = ["C_rate", "I", "Q", "Qn", "E"]
    x_to_plot = ["Q", "Qn"]
    y_to_plot = ["E"]
    if formation == 1:
        number_of_cycles = cycles_per_formations + cycles_per_Crates * (len(C_rates_list) - 1)
    if formation == 0:
        number_of_cycles = cycles_per_Crates * len(C_rates_list)

    if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
        pass
    else:
        number_of_cycles = 0

    # Going from "," to "."
    path = folder_path + "\\" + filename

    reading_file = open(path, "r")

    new_file_content = ""
    for raw in reading_file:
        stripped_raw = raw.strip()
        new_line = stripped_raw.replace(",", ".")
        new_file_content += new_line + "\n"
    reading_file.close()

    writing_file = open(path, "w")
    writing_file.write(new_file_content)
    writing_file.close()

    # Reading and recovering general features of the electrode
    with open(path, 'r') as f:
        lines = f.readlines()
        f.close()

    for line in range(1, 2):
        splitted_line = lines[line].split(' ')
        header = int(splitted_line[4])
        # print("header =", header)

    if check_mass == 0:
        for line in range(22, 23):
            splitted_line = lines[line].split(' ')
            mass_active_mg = float(splitted_line[5])
            mass_active_g = float(splitted_line[5]) / 1000
            print("mass active material = ", mass_active_g, " g")

    for line in range(30, 31):
        splitted_line = lines[line].split(' ')
        electrode_surface_cm2 = float(splitted_line[4])
        print("electrode surface = ", electrode_surface_cm2, "cm2")

    AM_loading = mass_active_mg / electrode_surface_cm2
    print("AM loading = ", AM_loading, " g cm-2")

    for line in range(header - 1, header):
        splitted_line = lines[line].split('	')
        # print(splitted_line)
        counter_temp = 0
        for i in splitted_line:
            if i == "mode":
                mode = int(counter_temp)
                # print("mode", counter_temp)
            if i == "Ns":
                Ns = int(counter_temp)
                # print("Ns", counter_temp)
            if i == "cycle number":
                N_full_cycle = int(counter_temp)
            if i == "Ecell/V":
                Ecell = int(counter_temp)
                # print("Ecell/V", counter_temp)
            if i == "<I>/mA":
                I_cell = int(counter_temp)
                # print("<I>/mA", counter_temp)
            if i == "Q discharge/mA.h":
                Q_dis_mAh = int(counter_temp)
                # print("Q discharge/mA.h", counter_temp)
            if i == "Q charge/mA.h":
                Q_c_mAh = int(counter_temp)
                # print("Q charge/mA.h", counter_temp)
            counter_temp += 1

    with open(path, 'r') as f:
        lines = f.readlines()[header:]
        f.close()

    test_1 = 0
    count_cycle = 0
    for line in lines:
        splitted_line = line.split("	")
        half_1 = int(splitted_line[int(Ns)])

        if test_1 == 1:
            if half_1 != half_0:
                count_cycle += 1

        half_0 = int(splitted_line[int(Ns)])
        test_1 = 1
    max_cycle_test = math.floor(count_cycle / 2)

    if max_cycle_test < number_of_cycles:
        finish_check = 0
    if max_cycle_test >= number_of_cycles:
        finish_check = 1

    # Rate capacility test
    if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
        list_n_cycles = []

        if finish_check == 0:
            max_cycle = max_cycle_test
            number_of_cycles = max_cycle - 1
            for i in range(number_of_cycles):
                a = i + 1
                list_n_cycles.append(a)

            if formation == 0:
                # print(max_cycle)
                for i in range(len(C_rates_list)):
                    num_rates = len(C_rates_list) - i
                    num_cycles_i = num_rates * cycles_per_Crates
                    if num_cycles_i >= max_cycle and num_cycles_i < max_cycle + cycles_per_Crates:
                        break
                del C_rates_list[num_rates:]
                # print(C_rates_list)
                cycle_per_last = max_cycle - ((len(C_rates_list) - 1) * cycles_per_Crates) - 1
                # print(cycle_per_last)

            if formation == 1:
                if max_cycle > cycles_per_formations:
                    max_cycle = max_cycle - cycles_per_formations
                    # print(max_cycle)
                    for i in range(len(C_rates_list)):
                        num_rates = len(C_rates_list) - i
                        num_cycles_i = num_rates * cycles_per_Crates
                        if num_cycles_i >= max_cycle and num_cycles_i < max_cycle + cycles_per_Crates:
                            break
                    del C_rates_list[num_rates + 1:]
                    cycle_per_last = max_cycle - ((len(C_rates_list) - 2) * cycles_per_Crates) - 1
                    # print(cycle_per_last)

                elif max_cycle <= cycles_per_formations:
                    cycles_per_formations = max_cycle - 1
                    del C_rates_list[1:]

        if finish_check == 1:
            for i in range(number_of_cycles):
                a = i + 1
                list_n_cycles.append(a)

    # Cycle life test
    num_cycles_dict = {}
    if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
        cycles_lct = 0
        cycles_max = number_of_cycles
        cycles_stop = []
        for i in range(len(C_rates_lct)):
            cycles_lct = Cycles_lct[i]
            cycles_max = cycles_max + cycles_lct
            if max_cycle_test < cycles_max:
                del C_rates_lct[i + 1:]
                cycles_stop.append(max_cycle_test - 1)
                if i == 0:
                    last_num_cycles = cycles_stop[i] - number_of_cycles
                else:
                    last_num_cycles = cycles_stop[i] - cycles_stop[i - 1]
                num_cycles_dict[C_rates_lct[i]] = []
                for n in range(last_num_cycles):
                    a = n + 1
                    num_cycles_dict[C_rates_lct[i]].append(a)
            elif max_cycle_test == cycles_max:
                del C_rates_lct[i + 1:]
                cycles_stop.append(max_cycle_test)
                if i == 0:
                    last_num_cycles = cycles_stop[i] - number_of_cycles
                else:
                    last_num_cycles = cycles_stop[i] - cycles_stop[i - 1]
                num_cycles_dict[C_rates_lct[i]] = []
                for n in range(last_num_cycles):
                    a = n + 1
                    num_cycles_dict[C_rates_lct[i]].append(a)
            else:
                cycles_stop.append(cycles_max)
                num_cycles_dict[C_rates_lct[i]] = []
                for n in range(cycles_lct):
                    a = n + 1
                    num_cycles_dict[C_rates_lct[i]].append(a)

                    # Initializing the dictionaries
    counters = {}

    # Rate capability test
    if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
        charges = {}
        if formation == 1:
            print("Formation accounted for charge")
            for i in range(cycles_per_formations):
                for j in list_variables:
                    list_temp_c = "list_charge_" + str(C_rates_list[0]) + "_" + str(j) + "_" + str(i + 1)
                    charges[list_temp_c] = []
                    counter_temp = "count_Charge_" + str(C_rates_list[0])
                    counters[counter_temp] = 1
            for a in range(1, len(C_rates_list)):

                if finish_check == 0 and C_rates_list[a] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last
                else:
                    cycles_to_count = cycles_per_Crates

                for i in range(cycles_to_count):
                    for j in list_variables:
                        list_temp_c = "list_charge_" + str(C_rates_list[a]) + "_" + str(j) + "_" + str(i + 1)
                        charges[list_temp_c] = []
                        counter_temp = "count_Charge_" + str(C_rates_list[a])
                        counters[counter_temp] = 1

        if formation == 0:
            print("No formation")
            for a in range(len(C_rates_list)):

                if finish_check == 0 and C_rates_list[a] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last
                else:
                    cycles_to_count = cycles_per_Crates

                for i in range(cycles_to_count):
                    for j in list_variables:
                        list_temp_c = "list_charge_" + str(C_rates_list[a]) + "_" + str(j) + "_" + str(i + 1)
                        charges[list_temp_c] = []
                        counter_temp = "count_Charge_" + str(C_rates_list[a])
                        counters[counter_temp] = 1

        # print(charges)
        # print(len(charges))

        discharges = {}
        if formation == 1:
            print("Formation accounted for discharge")
            for i in range(cycles_per_formations):
                for j in list_variables:
                    list_temp_d = "list_discharge_" + str(C_rates_list[0]) + "_" + str(j) + "_" + str(i + 1)
                    discharges[list_temp_d] = []
                    counter_temp = "count_Discharge_" + str(C_rates_list[0])
                    counters[counter_temp] = 1

            for a in range(1, len(C_rates_list)):

                if finish_check == 0 and C_rates_list[a] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last
                else:
                    cycles_to_count = cycles_per_Crates

                for i in range(cycles_to_count):
                    for j in list_variables:
                        list_temp_d = "list_discharge_" + str(C_rates_list[a]) + "_" + str(j) + "_" + str(i + 1)
                        discharges[list_temp_d] = []
                        counter_temp = "count_Discharge_" + str(C_rates_list[a])
                        counters[counter_temp] = 1

        if formation == 0:
            # print("No formation")
            for a in range(len(C_rates_list)):

                if finish_check == 0 and C_rates_list[a] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last
                else:
                    cycles_to_count = cycles_per_Crates

                for i in range(cycles_to_count):
                    for j in list_variables:
                        list_temp_d = "list_discharge_" + str(C_rates_list[a]) + "_" + str(j) + "_" + str(i + 1)
                        discharges[list_temp_d] = []
                        counter_temp = "count_Discharge_" + str(C_rates_list[a])
                        counters[counter_temp] = 1

    # Cycle life test(s)
    if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
        cycle_life_charges = {}
        cycle_life_discharges = {}
        for i in range(len(C_rates_lct)):
            for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                name_y = str(C_rates_lct[i]) + "_E_" + str(j)
                cycle_life_charges[name_y] = []
                cycle_life_discharges[name_y] = []
                for k in x_to_plot:
                    name_x = str(C_rates_lct[i]) + "_" + k + "_" + str(j)
                    cycle_life_charges[name_x] = []
                    cycle_life_discharges[name_x] = []

            counter_temp = "count_Discharge_CL_" + str(C_rates_lct[i])
            counters[counter_temp] = 1  # number_of_cycles + 1
            counter_temp = "count_Charge_CL_" + str(C_rates_lct[i])
            counters[counter_temp] = 1

    # Reading the data and filling the lists
    # Rate capability test
    if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
        if formation == 1:
            half_cycles_formation = cycles_per_formations * 2
            for line in lines:
                stripped_line = line.split('	')
                subprogram = int(stripped_line[int(mode)])
                half_cycle = int(stripped_line[int(Ns)])
                E_V = float(stripped_line[int(Ecell)])
                I_mA = float(stripped_line[int(I_cell)])
                Q_discharge_mAh = float(stripped_line[int(Q_dis_mAh)])
                # print(Q_discharge_mAh)
                Specific_Q_discharge_mAh_g = Q_discharge_mAh / mass_active_g
                Q_charge_mAh = float(stripped_line[int(Q_c_mAh)])
                Specific_Q_charge_mAh_g = Q_charge_mAh / mass_active_g
                if subprogram == 1:
                    if half_cycle <= half_cycles_formation:
                        C_rate = C_rates_list[0]
                        if (half_cycle % 2) != 0:
                            if counters["count_Charge_" + str(C_rate)] == counters[
                                "count_Discharge_" + str(C_rate)] + 1:
                                counters["count_Charge_" + str(C_rate)] += 1
                            for a in range(cycles_per_formations):
                                if (a + 1) == counters["count_Charge_" + str(C_rate)]:
                                    charges["list_charge_" + str(C_rates_list[0]) + "_" + "C_rate" + "_" + str(
                                        a + 1)].append(C_rate)
                                    charges[
                                        "list_charge_" + str(C_rates_list[0]) + "_" + "I" + "_" + str(a + 1)].append(
                                        I_mA)
                                    charges[
                                        "list_charge_" + str(C_rates_list[0]) + "_" + "Q" + "_" + str(a + 1)].append(
                                        Q_charge_mAh)
                                    charges[
                                        "list_charge_" + str(C_rates_list[0]) + "_" + "Qn" + "_" + str(a + 1)].append(
                                        Specific_Q_charge_mAh_g)
                                    charges[
                                        "list_charge_" + str(C_rates_list[0]) + "_" + "E" + "_" + str(a + 1)].append(
                                        E_V)

                        if (half_cycle % 2) == 0:
                            if counters["count_Charge_" + str(C_rate)] == counters["count_Discharge_" + str(C_rate)]:
                                counters["count_Charge_" + str(C_rate)] += 1
                            for a in range(cycles_per_formations):
                                if (a + 1) == counters["count_Discharge_" + str(C_rate)]:
                                    discharges["list_discharge_" + str(C_rates_list[0]) + "_" + "C_rate" + "_" + str(
                                        a + 1)].append(C_rate)
                                    discharges[
                                        "list_discharge_" + str(C_rates_list[0]) + "_" + "I" + "_" + str(a + 1)].append(
                                        I_mA)
                                    discharges[
                                        "list_discharge_" + str(C_rates_list[0]) + "_" + "Q" + "_" + str(a + 1)].append(
                                        Q_discharge_mAh)
                                    discharges["list_discharge_" + str(C_rates_list[0]) + "_" + "Qn" + "_" + str(
                                        a + 1)].append(Specific_Q_discharge_mAh_g)
                                    discharges[
                                        "list_discharge_" + str(C_rates_list[0]) + "_" + "E" + "_" + str(a + 1)].append(
                                        E_V)

                    for i in range(1, len(C_rates_list)):
                        C_rate = C_rates_list[i]

                        if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                            cycles_to_count = cycle_per_last
                        else:
                            cycles_to_count = cycles_per_Crates

                        limit_half_cycles_up = half_cycles_formation + i * 2
                        limit_half_cycles_down = limit_half_cycles_up - 2
                        if half_cycle > limit_half_cycles_down and half_cycle <= limit_half_cycles_up:
                            if (half_cycle % 2) != 0:
                                if counters["count_Charge_" + str(C_rate)] == counters[
                                    "count_Discharge_" + str(C_rate)] + 1:
                                    counters["count_Discharge_" + str(C_rate)] += 1
                                for a in range(cycles_to_count):
                                    if (a + 1) == counters["count_Charge_" + str(C_rate)]:
                                        charges["list_charge_" + str(C_rates_list[i]) + "_" + "C_rate" + "_" + str(
                                            a + 1)].append(C_rate)
                                        charges["list_charge_" + str(C_rates_list[i]) + "_" + "I" + "_" + str(
                                            a + 1)].append(I_mA)
                                        charges["list_charge_" + str(C_rates_list[i]) + "_" + "Q" + "_" + str(
                                            a + 1)].append(Q_charge_mAh)
                                        charges["list_charge_" + str(C_rates_list[i]) + "_" + "Qn" + "_" + str(
                                            a + 1)].append(Specific_Q_charge_mAh_g)
                                        charges["list_charge_" + str(C_rates_list[i]) + "_" + "E" + "_" + str(
                                            a + 1)].append(E_V)

                            if (half_cycle % 2) == 0:
                                if counters["count_Charge_" + str(C_rate)] == counters[
                                    "count_Discharge_" + str(C_rate)]:
                                    counters["count_Charge_" + str(C_rate)] += 1
                                for a in range(cycles_to_count):
                                    if (a + 1) == counters["count_Discharge_" + str(C_rate)]:
                                        discharges[
                                            "list_discharge_" + str(C_rates_list[i]) + "_" + "C_rate" + "_" + str(
                                                a + 1)].append(C_rate)
                                        discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "I" + "_" + str(
                                            a + 1)].append(I_mA)
                                        discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "Q" + "_" + str(
                                            a + 1)].append(Q_discharge_mAh)
                                        discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "Qn" + "_" + str(
                                            a + 1)].append(Specific_Q_discharge_mAh_g)
                                        discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "E" + "_" + str(
                                            a + 1)].append(E_V)

        if formation == 0:
            half_cycles_formation = cycles_per_formations * 2
            for line in lines:
                stripped_line = line.split('	')
                subprogram = int(stripped_line[int(mode)])
                half_cycle = int(stripped_line[int(Ns)])
                E_V = float(stripped_line[int(Ecell)])
                I_mA = float(stripped_line[int(I_cell)])
                Q_discharge_mAh = float(stripped_line[int(Q_dis_mAh)])
                # print(Q_discharge_mAh)
                Specific_Q_discharge_mAh_g = Q_discharge_mAh / mass_active_g
                Q_charge_mAh = float(stripped_line[int(Q_c_mAh)])
                Specific_Q_charge_mAh_g = Q_charge_mAh / mass_active_g
                if subprogram == 1:
                    for i in range(len(C_rates_list)):
                        C_rate = C_rates_list[i]

                        if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                            cycles_to_count = cycle_per_last
                        else:
                            cycles_to_count = cycles_per_Crates

                        limit_half_cycles_up = (i + 1) * 2
                        limit_half_cycles_down = limit_half_cycles_up - 2
                        if half_cycle > limit_half_cycles_down and half_cycle <= limit_half_cycles_up:
                            if (half_cycle % 2) != 0:
                                if counters["count_Charge_" + str(C_rate)] == counters[
                                    "count_Discharge_" + str(C_rate)] + 1:
                                    counters["count_Discharge_" + str(C_rate)] += 1
                                for a in range(cycles_to_count):
                                    if (a + 1) == counters["count_Charge_" + str(C_rate)]:
                                        charges["list_charge_" + str(C_rates_list[i]) + "_" + "C_rate" + "_" + str(
                                            a + 1)].append(C_rate)
                                        charges["list_charge_" + str(C_rates_list[i]) + "_" + "I" + "_" + str(
                                            a + 1)].append(I_mA)
                                        charges["list_charge_" + str(C_rates_list[i]) + "_" + "Q" + "_" + str(
                                            a + 1)].append(Q_charge_mAh)
                                        charges["list_charge_" + str(C_rates_list[i]) + "_" + "Qn" + "_" + str(
                                            a + 1)].append(Specific_Q_charge_mAh_g)
                                        charges["list_charge_" + str(C_rates_list[i]) + "_" + "E" + "_" + str(
                                            a + 1)].append(E_V)

                            if (half_cycle % 2) == 0:
                                if counters["count_Charge_" + str(C_rate)] == counters[
                                    "count_Discharge_" + str(C_rate)]:
                                    counters["count_Charge_" + str(C_rate)] += 1
                                for a in range(cycles_to_count):
                                    if (a + 1) == counters["count_Discharge_" + str(C_rate)]:
                                        discharges[
                                            "list_discharge_" + str(C_rates_list[i]) + "_" + "C_rate" + "_" + str(
                                                a + 1)].append(C_rate)
                                        discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "I" + "_" + str(
                                            a + 1)].append(I_mA)
                                        discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "Q" + "_" + str(
                                            a + 1)].append(Q_discharge_mAh)
                                        discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "Qn" + "_" + str(
                                            a + 1)].append(Specific_Q_discharge_mAh_g)
                                        discharges["list_discharge_" + str(C_rates_list[i]) + "_" + "E" + "_" + str(
                                            a + 1)].append(E_V)

    # Cycle life test(s)
    if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
        cycle = 1  # number_of_cycles + 1
        for line in lines:
            stripped_line = line.split('	')
            subprogram = int(stripped_line[int(mode)])
            half_cycle = int(stripped_line[int(Ns)])
            # full_cycle = float(splitted_line[N_full_cycle])+1
            E_V = float(stripped_line[int(Ecell)])
            I_mA = float(stripped_line[int(I_cell)])
            Q_discharge_mAh = float(stripped_line[int(Q_dis_mAh)])
            # print(Q_discharge_mAh)
            Specific_Q_discharge_mAh_g = Q_discharge_mAh / mass_active_g
            Q_charge_mAh = float(stripped_line[int(Q_c_mAh)])
            Specific_Q_charge_mAh_g = Q_charge_mAh / mass_active_g
            if subprogram == 1:
                for i in range(len(C_rates_lct)):
                    if (half_cycle % 2) != 0:
                        if counters["count_Charge_CL_" + str(C_rates_lct[i])] == counters[
                            "count_Discharge_CL_" + str(C_rates_lct[i])] + 1:
                            counters["count_Discharge_CL_" + str(C_rates_lct[i])] += 1
                            cycle = counters["count_Charge_CL_" + str(C_rates_lct[i])]
                        if cycle > number_of_cycles and cycle <= cycles_stop[i]:
                            cycle_CL = cycle - number_of_cycles
                            # num_cycle_life = math.ceil(full_cycle)-number_of_cycles-1
                            cycle_life_charges[str(C_rates_lct[i]) + "_E_" + str(cycle_CL - 1)].append(E_V)
                            cycle_life_charges[str(C_rates_lct[i]) + "_Q_" + str(cycle_CL - 1)].append(Q_charge_mAh)
                            cycle_life_charges[str(C_rates_lct[i]) + "_Qn_" + str(cycle_CL - 1)].append(
                                Specific_Q_charge_mAh_g)
                    if (half_cycle % 2) == 0:
                        if counters["count_Charge_CL_" + str(C_rates_lct[i])] == counters[
                            "count_Discharge_CL_" + str(C_rates_lct[i])]:
                            counters["count_Charge_CL_" + str(C_rates_lct[i])] += 1
                            cycle = counters["count_Discharge_CL_" + str(C_rates_lct[i])]
                        if cycle > number_of_cycles and cycle <= cycles_stop[i]:
                            cycle_CL = cycle - number_of_cycles
                            cycle_life_discharges[str(C_rates_lct[i]) + "_E_" + str(cycle_CL - 1)].append(E_V)
                            cycle_life_discharges[str(C_rates_lct[i]) + "_Q_" + str(cycle_CL - 1)].append(
                                Q_discharge_mAh)
                            cycle_life_discharges[str(C_rates_lct[i]) + "_Qn_" + str(cycle_CL - 1)].append(
                                Specific_Q_discharge_mAh_g)

# Common part
# Linear interpolation and averaging

# Rate capability test
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    limits = {}
    averages = {}
    std = {}
    averages_plus_std = {}
    averages_minus_std = {}
    maxs = {}
    mins = {}
    interpolated_discharges = {}
    interpolated_charges = {}
    for d in range(2):
        if d == 0:
            dictionary = discharges
            name_dict = "list_discharge_"
            direction = "discharge"
            interp = interpolated_discharges
        if d == 1:
            dictionary = charges
            name_dict = "list_charge_"
            direction = "charge"
            interp = interpolated_charges
        if formation == 1:
            for k in list_variables:
                list_min = []
                list_max = []
                for j in range(cycles_per_formations):
                    list_min.append(min(dictionary[name_dict + str(C_rates_list[0]) + "_" + str(k) + "_" + str(j + 1)]))
                    list_max.append(max(dictionary[name_dict + str(C_rates_list[0]) + "_" + str(k) + "_" + str(j + 1)]))
                    min_temp = max(list_min)
                    max_temp = min(list_max)
                    temp_n_min = str(C_rates_list[0]) + "_" + str(k) + "_min_" + direction
                    temp_n_max = str(C_rates_list[0]) + "_" + str(k) + "_max_" + direction
                    limits[temp_n_min] = min_temp
                    limits[temp_n_max] = max_temp

            for l in y_to_plot:
                for k in x_to_plot:
                    list_interpolated = []
                    Min = limits[str(C_rates_list[0]) + "_" + str(l) + "_min_" + direction]
                    Max = limits[str(C_rates_list[0]) + "_" + str(l) + "_max_" + direction]
                    # I = np.linspace(Min, Max, points)
                    if d == 0:
                        I = np.linspace(Max, Min, points)
                    if d == 1:
                        I = np.linspace(Min, Max, points)
                    # print(I)
                    for j in range(cycles_per_formations):
                        x_curve = dictionary[name_dict + str(C_rates_list[0]) + "_" + str(k) + "_" + str(j + 1)]
                        y_curve = dictionary[name_dict + str(C_rates_list[0]) + "_" + str(l) + "_" + str(j + 1)]
                        interpole = interp1d(y_curve, x_curve)
                        interpole_I = interpole(I)
                        list_interpolated.append(interpole_I)
                        if k == "Q":
                            name_int_E = str(C_rates_list[0]) + "_E_" + str(j + 1)
                            name_int_Q = str(C_rates_list[0]) + "_Q_" + str(j + 1)
                            interp[name_int_E] = I
                            interp[name_int_Q] = interpole_I
                    arrays = [np.array(x) for x in list_interpolated]
                    ave_x = [np.mean(k) for k in zip(*arrays)]
                    std_x = [np.std(k) for k in zip(*arrays)]
                    ave_plus_sd = [i + j for i, j in zip(ave_x, std_x)]
                    ave_minus_sd = [i - j for i, j in zip(ave_x, std_x)]
                    name_x = str(C_rates_list[0]) + "_" + str(k) + "_" + str(l) + "_x_" + direction
                    name_y = str(C_rates_list[0]) + "_" + str(k) + "_" + str(l) + "_y_" + direction
                    averages[name_x] = ave_x
                    averages[name_y] = I.tolist()
                    std[name_x] = std_x
                    averages_plus_std[name_x] = ave_plus_sd
                    averages_minus_std[name_x] = ave_minus_sd
                    max_x = max(averages[name_x])
                    min_x = min(averages[name_x])
                    max_y = max(averages[name_y])
                    min_y = min(averages[name_y])
                    maxs[name_x] = max_x
                    mins[name_x] = min_x
                    maxs[name_y] = max_y
                    mins[name_y] = min_y

            for i in range(1, len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for k in list_variables:
                    list_min = []
                    list_max = []
                    for j in range(cycles_to_count):
                        list_min.append(
                            min(dictionary[name_dict + str(C_rates_list[i]) + "_" + str(k) + "_" + str(j + 1)]))
                        list_max.append(
                            max(dictionary[name_dict + str(C_rates_list[i]) + "_" + str(k) + "_" + str(j + 1)]))
                    min_temp = max(list_min)
                    max_temp = min(list_max)
                    temp_n_min = str(C_rates_list[i]) + "_" + str(k) + "_min_" + direction
                    temp_n_max = str(C_rates_list[i]) + "_" + str(k) + "_max_" + direction
                    limits[temp_n_min] = min_temp
                    limits[temp_n_max] = max_temp

            for i in range(1, len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for l in y_to_plot:
                    for k in x_to_plot:
                        list_interpolated = []
                        Min = limits[str(C_rates_list[i]) + "_" + str(l) + "_min_" + direction]
                        Max = limits[str(C_rates_list[i]) + "_" + str(l) + "_max_" + direction]
                        if d == 0:
                            I = np.linspace(Max, Min, points)
                        if d == 1:
                            I = np.linspace(Min, Max, points)
                        # print(I)
                        for j in range(cycles_to_count):
                            x_curve = dictionary[name_dict + str(C_rates_list[i]) + "_" + str(k) + "_" + str(j + 1)]
                            y_curve = dictionary[name_dict + str(C_rates_list[i]) + "_" + str(l) + "_" + str(j + 1)]
                            # Not optimal procedure
                            if len(x_curve) > 3:
                                interpole = interp1d(y_curve, x_curve)
                                interpole_I = interpole(I)
                            else:
                                interpole_I = [float(max(x_curve))] * points
                                interpole_I = np.array(interpole_I)
                            list_interpolated.append(interpole_I)
                            if k == "Q":
                                name_int_E = str(C_rates_list[i]) + "_E_" + str(j + 1)
                                name_int_Q = str(C_rates_list[i]) + "_Q_" + str(j + 1)
                                interp[name_int_E] = I
                                interp[name_int_Q] = interpole_I
                        arrays = [np.array(x) for x in list_interpolated]
                        ave_x = [np.mean(k) for k in zip(*arrays)]
                        std_x = [np.std(k) for k in zip(*arrays)]
                        ave_plus_sd = [i + j for i, j in zip(ave_x, std_x)]
                        ave_minus_sd = [i - j for i, j in zip(ave_x, std_x)]
                        name_x = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_x_" + direction
                        name_y = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_y_" + direction
                        averages[name_x] = ave_x
                        averages[name_y] = I.tolist()
                        std[name_x] = std_x
                        averages_plus_std[name_x] = ave_plus_sd
                        averages_minus_std[name_x] = ave_minus_sd
                        max_x = max(averages[name_x])
                        min_x = min(averages[name_x])
                        max_y = max(averages[name_y])
                        min_y = min(averages[name_y])
                        maxs[name_x] = max_x
                        mins[name_x] = min_x
                        maxs[name_y] = max_y
                        mins[name_y] = min_y

        if formation == 0:
            for i in range(len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for k in list_variables:
                    list_min = []
                    list_max = []
                    for j in range(cycles_to_count):
                        list_min.append(
                            min(dictionary[name_dict + str(C_rates_list[i]) + "_" + str(k) + "_" + str(j + 1)]))
                        list_max.append(
                            max(dictionary[name_dict + str(C_rates_list[i]) + "_" + str(k) + "_" + str(j + 1)]))
                    min_temp = max(list_min)
                    max_temp = min(list_max)
                    temp_n_min = str(C_rates_list[i]) + "_" + str(k) + "_min_" + direction
                    temp_n_max = str(C_rates_list[i]) + "_" + str(k) + "_max_" + direction
                    limits[temp_n_min] = min_temp
                    limits[temp_n_max] = max_temp

            for i in range(len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for l in y_to_plot:
                    for k in x_to_plot:
                        list_interpolated = []
                        Min = limits[str(C_rates_list[i]) + "_" + str(l) + "_min_" + direction]
                        Max = limits[str(C_rates_list[i]) + "_" + str(l) + "_max_" + direction]
                        if d == 0:
                            I = np.linspace(Max, Min, points)
                        if d == 1:
                            I = np.linspace(Min, Max, points)
                        # print(I)
                        for j in range(cycles_to_count):
                            x_curve = dictionary[name_dict + str(C_rates_list[i]) + "_" + str(k) + "_" + str(j + 1)]
                            y_curve = dictionary[name_dict + str(C_rates_list[i]) + "_" + str(l) + "_" + str(j + 1)]
                            if len(x_curve) > 3:
                                interpole = interp1d(y_curve, x_curve)
                                interpole_I = interpole(I)
                            else:
                                interpole_I = [float(max(x_curve))] * points
                                interpole_I = np.array(interpole_I)
                            list_interpolated.append(interpole_I)
                            if k == "Q":
                                name_int_E = str(C_rates_list[i]) + "_E_" + str(j + 1)
                                name_int_Q = str(C_rates_list[i]) + "_Q_" + str(j + 1)
                                interp[name_int_E] = I
                                interp[name_int_Q] = interpole_I
                        arrays = [np.array(x) for x in list_interpolated]
                        ave_x = [np.mean(k) for k in zip(*arrays)]
                        std_x = [np.std(k) for k in zip(*arrays)]
                        ave_plus_sd = [i + j for i, j in zip(ave_x, std_x)]
                        ave_minus_sd = [i - j for i, j in zip(ave_x, std_x)]
                        name_x = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_x_" + direction
                        name_y = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_y_" + direction
                        averages[name_x] = ave_x
                        averages[name_y] = I.tolist()
                        std[name_x] = std_x
                        averages_plus_std[name_x] = ave_plus_sd
                        averages_minus_std[name_x] = ave_minus_sd
                        max_x = max(averages[name_x])
                        min_x = min(averages[name_x])
                        max_y = max(averages[name_y])
                        min_y = min(averages[name_y])
                        maxs[name_x] = max_x
                        mins[name_x] = min_x
                        maxs[name_y] = max_y
                        mins[name_y] = min_y

                    # Cycle life test(s)
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    interpolated_CL_discharges = {}
    interpolated_CL_charges = {}
    for d in range(2):
        if d == 0:
            dictionary = cycle_life_discharges
            interp = interpolated_CL_discharges
        if d == 1:
            dictionary = cycle_life_charges
            interp = interpolated_CL_charges
        for i in range(len(C_rates_lct)):
            for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                for k in x_to_plot:
                    name_x = str(C_rates_lct[i]) + "_" + k + "_" + str(j)
                    name_y = str(C_rates_lct[i]) + "_E_" + str(j)
                    x = dictionary[name_x]
                    y = dictionary[name_y]
                    Min = min(y)
                    Max = max(y)
                    if d == 0:
                        I = np.linspace(Max, Min, points)
                    if d == 1:
                        I = np.linspace(Min, Max, points)
                    interpole = interp1d(y, x)
                    interpole_I = interpole(I)
                    interp[name_x] = interpole_I
                    interp[name_y] = I


def round_down(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier) / multiplier


# Doing the derivatives for dQdV
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    dQdV_dict_d = {}
    dQdV_dict_c = {}
    # size_leveling = 15
    # filt = np.ones(size_leveling)/size_leveling
    # filt_minus = -round_down(size_leveling/2)
    # filt_minus = int(filt_minus)
    # filt_plus = round_down(size_leveling/2)
    # filt_plus = int(filt_plus)
    for d in range(2):
        y_smooth = []
        if d == 0:
            dictionary = interpolated_discharges
        if d == 1:
            dictionary = interpolated_charges
        if formation == 1:
            for j in range(cycles_per_formations):
                x_curve = dictionary[str(C_rates_list[0]) + "_Q_" + str(j + 1)]
                y_curve = dictionary[str(C_rates_list[0]) + "_E_" + str(j + 1)]
                y_smooth = savgol_filter(y_curve, sv_windows_lenght, sv_pol_ord)
                dQdV = np.gradient(x_curve, y_smooth)
                name_dict_y = str(C_rates_list[0]) + "_dQdV_" + str(j + 1)
                name_dict_x = str(C_rates_list[0]) + "_E_" + str(j + 1)
                if d == 0:
                    dQdV_dict_d[name_dict_x] = y_smooth
                    dQdV_dict_d[name_dict_y] = dQdV.tolist()
                if d == 1:
                    dQdV_dict_c[name_dict_x] = y_smooth
                    dQdV_dict_c[name_dict_y] = dQdV.tolist()

            for i in range(1, len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for j in range(cycles_to_count):
                    x_curve = dictionary[str(C_rates_list[i]) + "_Q_" + str(j + 1)]
                    y_curve = dictionary[str(C_rates_list[i]) + "_E_" + str(j + 1)]
                    y_smooth = savgol_filter(y_curve, sv_windows_lenght, sv_pol_ord)
                    dQdV = np.gradient(x_curve, y_smooth)
                    name_dict_y = str(C_rates_list[i]) + "_dQdV_" + str(j + 1)
                    name_dict_x = str(C_rates_list[i]) + "_E_" + str(j + 1)
                    if d == 0:
                        dQdV_dict_d[name_dict_x] = y_smooth
                        dQdV_dict_d[name_dict_y] = dQdV.tolist()
                    if d == 1:
                        dQdV_dict_c[name_dict_x] = y_smooth
                        dQdV_dict_c[name_dict_y] = dQdV.tolist()

        if formation == 0:
            for i in range(len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for j in range(cycles_to_count):
                    try:
                        x_curve = dictionary[str(C_rates_list[i]) + "_Q_" + str(j + 1)]
                        y_curve = dictionary[str(C_rates_list[i]) + "_E_" + str(j + 1)]
                        y_smooth = savgol_filter(y_curve, sv_windows_lenght, sv_pol_ord)
                        dQdV = np.gradient(x_curve, y_smooth)
                        name_dict_y = str(C_rates_list[i]) + "_dQdV_" + str(j + 1)
                        name_dict_x = str(C_rates_list[i]) + "_E_" + str(j + 1)
                        if d == 0:
                            dQdV_dict_d[name_dict_x] = y_smooth
                            dQdV_dict_d[name_dict_y] = dQdV.tolist()
                        if d == 1:
                            dQdV_dict_c[name_dict_x] = y_smooth
                            dQdV_dict_c[name_dict_y] = dQdV.tolist()
                    except:
                        continue

# Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    dQdv_CL_d = {}
    dQdv_CL_c = {}
    for d in range(2):
        if d == 0:
            dictionary = interpolated_CL_discharges
            dQdv_dict = dQdv_CL_d
        if d == 1:
            dictionary = interpolated_CL_charges
            dQdv_dict = dQdv_CL_c
        for i in range(len(C_rates_lct)):
            for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                try:
                    name_x = str(C_rates_lct[i]) + "_Q_" + str(j)
                    name_y = str(C_rates_lct[i]) + "_E_" + str(j)
                    x = dictionary[name_x]
                    y = dictionary[name_y]
                    y_smooth = savgol_filter(y, sv_windows_lenght, sv_pol_ord)
                    dQdV = np.gradient(x, y_smooth)
                    dQdv_dict[str(C_rates_lct[i]) + "_dQdV_" + str(j)] = dQdV.tolist()
                    dQdv_dict[name_y] = y_smooth
                except:
                    continue

# Doing the derivatives for dVdQ
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    dVdQ_dict_d = {}
    dVdQ_dict_c = {}
    # size_leveling = 15
    # filt = np.ones(size_leveling)/size_leveling
    # filt_minus = -round_down(size_leveling/2)
    # filt_minus = int(filt_minus)
    # filt_plus = round_down(size_leveling/2)
    # filt_plus = int(filt_plus)
    for d in range(2):
        y_smooth = []
        if d == 0:
            dictionary = interpolated_discharges
        if d == 1:
            dictionary = interpolated_charges
        if formation == 1:
            for j in range(cycles_per_formations):
                x_curve = dictionary[str(C_rates_list[0]) + "_Q_" + str(j + 1)]
                y_curve = dictionary[str(C_rates_list[0]) + "_E_" + str(j + 1)]
                y_smooth = savgol_filter(y_curve, sv_windows_lenght, sv_pol_ord)
                dVdQ = np.gradient(y_smooth, x_curve)
                name_dict_y = str(C_rates_list[0]) + "_dVdQ_" + str(j + 1)
                name_dict_x = str(C_rates_list[0]) + "_Q_" + str(j + 1)
                if d == 0:
                    dVdQ_dict_d[name_dict_x] = 100 * x_curve / max(x_curve)
                    dVdQ_dict_d[name_dict_y] = dVdQ.tolist()
                if d == 1:
                    dVdQ_dict_c[name_dict_x] = 100 * x_curve / max(x_curve)
                    dVdQ_dict_c[name_dict_y] = dVdQ.tolist()

            for i in range(1, len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for j in range(cycles_to_count):
                    x_curve = dictionary[str(C_rates_list[i]) + "_Q_" + str(j + 1)]
                    y_curve = dictionary[str(C_rates_list[i]) + "_E_" + str(j + 1)]
                    y_smooth = savgol_filter(y_curve, sv_windows_lenght, sv_pol_ord)
                    dVdQ = np.gradient(y_smooth, x_curve)
                    name_dict_y = str(C_rates_list[i]) + "_dVdQ_" + str(j + 1)
                    name_dict_x = str(C_rates_list[i]) + "_Q_" + str(j + 1)
                    if d == 0:
                        dVdQ_dict_d[name_dict_x] = 100 * x_curve / max(x_curve)
                        dVdQ_dict_d[name_dict_y] = dVdQ.tolist()
                    if d == 1:
                        dVdQ_dict_c[name_dict_x] = 100 * x_curve / max(x_curve)
                        dVdQ_dict_c[name_dict_y] = dVdQ.tolist()

        if formation == 0:
            for i in range(len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last

                else:
                    cycles_to_count = cycles_per_Crates

                for j in range(cycles_to_count):
                    try:
                        x_curve = dictionary[str(C_rates_list[i]) + "_Q_" + str(j + 1)]
                        y_curve = dictionary[str(C_rates_list[i]) + "_E_" + str(j + 1)]
                        y_smooth = savgol_filter(y_curve, sv_windows_lenght, sv_pol_ord)
                        dVdQ = np.gradient(y_smooth, x_curve)
                        name_dict_y = str(C_rates_list[i]) + "_dVdQ_" + str(j + 1)
                        name_dict_x = str(C_rates_list[i]) + "_Q_" + str(j + 1)
                        if d == 0:
                            dVdQ_dict_d[name_dict_x] = 100 * x_curve / max(x_curve)
                            dVdQ_dict_d[name_dict_y] = dVdQ.tolist()
                        if d == 1:
                            dVdQ_dict_c[name_dict_x] = 100 * x_curve / max(x_curve)
                            dVdQ_dict_c[name_dict_y] = dVdQ.tolist()
                    except:
                        continue

# Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    dVdQ_CL_d = {}
    dVdQ_CL_c = {}
    for d in range(2):
        if d == 0:
            dictionary = interpolated_CL_discharges
            dVdQ_dict = dVdQ_CL_d
        if d == 1:
            dictionary = interpolated_CL_charges
            dVdQ_dict = dVdQ_CL_c
        for i in range(len(C_rates_lct)):
            for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                try:
                    name_x = str(C_rates_lct[i]) + "_Q_" + str(j)
                    name_y = str(C_rates_lct[i]) + "_E_" + str(j)
                    x = dictionary[name_x]
                    y = dictionary[name_y]
                    y_smooth = savgol_filter(y, sv_windows_lenght, sv_pol_ord)
                    # x_smooth = savgol_filter(y, sv_windows_lenght, sv_pol_ord)
                    dVdQ = np.gradient(y_smooth, x)
                    dVdQ_dict[str(C_rates_lct[i]) + "_dVdQ_" + str(j)] = dVdQ.tolist()
                    dVdQ_dict[name_x] = 100 * x / max(x)
                except:
                    continue
# Report averaged (dis)charge curves in DataFame for easier printing
# Averages
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    overall_dic = {}
    overall_dic["charge"] = {}
    overall_dic["discharge"] = {}
    for key in averages:
        splitted_key = key.split("_")
        splitted_key2 = key.split("_Q")
        splitted_key3 = key.split("_E")
        key_dir = splitted_key[-1]
        if splitted_key[-2] == "x":
            new_key_ave = splitted_key2[0] + "_" + splitted_key3[0].split("_")[-1] + "_ave"
            new_key_std = splitted_key2[0] + "_" + splitted_key3[0].split("_")[-1] + "_std"
            overall_dic[key_dir][new_key_ave] = averages[key]
            overall_dic[key_dir][new_key_std] = std[key]
        if splitted_key[-2] == "y":
            new_key_ave = splitted_key2[0] + "_" + splitted_key[-3] + "_ave"
            overall_dic[key_dir][new_key_ave] = averages[key]

    ave_results_charge = pd.DataFrame(overall_dic["charge"])
    ave_results_discharge = pd.DataFrame(overall_dic["discharge"])

# Recovering end of (dis)charge information
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    if formation == 1:
        ends_d = {}
        ends_c = {}
        for d in range(2):
            for l in y_to_plot:
                for k in x_to_plot:
                    end_of_discharge = []
                    end_of_charge = []
                    for j in range(cycles_per_formations):
                        if d == 0:
                            x_curve = discharges[
                                "list_discharge_" + str(C_rates_list[0]) + "_" + str(k) + "_" + str(j + 1)]
                            end_of_discharge.append(max(x_curve))
                            ends_d[k] = end_of_discharge
                        if d == 1:
                            x_curve = charges["list_charge_" + str(C_rates_list[0]) + "_" + str(k) + "_" + str(j + 1)]
                            end_of_charge.append(max(x_curve))
                            ends_c[k] = end_of_charge
                    for i in range(1, len(C_rates_list)):

                        if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                            cycles_to_count = cycle_per_last

                        else:
                            cycles_to_count = cycles_per_Crates

                        for j in range(cycles_to_count):
                            if d == 0:
                                x_curve = discharges[
                                    "list_discharge_" + str(C_rates_list[i]) + "_" + str(k) + "_" + str(j + 1)]
                                end_of_discharge.append(max(x_curve))
                                ends_d[k] = end_of_discharge
                            if d == 1:
                                x_curve = charges[
                                    "list_charge_" + str(C_rates_list[i]) + "_" + str(k) + "_" + str(j + 1)]
                                end_of_charge.append(max(x_curve))
                                ends_c[k] = end_of_charge

    if formation == 0:
        ends_d = {}
        ends_c = {}
        for d in range(2):
            for l in y_to_plot:
                for k in x_to_plot:
                    end_of_discharge = []
                    end_of_charge = []
                    for i in range(len(C_rates_list)):

                        if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                            cycles_to_count = cycle_per_last

                        else:
                            cycles_to_count = cycles_per_Crates

                        for j in range(cycles_to_count):
                            if d == 0:
                                x_curve = discharges[
                                    "list_discharge_" + str(C_rates_list[i]) + "_" + str(k) + "_" + str(j + 1)]
                                if x_curve[-1] < max(x_curve) - 0.02 * max(x_curve):
                                    Q_to_append = x_curve[-1] + max(x_curve)
                                else:
                                    Q_to_append = x_curve[-1]
                                end_of_discharge.append(Q_to_append)
                                ends_d[k] = end_of_discharge
                            if d == 1:
                                x_curve = charges[
                                    "list_charge_" + str(C_rates_list[i]) + "_" + str(k) + "_" + str(j + 1)]
                                if x_curve[-1] < max(x_curve) - 0.02 * max(x_curve):
                                    Q_to_append = x_curve[-1] + max(x_curve)
                                else:
                                    Q_to_append = x_curve[-1]
                                end_of_charge.append(Q_to_append)
                                ends_c[k] = end_of_charge

                            # Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    ends_CL_d = {}
    ends_CL_c = {}
    CV_Qd = []
    CV_Qc = []
    for d in range(2):
        if d == 0:
            dictionary = cycle_life_discharges
            end_dict = ends_CL_d
            CV_Q = CV_Qd
        if d == 1:
            dictionary = cycle_life_charges
            end_dict = ends_CL_c
            CV_Q = CV_Qc
        for i in range(len(C_rates_lct)):
            for k in x_to_plot:
                name_end = str(C_rates_lct[i]) + "_" + str(k)
                end_dict[name_end] = []
                for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                    name_x = str(C_rates_lct[i]) + "_" + k + "_" + str(j)
                    x = dictionary[name_x]
                    if x[-1] < max(x) - 0.02 * max(x):
                        Q_to_append = x[-1] + max(x)
                        if k == "Q":
                            CV_Q_i = (x[-1] / max(x)) * 100
                            CV_Q.append(CV_Q_i)
                    else:
                        Q_to_append = x[-1]

                    end_dict[name_end].append(Q_to_append)

# Getting average potential of each cycle
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    if formation == 1:
        ave_V_discharge = []
        std_V_discharge = []
        ave_V_charge = []
        std_V_charge = []
        for d in range(2):
            for l in y_to_plot:
                for j in range(cycles_per_formations):
                    if d == 0:
                        x_curve = discharges["list_discharge_" + str(C_rates_list[0]) + "_" + str(l) + "_" + str(j + 1)]
                        ave_V_discharge.append(np.average(x_curve))
                        std_V_discharge.append(np.std(x_curve))
                    if d == 1:
                        x_curve = charges["list_charge_" + str(C_rates_list[0]) + "_" + str(l) + "_" + str(j + 1)]
                        ave_V_charge.append(np.average(x_curve))
                        std_V_charge.append(np.std(x_curve))

                for i in range(1, len(C_rates_list)):
                    if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                        cycles_to_count = cycle_per_last

                    else:
                        cycles_to_count = cycles_per_Crates

                    for j in range(cycles_to_count):
                        if d == 0:
                            x_curve = discharges[
                                "list_discharge_" + str(C_rates_list[i]) + "_" + str(l) + "_" + str(j + 1)]
                            ave_V_discharge.append(np.average(x_curve))
                            std_V_discharge.append(np.std(x_curve))
                        if d == 1:
                            x_curve = charges["list_charge_" + str(C_rates_list[i]) + "_" + str(l) + "_" + str(j + 1)]
                            ave_V_charge.append(np.average(x_curve))
                            std_V_charge.append(np.std(x_curve))

    if formation == 0:
        ave_V_discharge = []
        std_V_discharge = []
        ave_V_charge = []
        std_V_charge = []
        for d in range(2):
            for l in y_to_plot:
                for i in range(len(C_rates_list)):

                    if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                        cycles_to_count = cycle_per_last

                    else:
                        cycles_to_count = cycles_per_Crates

                    for j in range(cycles_to_count):
                        if d == 0:
                            x_curve = discharges[
                                "list_discharge_" + str(C_rates_list[i]) + "_" + str(l) + "_" + str(j + 1)]
                            ave_V_discharge.append(np.average(x_curve))
                            std_V_discharge.append(np.std(x_curve))
                        if d == 1:
                            x_curve = charges["list_charge_" + str(C_rates_list[i]) + "_" + str(l) + "_" + str(j + 1)]
                            ave_V_charge.append(np.average(x_curve))
                            std_V_charge.append(np.std(x_curve))

                        # Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    ave_V_CL_d = []
    std_V_CL_d = []
    ave_V_CL_c = []
    std_V_CL_c = []
    max_E = 0
    min_E = 10000
    for d in range(2):
        if d == 0:
            dictionary = cycle_life_discharges
            ave_list = ave_V_CL_d
            std_list = std_V_CL_d
        if d == 1:
            dictionary = cycle_life_charges
            ave_list = ave_V_CL_c
            std_list = std_V_CL_c
        for i in range(len(C_rates_lct)):
            for l in y_to_plot:
                # name_end = str(C_rates_lct[i]) + "_" + str(k)
                # end_dict[name_end] = []
                for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                    name_x = str(C_rates_lct[i]) + "_" + l + "_" + str(j)
                    x = dictionary[name_x]
                    min_E_i = min(x)
                    if min_E_i < min_E:
                        min_E = min_E_i
                    max_E_i = max(x)
                    if max_E_i > max_E:
                        max_E = max_E_i
                    ave_list.append(np.average(x))
                    std_list.append(np.std(x))

# Common interpolation
# Cycle life test(s)
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    Mins_CL = {}
    Maxs_CL = {}
    for d in range(2):
        if d == 0:
            dictionary = cycle_life_discharges
            name = "charge_"
        if d == 1:
            dictionary = cycle_life_charges
            name = "discharge_"
        for i in range(len(C_rates_lct)):
            name = name + str(i)
            Mins_CL[name] = []
            Maxs_CL[name] = []
            for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                name_x = str(C_rates_lct[i]) + "_Q_" + str(j)
                name_y = str(C_rates_lct[i]) + "_E_" + str(j)
                x = dictionary[name_x]
                y = dictionary[name_y]
                Min_i = min(y)
                Mins_CL[name].append(Min_i)
                Max_i = max(y)
                Maxs_CL[name].append(Max_i)

    points_common = points * 2
    interpolated_CL_discharges_common = {}
    interpolated_CL_charges_common = {}
    for d in range(2):
        if d == 0:
            dictionary = cycle_life_discharges
            interp = interpolated_CL_discharges_common
            name = "charge_"
        if d == 1:
            dictionary = cycle_life_charges
            interp = interpolated_CL_charges_common
            name = "discharge_"
        for i in range(len(C_rates_lct)):
            name = name + str(i)
            Max = min(Maxs_CL[name])
            Min = max(Mins_CL[name])
            for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                name_x = str(C_rates_lct[i]) + "_Q_" + str(j)
                name_y = str(C_rates_lct[i]) + "_E_" + str(j)
                x = dictionary[name_x]
                y = dictionary[name_y]
                if d == 0:
                    I = np.linspace(Max, Min, points_common)
                if d == 1:
                    I = np.linspace(Min, Max, points_common)
                interpole = interp1d(y, x)
                interpole_I = interpole(I)
                interp[name_x] = interpole_I
                interp[name_y] = I

# Getting variance of the capacity differences
# Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    var_Q_discharge_CL1 = []
    var_Q_charge_CL1 = []
    var_Q_discharge_CL2 = []
    var_Q_charge_CL2 = []
    for d in range(2):
        if d == 0:
            dictionary = interpolated_CL_discharges_common
            var_list1 = var_Q_discharge_CL1
            var_list2 = var_Q_discharge_CL2
        if d == 1:
            dictionary = interpolated_CL_charges_common
            var_list1 = var_Q_charge_CL1
            var_list2 = var_Q_charge_CL2
        for i in range(len(C_rates_lct)):
            c_ref1 = 0  # First cycle
            c_ref2 = 9  # 10 th cycle
            name_ref1 = str(C_rates_lct[i]) + "_Q_" + str(c_ref1)
            name_ref2 = str(C_rates_lct[i]) + "_Q_" + str(c_ref2)
            for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                if j >= c_ref1:
                    name_x = str(C_rates_lct[i]) + "_Q_" + str(j)
                    x = dictionary[name_x]
                    x_ref1 = dictionary[name_ref1]
                    x_diff1 = [x[i] - x_ref1[i] for i in range(len(x_ref1))]
                    var_list1.append(statistics.variance(x_diff1))

                if j >= c_ref2:
                    name_x = str(C_rates_lct[i]) + "_Q_" + str(j)
                    x = dictionary[name_x]
                    x_ref2 = dictionary[name_ref2]
                    x_diff2 = [x[i] - x_ref2[i] for i in range(len(x_ref2))]
                    var_list2.append(statistics.variance(x_diff2))

# Averaging end of (dis)charge values
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    aves_end_d = {}
    aves_end_c = {}
    stds_end_d = {}
    stds_end_c = {}
    if formation == 1:
        for d in range(2):
            for k in x_to_plot:
                count_cycle = 0
                ends_of_discharge = []
                ends_of_charge = []
                ave_end_of_discharge = []
                ave_end_of_charge = []
                std_end_of_discharge = []
                std_end_of_charge = []
                for j in range(cycles_per_formations):
                    if d == 0:
                        ends_of_discharge.append(ends_d[k][count_cycle])
                        count_cycle += 1
                    if d == 1:
                        ends_of_charge.append(ends_c[k][count_cycle])
                        count_cycle += 1
                if d == 0:
                    ave = np.mean(ends_of_discharge)
                    std = np.std(ends_of_discharge)
                    ave_end_of_discharge.append(ave)
                    std_end_of_discharge.append(std)
                if d == 1:
                    ave = np.mean(ends_of_charge)
                    std = np.std(ends_of_charge)
                    ave_end_of_charge.append(ave)
                    std_end_of_charge.append(std)
                for i in range(1, len(C_rates_list)):
                    ends_of_discharge = []
                    ends_of_charge = []
                    if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                        cycles_to_count = cycle_per_last
                    else:
                        cycles_to_count = cycles_per_Crates

                    for j in range(cycles_to_count):
                        if d == 0:
                            ends_of_discharge.append(ends_d[k][count_cycle])
                            # print(ends_of_discharge, C_rates_list[i])
                            count_cycle += 1
                        if d == 1:
                            ends_of_charge.append(ends_c[k][count_cycle])
                            count_cycle += 1
                    if d == 0:
                        ave = np.mean(ends_of_discharge)
                        std = np.std(ends_of_discharge)
                        ave_end_of_discharge.append(ave)
                        std_end_of_discharge.append(std)
                        aves_end_d[k] = ave_end_of_discharge
                        stds_end_d[k] = std_end_of_discharge
                    if d == 1:
                        ave = np.mean(ends_of_charge)
                        std = np.std(ends_of_charge)
                        ave_end_of_charge.append(ave)
                        std_end_of_charge.append(std)
                        aves_end_c[k] = ave_end_of_charge
                        stds_end_c[k] = std_end_of_charge

    if formation == 0:
        for d in range(2):
            for k in x_to_plot:
                count_cycle = 0
                ends_of_discharge = []
                ends_of_charge = []
                ave_end_of_discharge = []
                ave_end_of_charge = []
                std_end_of_discharge = []
                std_end_of_charge = []
                for i in range(len(C_rates_list)):
                    ends_of_discharge = []
                    ends_of_charge = []

                    if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                        cycles_to_count = cycle_per_last
                    else:
                        cycles_to_count = cycles_per_Crates

                    for j in range(cycles_to_count):
                        if d == 0:
                            ends_of_discharge.append(ends_d[k][count_cycle])
                            count_cycle += 1
                        if d == 1:
                            ends_of_charge.append(ends_c[k][count_cycle])
                            count_cycle += 1
                    if d == 0:
                        ave = np.mean(ends_of_discharge)
                        std = np.std(ends_of_discharge)
                        ave_end_of_discharge.append(ave)
                        std_end_of_discharge.append(std)
                        aves_end_d[k] = ave_end_of_discharge
                        stds_end_d[k] = std_end_of_discharge
                    if d == 1:
                        ave = np.mean(ends_of_charge)
                        std = np.std(ends_of_charge)
                        ave_end_of_charge.append(ave)
                        std_end_of_charge.append(std)
                        aves_end_c[k] = ave_end_of_charge
                        stds_end_c[k] = std_end_of_charge

# Defining limits for the plottings
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    limits_plotting = {}
    for d in range(2):
        if d == 0:
            direction = "discharge"
        if d == 1:
            direction = "charge"
        for l in y_to_plot:
            limit_up_y = []
            limit_down_y = []
            for k in x_to_plot:
                limit_up_x = []
                limit_down_x = []
                for i in range(len(C_rates_list)):
                    name_x = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_x_" + direction
                    name_y = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_y_" + direction
                    limit_up_x.append(maxs[name_x])
                    limit_down_x.append(mins[name_x])
                    if k == x_to_plot[0]:
                        limit_up_y.append(maxs[name_y])
                        limit_down_y.append(mins[name_y])
                name_x_lim_up = str(k) + "_up_" + direction
                name_x_lim_down = str(k) + "_down_" + direction
                limits_plotting[name_x_lim_up] = limit_up_x
                limits_plotting[name_x_lim_down] = limit_down_x
            name_y_lim_up = str(l) + "_up_" + direction
            name_y_lim_down = str(l) + "_down_" + direction
            limits_plotting[name_y_lim_up] = limit_up_y
            limits_plotting[name_y_lim_down] = limit_down_y

        # Cycle life test(s)
## Problem here is that I am defining a unique series of upper and lower limits for E, Q, and Qn -
# and not n series of limits (where n is len(C_rates_lct)) + I am not distinguishing between charge and discharge
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    maxs_E = []
    mins_E = []
    maxs_Q = []
    mins_Q = []
    maxs_Qn = []
    mins_Qn = []
    for d in range(2):
        if d == 0:
            dictionary = cycle_life_discharges
        if d == 1:
            dictionary = cycle_life_charges
        for i in range(len(C_rates_lct)):
            for j in range(len(num_cycles_dict[C_rates_lct[i]])):
                for k in x_to_plot:
                    name_x = str(C_rates_lct[i]) + "_" + k + "_" + str(j)
                    name_y = str(C_rates_lct[i]) + "_E_" + str(j)
                    x = dictionary[name_x]
                    y = dictionary[name_y]
                    maxs_E.append(max(y))
                    mins_E.append(min(y))
                    if k == "Q":
                        maxs_Q.append(max(x))
                        mins_Q.append(min(x))
                    if k == "Qn":
                        maxs_Qn.append(max(x))
                        mins_Qn.append(min(x))

# Plotting averages - half cycles
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    for d in range(2):
        if d == 0:
            direction = "discharge"
            title = 'Discharges_Csp_half_'
        if d == 1:
            direction = "charge"
            title = 'Charges_Csp_half_'
        for l in y_to_plot:
            for k in x_to_plot:
                titlef = title + str(k)
                fig = plt.figure(figsize=(13, 11))
                plt.title(titlef)
                for i in range(len(C_rates_list)):
                    j = (i + 1) / (len(C_rates_list))
                    x_axis = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_x_" + direction
                    x_axis_minus = averages_minus_std[x_axis]
                    x_axis_plus = averages_plus_std[x_axis]
                    x_axis = averages[x_axis]
                    y_axis = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_y_" + direction
                    y_axis = averages[y_axis]
                    plt.plot(x_axis, y_axis, linewidth=2, color=plt.cm.cool(j), label=str(C_rates_list[i]))
                    plt.fill(x_axis_minus + x_axis_plus[::-1], y_axis + y_axis[::-1], color=plt.cm.cool(j),
                             alpha=0.35)  # seismic
                    plt.legend(fontsize=22)  # loc='lower left'
                    plt.xticks(fontsize=22)
                    plt.yticks(fontsize=22)
                    if k == "Q":
                        plt.xlabel(k + " (mAh)", fontsize=28)  # fontweight="bold"
                        plt.ylabel(l + " (V)", fontsize=28)
                    if k == "Qn":
                        plt.xlabel(k + " (mAh g$^{-1}$)", fontsize=28)  # fontweight="bold"
                        plt.ylabel(l + " (V)", fontsize=28)
                    name_x_lim_up = str(k) + "_up_" + direction
                    name_x_lim_down = str(k) + "_down_" + direction
                    name_y_lim_up = str(l) + "_up_" + direction
                    name_y_lim_down = str(l) + "_down_" + direction
                    x_max = max(limits_plotting[name_x_lim_up])
                    x_min = min(limits_plotting[name_x_lim_down])
                    y_max = max(limits_plotting[name_y_lim_up])
                    y_min = min(limits_plotting[name_y_lim_down])
                    x_min = x_min - (0.05 * x_min)
                    x_max = x_max + (0.05 * x_max)
                    # y_min = y_min - (0.01*y_min)
                    # y_max = y_max + (0.01*y_max)
                    plt.xlim([0, x_max])
                    plt.ylim([y_min, y_max])

                    plot_save = directory_res + "\\" + titlef
                    plt.savefig(plot_save)

# Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for d in range(2):
        if d == 0:
            dictionary = cycle_life_discharges
            title = 'Discharges_Csp_half_'
        if d == 1:
            dictionary = cycle_life_charges
            title = 'Charges_Csp_half_'
        for k in x_to_plot:
            titlef = title + str(k)
            fig = plt.figure(figsize=(13, 11))
            plt.title(titlef)
            for i in range(len(C_rates_lct)):
                range_cycles = len(num_cycles_dict[C_rates_lct[i]])
                every_plot = 5
                if range_cycles >= 50:
                    every_plot = 10
                if range_cycles >= 100:
                    every_plot = 20
                if range_cycles >= 200:
                    every_plot = 30
                if range_cycles >= 300:
                    every_plot = 50
                if range_cycles >= 1000:
                    every_plot = 100
                for j in range(range_cycles):
                    c = j / (range_cycles)
                    name_x = str(C_rates_lct[i]) + "_" + k + "_" + str(j)
                    name_y = str(C_rates_lct[i]) + "_E_" + str(j)
                    x = dictionary[name_x]
                    y = dictionary[name_y]
                    if j == 0:
                        plt.plot(x, y, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.legend(fontsize=22)  # loc='lower left'
                        plt.xticks(fontsize=22)
                        plt.yticks(fontsize=22)
                    for n in range(100):
                        if j + 1 == every_plot * n:
                            plt.plot(x, y, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                            plt.legend(fontsize=22)  # loc='lower left'
                            plt.xticks(fontsize=22)
                            plt.yticks(fontsize=22)
                    if j == range_cycles - 1:
                        plt.plot(x, y, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.legend(fontsize=22)  # loc='lower left'
                        plt.xticks(fontsize=22)
                        plt.yticks(fontsize=22)
            if k == "Q":
                plt.xlabel(k + " (mAh)", fontsize=28)  # fontweight="bold"
                plt.ylabel("E (V)", fontsize=28)
                plt.xlim([0, max(maxs_Q)])
            if k == "Qn":
                plt.xlabel(k + " (mAh g$^{-1}$)", fontsize=28)
                plt.ylabel("E (V)", fontsize=28)
                plt.xlim([0, max(maxs_Qn)])
            plt.ylim([min(mins_E), max(maxs_E)])
            plot_save = directory_res_CL + "\\" + titlef
            plt.savefig(plot_save)

        # Plotting averages - full cycles
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    for l in y_to_plot:
        for k in x_to_plot:
            title = 'Full_cycle_'
            titlef = title + str(k)
            fig = plt.figure(figsize=(13, 11))
            for i in range(len(C_rates_list)):
                j = (i + 1) / (len(C_rates_list))
                x_axis_d = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_x_discharge"
                x_axis_minus_d = averages_minus_std[x_axis_d]
                x_axis_plus_d = averages_plus_std[x_axis_d]
                x_axis_d = averages[x_axis_d]
                y_axis_d = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_y_discharge"
                y_axis_d = averages[y_axis_d]
                plt.plot(x_axis_d, y_axis_d, linewidth=2, color=plt.cm.cool(j), label=str(C_rates_list[i]))
                plt.fill(x_axis_minus_d + x_axis_plus_d[::-1], y_axis_d + y_axis_d[::-1], color=plt.cm.cool(j),
                         alpha=0.35)  # seismic
                x_axis_c = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_x_charge"
                x_axis_minus_c = averages_minus_std[x_axis_c]
                x_axis_plus_c = averages_plus_std[x_axis_c]
                x_axis_c = averages[x_axis_c]
                y_axis_c = str(C_rates_list[i]) + "_" + str(k) + "_" + str(l) + "_y_charge"
                y_axis_c = averages[y_axis_c]
                plt.plot(x_axis_c, y_axis_c, linewidth=2, color=plt.cm.cool(j))
                plt.fill(x_axis_minus_c + x_axis_plus_c[::-1], y_axis_c + y_axis_c[::-1], color=plt.cm.cool(j),
                         alpha=0.35)  # seismic
                plt.legend(fontsize=22)  # loc='lower left'
                plt.xticks(fontsize=22)
                plt.yticks(fontsize=22)
                if k == "Q":
                    plt.xlabel(k + " (mAh)", fontsize=28)  # fontweight="bold"
                    plt.ylabel(l + " (V)", fontsize=28)
                if k == "Qn":
                    plt.xlabel(k + " (mAh g$^{-1}$)", fontsize=28)  # fontweight="bold"
                    plt.ylabel(l + " (V)", fontsize=28)
                name_x_lim_up_d = str(k) + "_up_discharge"
                name_x_lim_down_d = str(k) + "_down_discharge"
                name_x_lim_up_c = str(k) + "_up_charge"
                name_x_lim_down_c = str(k) + "_down_charge"
                name_y_lim_up = str(l) + "_up_discharge"
                name_y_lim_down = str(l) + "_down_discharge"
                x_max_d = max(limits_plotting[name_x_lim_up_d])
                x_min_d = min(limits_plotting[name_x_lim_down_d])
                x_max_c = max(limits_plotting[name_x_lim_up_c])
                x_min_c = min(limits_plotting[name_x_lim_down_c])
                x_max = max(x_max_d, x_max_c)
                x_min = min(x_min_d, x_min_c)
                y_max = max(limits_plotting[name_y_lim_up])
                y_min = min(limits_plotting[name_y_lim_down])
                x_min = x_min - (0.05 * x_min)
                x_max = x_max + (0.05 * x_max)
                # y_min = y_min - (0.01*y_min)
                # y_max = y_max + (0.01*y_max)
                plt.xlim([0, x_max])
                plt.ylim([y_min, y_max])
                plt.title(titlef)

                plot_save = directory_res + "\\" + titlef
                plt.savefig(plot_save)

# Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    title = 'Full_cycle_'
    for k in x_to_plot:
        titlef = title + str(k)
        fig = plt.figure(figsize=(13, 11))
        plt.title(titlef)
        for i in range(len(C_rates_lct)):
            range_cycles = len(num_cycles_dict[C_rates_lct[i]])
            every_plot = 5
            if range_cycles >= 50:
                every_plot = 10
            if range_cycles >= 100:
                every_plot = 20
            if range_cycles >= 200:
                every_plot = 30
            if range_cycles >= 300:
                every_plot = 50
            if range_cycles >= 1000:
                every_plot = 100
            for j in range(range_cycles):
                c = j / (range_cycles)
                name_x = str(C_rates_lct[i]) + "_" + k + "_" + str(j)
                name_y = str(C_rates_lct[i]) + "_E_" + str(j)
                xd = cycle_life_discharges[name_x]
                yd = cycle_life_discharges[name_y]
                xc = cycle_life_charges[name_x]
                yc = cycle_life_charges[name_y]
                if j == 0:
                    plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                    plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                    plt.legend(fontsize=22)  # loc='lower left'
                    plt.xticks(fontsize=22)
                    plt.yticks(fontsize=22)
                for n in range(100):
                    if j + 1 == every_plot * n:
                        plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                        plt.legend(fontsize=22)  # loc='lower left'
                        plt.xticks(fontsize=22)
                        plt.yticks(fontsize=22)
                if j == range_cycles - 1:
                    plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                    plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                    plt.legend(fontsize=22)  # loc='lower left'
                    plt.xticks(fontsize=22)
                    plt.yticks(fontsize=22)
        if k == "Q":
            plt.xlabel(k + " (mAh)", fontsize=28)  # fontweight="bold"
            plt.ylabel("E (V)", fontsize=28)
            plt.xlim([0, max(maxs_Q)])
        if k == "Qn":
            plt.xlabel(k + " (mAh g$^{-1}$)", fontsize=28)
            plt.ylabel("E (V)", fontsize=28)
            plt.xlim([0, max(maxs_Qn)])
        plt.ylim([min(mins_E), max(maxs_E)])
        plot_save = directory_res_CL + "\\" + titlef
        plt.savefig(plot_save)

# Cutting the extremes in the dQdV curves
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    safe_remove = 0.05
    max_up_E_charge = max(limits_plotting["E_up_charge"])
    max_up_E_discharge = max(limits_plotting["E_up_discharge"])
    max_up_E_common = min(max_up_E_charge, max_up_E_discharge) - safe_remove
    min_down_E_charge = min(limits_plotting["E_down_charge"])
    min_down_E_discharge = min(limits_plotting["E_down_discharge"])
    min_down_E_common = max(min_down_E_charge, min_down_E_discharge) + safe_remove

    temp_E = []
    temp_dQdV = []
    if formation == 0:
        for i in range(len(C_rates_list)):

            if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                cycles_to_count = cycle_per_last
            else:
                cycles_to_count = cycles_per_Crates

            for j in range(cycles_to_count):
                name_dict_y = str(C_rates_list[i]) + "_dQdV_" + str(j + 1)
                name_dict_x = str(C_rates_list[i]) + "_E_" + str(j + 1)

                for k in range(len(dQdV_dict_c[name_dict_x])):
                    tester = dQdV_dict_c[name_dict_x][k]
                    if tester > min_down_E_common and tester < max_up_E_common:
                        temp_dQdV.append(dQdV_dict_c[name_dict_y][k])
                        temp_E.append(dQdV_dict_c[name_dict_x][k])
                dQdV_dict_c[name_dict_y] = temp_dQdV
                dQdV_dict_c[name_dict_x] = temp_E
                temp_E = []
                temp_dQdV = []

                for h in range(len(dQdV_dict_d[name_dict_x])):
                    tester = dQdV_dict_d[name_dict_x][h]
                    if tester > min_down_E_common and tester < max_up_E_common:
                        temp_dQdV.append(dQdV_dict_d[name_dict_y][h])
                        temp_E.append(dQdV_dict_d[name_dict_x][h])
                dQdV_dict_d[name_dict_y] = temp_dQdV
                dQdV_dict_d[name_dict_x] = temp_E
                temp_E = []
                temp_dQdV = []

    if formation == 1:
        for i in range(len(C_rates_list)):

            if C_rates_list[i] == C_rates_list[0]:
                cycles_to_count = cycles_per_formations
            elif finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                cycles_to_count = cycle_per_last
            else:
                cycles_to_count = cycles_per_Crates

            for j in range(cycles_to_count):
                name_dict_y = str(C_rates_list[i]) + "_dQdV_" + str(j + 1)
                name_dict_x = str(C_rates_list[i]) + "_E_" + str(j + 1)

                for k in range(len(dQdV_dict_c[name_dict_x])):
                    tester = dQdV_dict_c[name_dict_x][k]
                    if tester > min_down_E_common and tester < max_up_E_common:
                        temp_dQdV.append(dQdV_dict_c[name_dict_y][k])
                        temp_E.append(dQdV_dict_c[name_dict_x][k])
                dQdV_dict_c[name_dict_y] = temp_dQdV
                dQdV_dict_c[name_dict_x] = temp_E
                temp_E = []
                temp_dQdV = []

                for h in range(len(dQdV_dict_d[name_dict_x])):
                    tester = dQdV_dict_d[name_dict_x][h]
                    if tester > min_down_E_common and tester < max_up_E_common:
                        temp_dQdV.append(dQdV_dict_d[name_dict_y][h])
                        temp_E.append(dQdV_dict_d[name_dict_x][h])
                dQdV_dict_d[name_dict_y] = temp_dQdV
                dQdV_dict_d[name_dict_x] = temp_E
                temp_E = []
                temp_dQdV = []

            # Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for i in range(len(C_rates_lct)):
        range_cycles = len(num_cycles_dict[C_rates_lct[i]])
        for j in range(range_cycles):
            name_x = str(C_rates_lct[i]) + "_E_" + str(j)
            name_y = str(C_rates_lct[i]) + "_dQdV_" + str(j)
            for l in range(2):
                temp_E = []
                temp_dQdV = []
                if l == 0:
                    dictionary = dQdv_CL_c
                if l == 1:
                    dictionary = dQdv_CL_d
                for k in range(len(dictionary[name_x])):
                    tester = dictionary[name_x][k]
                    if tester > min_down_E_common and tester < max_up_E_common:
                        temp_dQdV.append(dictionary[name_y][k])
                        temp_E.append(dictionary[name_x][k])
                dictionary[name_y] = temp_dQdV
                dictionary[name_x] = temp_E

# Cutting the extremes in the dVdQ curves
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    safe_remove = 4
    max_up_Q = 100 - safe_remove / 2
    min_up_Q = 0 + safe_remove / 2

    temp_Q = []
    temp_dVdQ = []
    if formation == 0:
        for i in range(len(C_rates_list)):

            if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                cycles_to_count = cycle_per_last
            else:
                cycles_to_count = cycles_per_Crates

            for j in range(cycles_to_count):
                name_dict_y = str(C_rates_list[i]) + "_dVdQ_" + str(j + 1)
                name_dict_x = str(C_rates_list[i]) + "_Q_" + str(j + 1)

                for k in range(len(dVdQ_dict_c[name_dict_x])):
                    tester = dVdQ_dict_c[name_dict_x][k]
                    if tester > min_up_Q and tester < max_up_Q:
                        temp_dVdQ.append(dVdQ_dict_c[name_dict_y][k])
                        temp_Q.append(dVdQ_dict_c[name_dict_x][k])
                dVdQ_dict_c[name_dict_y] = temp_dVdQ
                dVdQ_dict_c[name_dict_x] = temp_Q
                temp_Q = []
                temp_dVdQ = []

                for h in range(len(dVdQ_dict_d[name_dict_x])):
                    tester = dVdQ_dict_d[name_dict_x][h]
                    if tester > min_up_Q and tester < max_up_Q:
                        temp_dVdQ.append(dVdQ_dict_d[name_dict_y][h])
                        temp_Q.append(dVdQ_dict_d[name_dict_x][h])
                dVdQ_dict_d[name_dict_y] = temp_dVdQ
                dVdQ_dict_d[name_dict_x] = temp_Q
                temp_Q = []
                temp_dVdQ = []

    if formation == 1:
        for i in range(len(C_rates_list)):

            if C_rates_list[i] == C_rates_list[0]:
                cycles_to_count = cycles_per_formations
            elif finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                cycles_to_count = cycle_per_last
            else:
                cycles_to_count = cycles_per_Crates

            for j in range(cycles_to_count):
                name_dict_y = str(C_rates_list[i]) + "_dVdQ_" + str(j + 1)
                name_dict_x = str(C_rates_list[i]) + "_Q_" + str(j + 1)

                for k in range(len(dVdQ_dict_c[name_dict_x])):
                    tester = dVdQ_dict_c[name_dict_x][k]
                    if tester > min_up_Q and tester < max_up_Q:
                        temp_dVdQ.append(dVdQ_dict_c[name_dict_y][k])
                        temp_Q.append(dVdQ_dict_c[name_dict_x][k])
                dVdQ_dict_c[name_dict_y] = temp_dVdQ
                dVdQ_dict_c[name_dict_x] = temp_Q
                temp_Q = []
                temp_dVdQ = []

                for h in range(len(dVdQ_dict_d[name_dict_x])):
                    tester = dVdQ_dict_d[name_dict_x][h]
                    if tester > min_up_Q and tester < max_up_Q:
                        temp_dVdQ.append(dVdQ_dict_d[name_dict_y][h])
                        temp_Q.append(dVdQ_dict_d[name_dict_x][h])
                dVdQ_dict_d[name_dict_y] = temp_dVdQ
                dVdQ_dict_d[name_dict_x] = temp_Q
                temp_Q = []
                temp_dVdQ = []

            # Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for i in range(len(C_rates_lct)):
        range_cycles = len(num_cycles_dict[C_rates_lct[i]])
        for j in range(range_cycles):
            name_x = str(C_rates_lct[i]) + "_Q_" + str(j)
            name_y = str(C_rates_lct[i]) + "_dVdQ_" + str(j)
            for l in range(2):
                temp_Q = []
                temp_dVdQ = []
                if l == 0:
                    dictionary = dVdQ_CL_c
                if l == 1:
                    dictionary = dVdQ_CL_d
                for k in range(len(dictionary[name_x])):
                    tester = dictionary[name_x][k]
                    if tester > min_up_Q and tester < max_up_Q:
                        temp_dVdQ.append(dictionary[name_y][k])
                        temp_Q.append(dictionary[name_x][k])
                dictionary[name_y] = temp_dVdQ
                dictionary[name_x] = temp_Q

# Plotting dQdV curves with peaks' position
# Rate capability test
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    peaks_der_c = {}
    peaks_der_c_w = {}
    peaks_der_d = {}
    peaks_der_d_w = {}
    fact_height = 1000  # 1000
    fact_prominence = 20  # 20
    fact_distance = 10  # 20
    min_width = 2
    for plots in range(2):
        if plots == 0:
            titlef = 'dQdV_curves'
        if plots == 1:
            titlef = 'dQdV_curves + peaks'
        fig = plt.figure(figsize=(13, 11))
        counter = 1
        # y_max = 1000000
        # y_min = -1000000
        # x_max = 1000000
        # x_min = -1000000
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        plt.xlabel("E (V)", fontsize=28)
        plt.ylabel("dQ/dV (mAh V$^{-1}$)", fontsize=28)
        plt.title(titlef, fontsize=24)

        if formation == 1:
            for j in range(cycles_per_formations):
                col_index = counter / number_of_cycles
                name_dict_y = str(C_rates_list[0]) + "_dQdV_" + str(j + 1)
                name_dict_x = str(C_rates_list[0]) + "_E_" + str(j + 1)
                xc, xd = dQdV_dict_c[name_dict_x], dQdV_dict_d[name_dict_x]
                yc, yd = dQdV_dict_c[name_dict_y], dQdV_dict_d[name_dict_y]
                lab = str(C_rates_list[0]) + "_" + str(counter) + " [f]"
                if plots == 0:
                    plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                    plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                if plots == 1:
                    # Peaks identification on charge
                    peaks_c = find_peaks(yc, height=0.01)  # prominence=100
                    # print(peaks_c[0])
                    try:
                        x_pos_c = [xc[i] for i in peaks_c[0]]
                        x_heights_c = peaks_c[1]["peak_heights"]
                        min_height_c = max(x_heights_c) / fact_height
                        min_prominence_c = max(x_heights_c) / fact_prominence
                        min_disance_c = fact_distance  # len(yc)/fact_distance
                        peaks_c = find_peaks(yc, height=min_height_c, distance=min_disance_c,
                                             prominence=min_prominence_c, width=min_width)
                    except:
                        continue
                    # Peaks identification on discharge
                    yd_rev = [i * -1 for i in yd]
                    peaks_d = find_peaks(yd_rev, height=0.01)  # prominence=100
                    try:
                        x_pos_d = [xd[i] for i in peaks_d[0]]
                        x_heights_d = peaks_d[1]["peak_heights"]
                        min_height_d = max(x_heights_d) / fact_height
                        min_prominence_d = max(x_heights_d) / fact_prominence
                        min_disance_d = fact_distance  # len(yc)/fact_distance
                        peaks_d = find_peaks(yd_rev, height=min_height_d, distance=min_disance_d,
                                             prominence=min_prominence_d, width=min_width)
                    except:
                        continue
                    # Saving in dictionaries
                    peaks_der_c[name_dict_x] = [xc[i] for i in peaks_c[0]]
                    peaks_der_c[name_dict_y] = [yc[i] for i in peaks_c[0]]
                    peaks_der_c_w[name_dict_y] = [i for i in peaks_c[1]['widths']]
                    peaks_der_d[name_dict_x] = [xd[i] for i in peaks_d[0]]
                    peaks_der_d[name_dict_y] = [yd[i] for i in peaks_d[0]]
                    peaks_der_d_w[name_dict_y] = [i for i in peaks_d[1]['widths']]
                    # Plotting
                    lab = str(C_rates_list[0]) + "_" + str(counter) + " [f]"
                    plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                    plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                    plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o',
                             color=plt.cm.cool(col_index), markersize=7)
                    plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o',
                             color=plt.cm.cool(col_index), markersize=7)
                counter += 1

            for i in range(1, len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last
                else:
                    cycles_to_count = cycles_per_Crates

                for j in range(cycles_to_count):
                    name_dict_y = str(C_rates_list[i]) + "_dQdV_" + str(j + 1)
                    name_dict_x = str(C_rates_list[i]) + "_E_" + str(j + 1)
                    xc, xd = dQdV_dict_c[name_dict_x], dQdV_dict_d[name_dict_x]
                    yc, yd = dQdV_dict_c[name_dict_y], dQdV_dict_d[name_dict_y]
                    # Peaks identification on charge
                    peaks_c = find_peaks(yc, height=0.01)  # prominence=100
                    # print(peaks_c[0])
                    try:
                        x_pos_c = [xc[i] for i in peaks_c[0]]
                        x_heights_c = peaks_c[1]["peak_heights"]
                        min_height_c = max(x_heights_c) / fact_height
                        min_prominence_c = max(x_heights_c) / fact_prominence
                        min_disance_c = fact_distance  # len(yc)/fact_distance
                        peaks_c = find_peaks(yc, height=min_height_c, distance=min_disance_c,
                                             prominence=min_prominence_c, width=min_width)
                    except:
                        continue
                    # Peaks identification on discharge
                    yd_rev = [i * -1 for i in yd]
                    peaks_d = find_peaks(yd_rev, height=0.01)  # prominence=100
                    try:
                        x_pos_d = [xd[i] for i in peaks_d[0]]
                        x_heights_d = peaks_d[1]["peak_heights"]
                        min_height_d = max(x_heights_d) / fact_height
                        min_prominence_d = max(x_heights_d) / fact_prominence
                        min_disance_d = fact_distance  # len(yc)/fact_distance
                        peaks_d = find_peaks(yd_rev, height=min_height_d, distance=min_disance_d,
                                             prominence=min_prominence_d, width=min_width)
                    except:
                        continue
                    # Saving in dictionaries
                    peaks_der_c[name_dict_x] = [xc[i] for i in peaks_c[0]]
                    peaks_der_c[name_dict_y] = [yc[i] for i in peaks_c[0]]
                    peaks_der_c_w[name_dict_y] = [i for i in peaks_c[1]['widths']]
                    peaks_der_d[name_dict_x] = [xd[i] for i in peaks_d[0]]
                    peaks_der_d[name_dict_y] = [yd[i] for i in peaks_d[0]]
                    peaks_der_d_w[name_dict_y] = [i for i in peaks_d[1]['widths']]
                    if j == 0 or j == cycles_to_count - 1:
                        col_index = counter / number_of_cycles
                        # print(col_index)
                        xc, xd = dQdV_dict_c[name_dict_x], dQdV_dict_d[name_dict_x]
                        yc, yd = dQdV_dict_c[name_dict_y], dQdV_dict_d[name_dict_y]
                        lab = str(C_rates_list[i]) + "_" + str(counter)
                        if plots == 0:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                        if plots == 1:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                            plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o',
                                     color=plt.cm.cool(col_index), markersize=7)
                            plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o',
                                     color=plt.cm.cool(col_index), markersize=7)
                    counter += 1

        if formation == 0:
            for i in range(len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last
                else:
                    cycles_to_count = cycles_per_Crates

                for j in range(cycles_to_count):
                    name_dict_y = str(C_rates_list[i]) + "_dQdV_" + str(j + 1)
                    name_dict_x = str(C_rates_list[i]) + "_E_" + str(j + 1)
                    xc, xd = dQdV_dict_c[name_dict_x], dQdV_dict_d[name_dict_x]
                    yc, yd = dQdV_dict_c[name_dict_y], dQdV_dict_d[name_dict_y]
                    # Peaks identification on charge
                    peaks_c = find_peaks(yc, height=0.01)  # prominence=100
                    # print(peaks_c[0])
                    try:
                        x_pos_c = [xc[i] for i in peaks_c[0]]
                        x_heights_c = peaks_c[1]["peak_heights"]
                        min_height_c = max(x_heights_c) / fact_height
                        min_prominence_c = max(x_heights_c) / fact_prominence
                        min_disance_c = fact_distance  # len(yc)/fact_distance
                        peaks_c = find_peaks(yc, height=min_height_c, distance=min_disance_c,
                                             prominence=min_prominence_c, width=min_width)
                    except:
                        continue
                    # Peaks identification on discharge
                    yd_rev = [i * -1 for i in yd]
                    peaks_d = find_peaks(yd_rev, height=0.01)  # prominence=100
                    try:
                        x_pos_d = [xd[i] for i in peaks_d[0]]
                        x_heights_d = peaks_d[1]["peak_heights"]
                        min_height_d = max(x_heights_d) / fact_height
                        min_prominence_d = max(x_heights_d) / fact_prominence
                        min_disance_d = fact_distance  # len(yc)/fact_distance
                        peaks_d = find_peaks(yd_rev, height=min_height_d, distance=min_disance_d,
                                             prominence=min_prominence_d, width=min_width)
                    except:
                        continue
                    # Saving in dictionaries
                    peaks_der_c[name_dict_x] = [xc[i] for i in peaks_c[0]]
                    peaks_der_c[name_dict_y] = [yc[i] for i in peaks_c[0]]
                    peaks_der_c_w[name_dict_y] = [i for i in peaks_c[1]['widths']]
                    peaks_der_d[name_dict_x] = [xd[i] for i in peaks_d[0]]
                    peaks_der_d[name_dict_y] = [yd[i] for i in peaks_d[0]]
                    peaks_der_d_w[name_dict_y] = [i for i in peaks_d[1]['widths']]
                    if j == 0 or j == cycles_to_count - 1:
                        col_index = counter / number_of_cycles
                        xc, xd = dQdV_dict_c[name_dict_x], dQdV_dict_d[name_dict_x]
                        yc, yd = dQdV_dict_c[name_dict_y], dQdV_dict_d[name_dict_y]
                        lab = str(C_rates_list[i]) + "_" + str(counter)
                        if plots == 0:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                        if plots == 1:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                            plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o',
                                     color=plt.cm.cool(col_index), markersize=7)
                            plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o',
                                     color=plt.cm.cool(col_index), markersize=7)
                    counter += 1

        plt.legend(fontsize=16)  # bbox_to_anchor=(0.986,1)
        plt.xlim([min_down_E_common, max_up_E_common])
        # plt.ylim([-3,3])

        plot_save = directory_res + "\\" + titlef
        plt.savefig(plot_save)

# Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    peaks_derCL_c = {}
    peaks_derCL_c_w = {}
    peaks_derCL_d = {}
    peaks_derCL_d_w = {}
    fact_height = 1000  # 1000
    fact_prominence = 20  # 20
    fact_distance = 10  # 20
    min_width = 2
    for plots in range(2):
        if plots == 0:
            titlef = 'dQdV_curves'
        if plots == 1:
            titlef = 'dQdV_curves + peaks'
        fig = plt.figure(figsize=(13, 11))
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        plt.xlabel("E (V)", fontsize=28)
        plt.ylabel("dQ/dV (mAh V$^{-1}$)", fontsize=28)
        plt.title(titlef, fontsize=24)
        for i in range(len(C_rates_lct)):
            range_cycles = len(num_cycles_dict[C_rates_lct[i]])
            every_plot = 5
            if range_cycles >= 50:
                every_plot = 10
            if range_cycles >= 100:
                every_plot = 20
            if range_cycles >= 200:
                every_plot = 30
            if range_cycles >= 300:
                every_plot = 50
            if range_cycles >= 1000:
                every_plot = 100
            for j in range(range_cycles):
                c = j / (range_cycles)
                name_x = str(C_rates_lct[i]) + "_E_" + str(j)
                name_y = str(C_rates_lct[i]) + "_dQdV_" + str(j)
                xc, xd = dQdv_CL_c[name_x], dQdv_CL_d[name_x]
                yc, yd = dQdv_CL_c[name_y], dQdv_CL_d[name_y]
                if plots == 0:
                    if j == 0:
                        plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                    for n in range(100):
                        if j + 1 == every_plot * n:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                    if j == range_cycles - 1:
                        plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                if plots == 1:
                    # Peaks identification on charge
                    peaks_c = find_peaks(yc, height=0.01)  # prominence=100
                    # print(peaks_c[0])
                    try:
                        x_pos_c = [xc[i] for i in peaks_c[0]]
                        x_heights_c = peaks_c[1]["peak_heights"]
                        min_height_c = max(x_heights_c) / fact_height
                        min_prominence_c = max(x_heights_c) / fact_prominence
                        min_disance_c = fact_distance  # len(yc)/fact_distance
                        peaks_c = find_peaks(yc, height=min_height_c, distance=min_disance_c,
                                             prominence=min_prominence_c, width=min_width)
                    except:
                        continue
                    # Peaks identification on discharge
                    yd_rev = [i * -1 for i in yd]
                    peaks_d = find_peaks(yd_rev, height=0.01)  # prominence=100
                    try:
                        x_pos_d = [xd[i] for i in peaks_d[0]]
                        x_heights_d = peaks_d[1]["peak_heights"]
                        min_height_d = max(x_heights_d) / fact_height
                        min_prominence_d = max(x_heights_d) / fact_prominence
                        min_disance_d = fact_distance  # len(yc)/fact_distance
                        peaks_d = find_peaks(yd_rev, height=min_height_d, distance=min_disance_d,
                                             prominence=min_prominence_d, width=min_width)
                    except:
                        continue
                    # Saving in dictionaries
                    peaks_derCL_c[name_x] = [xc[i] for i in peaks_c[0]]
                    peaks_derCL_c[name_y] = [yc[i] for i in peaks_c[0]]
                    peaks_derCL_c_w[name_y] = [i for i in peaks_c[1]['widths']]
                    peaks_derCL_d[name_x] = [xd[i] for i in peaks_d[0]]
                    peaks_derCL_d[name_y] = [yd[i] for i in peaks_d[0]]
                    peaks_derCL_d_w[name_y] = [i for i in peaks_d[1]['widths']]
                    # Plotting
                    if j == 0:
                        plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                        plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o', color=plt.cm.cool(c),
                                 markersize=7)
                        plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o', color=plt.cm.cool(c),
                                 markersize=7)
                    for n in range(100):
                        if j + 1 == every_plot * n:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                            plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o',
                                     color=plt.cm.cool(c), markersize=7)
                            plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o',
                                     color=plt.cm.cool(c), markersize=7)
                    if j == range_cycles - 1:
                        plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                        plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o', color=plt.cm.cool(c),
                                 markersize=7)
                        plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o', color=plt.cm.cool(c),
                                 markersize=7)

            # plt.ylim([min(mins_E)+safe_remove, max(maxs_E)-safe_remove])
            plt.legend(fontsize=16)
            plot_save = directory_res_CL + "\\" + titlef
            plt.savefig(plot_save)

# Plotting dVdQ curves with peaks' position
# Rate capability test
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    peaks_der_c1 = {}
    peaks_der_c1_w = {}
    peaks_der_d1 = {}
    peaks_der_d1_w = {}
    fact_height = 1000  # 1000
    fact_prominence = 20  # 20
    fact_distance = 10  # 20
    min_width = 2
    for plots in range(2):
        if plots == 0:
            titlef = 'dVdQ_curves'
        if plots == 1:
            titlef = 'dVdQ_curves + peaks'
        fig = plt.figure(figsize=(13, 11))
        counter = 1
        # y_max = 1000000
        # y_min = -1000000
        # x_max = 1000000
        # x_min = -1000000
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        plt.xlabel("Q (mAh)", fontsize=28)
        plt.ylabel("dV/dQ (V mAh$^{-1}$)", fontsize=28)
        plt.title(titlef, fontsize=24)

        if formation == 1:
            for j in range(cycles_per_formations):
                col_index = counter / number_of_cycles
                name_dict_y = str(C_rates_list[0]) + "_dVdQ_" + str(j + 1)
                name_dict_x = str(C_rates_list[0]) + "_Q_" + str(j + 1)
                xc, xd = dVdQ_dict_c[name_dict_x], dVdQ_dict_d[name_dict_x]
                yc, yd = dVdQ_dict_c[name_dict_y], dVdQ_dict_d[name_dict_y]
                lab = str(C_rates_list[0]) + "_" + str(counter) + " [f]"
                if plots == 0:
                    plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                    plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                if plots == 1:
                    # Peaks identification on charge
                    peaks_c = find_peaks(yc, height=0.01)  # prominence=100
                    # print(peaks_c[0])
                    try:
                        x_pos_c = [xc[i] for i in peaks_c[0]]
                        x_heights_c = peaks_c[1]["peak_heights"]
                        min_height_c = max(x_heights_c) / fact_height
                        min_prominence_c = max(x_heights_c) / fact_prominence
                        min_disance_c = fact_distance  # len(yc)/fact_distance
                        peaks_c = find_peaks(yc, height=min_height_c, distance=min_disance_c,
                                             prominence=min_prominence_c, width=min_width)
                    except:
                        continue
                    # Peaks identification on discharge
                    yd_rev = [i * -1 for i in yd]
                    peaks_d = find_peaks(yd_rev, height=0.01)  # prominence=100
                    try:
                        x_pos_d = [xd[i] for i in peaks_d[0]]
                        x_heights_d = peaks_d[1]["peak_heights"]
                        min_height_d = max(x_heights_d) / fact_height
                        min_prominence_d = max(x_heights_d) / fact_prominence
                        min_disance_d = fact_distance  # len(yc)/fact_distance
                        peaks_d = find_peaks(yd_rev, height=min_height_d, distance=min_disance_d,
                                             prominence=min_prominence_d, width=min_width)
                    except:
                        continue
                    # Saving in dictionaries
                    peaks_der_c1[name_dict_x] = [xc[i] for i in peaks_c[0]]
                    peaks_der_c1[name_dict_y] = [yc[i] for i in peaks_c[0]]
                    peaks_der_c1_w[name_dict_y] = [i for i in peaks_c[1]['widths']]
                    peaks_der_d1[name_dict_x] = [xd[i] for i in peaks_d[0]]
                    peaks_der_d1[name_dict_y] = [yd[i] for i in peaks_d[0]]
                    peaks_der_d1_w[name_dict_y] = [i for i in peaks_d[1]['widths']]
                    # Plotting
                    lab = str(C_rates_list[0]) + "_" + str(counter) + " [f]"
                    plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                    plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                    plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o',
                             color=plt.cm.cool(col_index), markersize=7)
                    plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o',
                             color=plt.cm.cool(col_index), markersize=7)
                counter += 1

            for i in range(1, len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last
                else:
                    cycles_to_count = cycles_per_Crates

                for j in range(cycles_to_count):
                    name_dict_y = str(C_rates_list[i]) + "_dVdQ_" + str(j + 1)
                    name_dict_x = str(C_rates_list[i]) + "_Q_" + str(j + 1)
                    xc, xd = dVdQ_dict_c[name_dict_x], dVdQ_dict_d[name_dict_x]
                    yc, yd = dVdQ_dict_c[name_dict_y], dVdQ_dict_d[name_dict_y]
                    # Peaks identification on charge
                    peaks_c = find_peaks(yc, height=0.01)  # prominence=100
                    # print(peaks_c[0])
                    try:
                        x_pos_c = [xc[i] for i in peaks_c[0]]
                        x_heights_c = peaks_c[1]["peak_heights"]
                        min_height_c = max(x_heights_c) / fact_height
                        min_prominence_c = max(x_heights_c) / fact_prominence
                        min_disance_c = fact_distance  # len(yc)/fact_distance
                        peaks_c = find_peaks(yc, height=min_height_c, distance=min_disance_c,
                                             prominence=min_prominence_c, width=min_width)
                    except:
                        continue
                    # Peaks identification on discharge
                    yd_rev = [i * -1 for i in yd]
                    peaks_d = find_peaks(yd_rev, height=0.01)  # prominence=100
                    try:
                        x_pos_d = [xd[i] for i in peaks_d[0]]
                        x_heights_d = peaks_d[1]["peak_heights"]
                        min_height_d = max(x_heights_d) / fact_height
                        min_prominence_d = max(x_heights_d) / fact_prominence
                        min_disance_d = fact_distance  # len(yc)/fact_distance
                        peaks_d = find_peaks(yd_rev, height=min_height_d, distance=min_disance_d,
                                             prominence=min_prominence_d, width=min_width)
                    except:
                        continue
                    # Saving in dictionaries
                    peaks_der_c1[name_dict_x] = [xc[i] for i in peaks_c[0]]
                    peaks_der_c1[name_dict_y] = [yc[i] for i in peaks_c[0]]
                    peaks_der_c1_w[name_dict_y] = [i for i in peaks_c[1]['widths']]
                    peaks_der_d1[name_dict_x] = [xd[i] for i in peaks_d[0]]
                    peaks_der_d1[name_dict_y] = [yd[i] for i in peaks_d[0]]
                    peaks_der_d1_w[name_dict_y] = [i for i in peaks_d[1]['widths']]
                    if j == 0 or j == cycles_to_count - 1:
                        col_index = counter / number_of_cycles
                        # print(col_index)
                        xc, xd = dVdQ_dict_c[name_dict_x], dVdQ_dict_d[name_dict_x]
                        yc, yd = dVdQ_dict_c[name_dict_y], dVdQ_dict_d[name_dict_y]
                        lab = str(C_rates_list[i]) + "_" + str(counter)
                        if plots == 0:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                        if plots == 1:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                            plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o',
                                     color=plt.cm.cool(col_index), markersize=7)
                            plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o',
                                     color=plt.cm.cool(col_index), markersize=7)
                    counter += 1

        if formation == 0:
            for i in range(len(C_rates_list)):

                if finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                    cycles_to_count = cycle_per_last
                else:
                    cycles_to_count = cycles_per_Crates

                for j in range(cycles_to_count):
                    name_dict_y = str(C_rates_list[i]) + "_dVdQ_" + str(j + 1)
                    name_dict_x = str(C_rates_list[i]) + "_Q_" + str(j + 1)
                    xc, xd = dQdV_dict_c[name_dict_x], dQdV_dict_d[name_dict_x]
                    yc, yd = dQdV_dict_c[name_dict_y], dQdV_dict_d[name_dict_y]
                    # Peaks identification on charge
                    peaks_c = find_peaks(yc, height=0.01)  # prominence=100
                    # print(peaks_c[0])
                    try:
                        x_pos_c = [xc[i] for i in peaks_c[0]]
                        x_heights_c = peaks_c[1]["peak_heights"]
                        min_height_c = max(x_heights_c) / fact_height
                        min_prominence_c = max(x_heights_c) / fact_prominence
                        min_disance_c = fact_distance  # len(yc)/fact_distance
                        peaks_c = find_peaks(yc, height=min_height_c, distance=min_disance_c,
                                             prominence=min_prominence_c, width=min_width)
                    except:
                        continue
                    # Peaks identification on discharge
                    yd_rev = [i * -1 for i in yd]
                    peaks_d = find_peaks(yd_rev, height=0.01)  # prominence=100
                    try:
                        x_pos_d = [xd[i] for i in peaks_d[0]]
                        x_heights_d = peaks_d[1]["peak_heights"]
                        min_height_d = max(x_heights_d) / fact_height
                        min_prominence_d = max(x_heights_d) / fact_prominence
                        min_disance_d = fact_distance  # len(yc)/fact_distance
                        peaks_d = find_peaks(yd_rev, height=min_height_d, distance=min_disance_d,
                                             prominence=min_prominence_d, width=min_width)
                    except:
                        continue
                    # Saving in dictionaries
                    peaks_der_c1[name_dict_x] = [xc[i] for i in peaks_c[0]]
                    peaks_der_c1[name_dict_y] = [yc[i] for i in peaks_c[0]]
                    peaks_der_c1_w[name_dict_y] = [i for i in peaks_c[1]['widths']]
                    peaks_der_d1[name_dict_x] = [xd[i] for i in peaks_d[0]]
                    peaks_der_d1[name_dict_y] = [yd[i] for i in peaks_d[0]]
                    peaks_der_d1_w[name_dict_y] = [i for i in peaks_d[1]['widths']]
                    if j == 0 or j == cycles_to_count - 1:
                        col_index = counter / number_of_cycles
                        xc, xd = dVdQ_dict_c[name_dict_x], dVdQ_dict_d[name_dict_x]
                        yc, yd = dVdQ_dict_c[name_dict_y], dVdQ_dict_d[name_dict_y]
                        lab = str(C_rates_list[i]) + "_" + str(counter)
                        if plots == 0:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                        if plots == 1:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(col_index), label=lab)
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(col_index))
                            plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o',
                                     color=plt.cm.cool(col_index), markersize=7)
                            plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o',
                                     color=plt.cm.cool(col_index), markersize=7)
                    counter += 1

        plt.legend(fontsize=16)  # bbox_to_anchor=(0.986,1)
        plt.xlim([min_up_Q, max_up_Q])
        # plt.ylim([-3,3])

        plot_save = directory_res + "\\" + titlef
        plt.savefig(plot_save)

# Cycle life test
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    peaks_derCL1_c = {}
    peaks_derCL1_c_w = {}
    peaks_derCL1_d = {}
    peaks_derCL1_d_w = {}
    fact_height = 1000  # 1000
    fact_prominence = 20  # 20
    fact_distance = 10  # 20
    min_width = 2
    for plots in range(2):
        if plots == 0:
            titlef = 'dVdQ_curves'
        if plots == 1:
            titlef = 'dVdQ_curves + peaks'
        fig = plt.figure(figsize=(13, 11))
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        plt.xlabel("Q (%)", fontsize=28)
        plt.ylabel("dV/dQ (V mAh$^{-1}$)", fontsize=28)
        plt.title(titlef, fontsize=24)
        for i in range(len(C_rates_lct)):
            range_cycles = len(num_cycles_dict[C_rates_lct[i]])
            every_plot = 5
            if range_cycles >= 50:
                every_plot = 10
            if range_cycles >= 100:
                every_plot = 20
            if range_cycles >= 200:
                every_plot = 30
            if range_cycles >= 300:
                every_plot = 50
            if range_cycles >= 1000:
                every_plot = 100
            for j in range(range_cycles):
                c = j / (range_cycles)
                name_x = str(C_rates_lct[i]) + "_Q_" + str(j)
                name_y = str(C_rates_lct[i]) + "_dVdQ_" + str(j)
                xc, xd = dVdQ_CL_c[name_x], dVdQ_CL_d[name_x]
                yc, yd = dVdQ_CL_c[name_y], dVdQ_CL_d[name_y]
                if plots == 0:
                    if j == 0:
                        plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                    for n in range(100):
                        if j + 1 == every_plot * n:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                    if j == range_cycles - 1:
                        plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                if plots == 1:
                    # Peaks identification on charge
                    peaks_c = find_peaks(yc, height=0.01)  # prominence=100
                    # print(peaks_c[0])
                    try:
                        x_pos_c = [xc[i] for i in peaks_c[0]]
                        x_heights_c = peaks_c[1]["peak_heights"]
                        min_height_c = max(x_heights_c) / fact_height
                        min_prominence_c = max(x_heights_c) / fact_prominence
                        min_disance_c = fact_distance  # len(yc)/fact_distance
                        peaks_c = find_peaks(yc, height=min_height_c, distance=min_disance_c,
                                             prominence=min_prominence_c, width=min_width)
                    except:
                        continue
                    # Peaks identification on discharge
                    yd_rev = [i * -1 for i in yd]
                    peaks_d = find_peaks(yd_rev, height=0.01)  # prominence=100
                    try:
                        x_pos_d = [xd[i] for i in peaks_d[0]]
                        x_heights_d = peaks_d[1]["peak_heights"]
                        min_height_d = max(x_heights_d) / fact_height
                        min_prominence_d = max(x_heights_d) / fact_prominence
                        min_disance_d = fact_distance  # len(yc)/fact_distance
                        peaks_d = find_peaks(yd_rev, height=min_height_d, distance=min_disance_d,
                                             prominence=min_prominence_d, width=min_width)
                    except:
                        continue
                    # Saving in dictionaries
                    peaks_derCL1_c[name_x] = [xc[i] for i in peaks_c[0]]
                    peaks_derCL1_c[name_y] = [yc[i] for i in peaks_c[0]]
                    peaks_derCL1_c_w[name_y] = [i for i in peaks_c[1]['widths']]
                    peaks_derCL1_d[name_x] = [xd[i] for i in peaks_d[0]]
                    peaks_derCL1_d[name_y] = [yd[i] for i in peaks_d[0]]
                    peaks_derCL1_d_w[name_y] = [i for i in peaks_d[1]['widths']]
                    # Plotting
                    if j == 0:
                        plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                        plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o', color=plt.cm.cool(c),
                                 markersize=7)
                        plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o', color=plt.cm.cool(c),
                                 markersize=7)
                    for n in range(100):
                        if j + 1 == every_plot * n:
                            plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                            plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                            plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o',
                                     color=plt.cm.cool(c), markersize=7)
                            plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o',
                                     color=plt.cm.cool(c), markersize=7)
                    if j == range_cycles - 1:
                        plt.plot(xc, yc, linewidth=2, color=plt.cm.cool(c), label=str(j + 1))
                        plt.plot(xd, yd, linewidth=2, color=plt.cm.cool(c))
                        plt.plot([xc[i] for i in peaks_c[0]], [yc[i] for i in peaks_c[0]], 'o', color=plt.cm.cool(c),
                                 markersize=7)
                        plt.plot([xd[i] for i in peaks_d[0]], [yd[i] for i in peaks_d[0]], 'o', color=plt.cm.cool(c),
                                 markersize=7)

            # plt.ylim([min(mins_E)+safe_remove, max(maxs_E)-safe_remove])
            plt.xlim([min_up_Q, max_up_Q])
            # plt.ylim([-10,10])
            plt.legend(fontsize=16)
            plot_save = directory_res_CL + "\\" + titlef
            plt.savefig(plot_save)

# Plotting capacity as a function of the cycle number
# Rate capability test
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    lines = {}
    for k in x_to_plot:
        if len(list_n_cycles) > 10:
            fig = plt.figure(figsize=(len(list_n_cycles), 9))
        else:
            fig = plt.figure(figsize=(11, 9))
        title = "End of (dis)charge capacity as a function of the cycle number - "
        titlef = title + str(k)
        plt.title(titlef)
        plt.scatter(list_n_cycles, ends_c[k], s=150, color='red', label="charge", marker="^")
        plt.scatter(list_n_cycles, ends_d[k], s=150, color='blue', label="discharge", marker="o")
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        plt.xlabel("Cycle number", fontsize=28)  # fontweight="bold"
        if k == "Q":
            plt.ylabel(k + " (mAh)", fontsize=28)  # fontweight="bold"
        if k == "Qn":
            plt.ylabel(k + " (mAh g$^{-1}$)", fontsize=28)  # fontweight="bold"
        name_y_lim_up_d = str(k) + "_up_discharge"
        name_y_lim_up_c = str(k) + "_up_charge"
        y_max_d = max(ends_d[k])
        y_max_c = max(ends_c[k])
        y_max = max(y_max_d, y_max_c)
        y_max = y_max + (0.05 * y_max)
        plt.xlim([0, len(list_n_cycles) + 0.5])
        plt.ylim([0, y_max])
        x_line_i = 0
        if formation == 1:
            y_line = [-10, y_max + 0.1 * y_max]
            x_line_i = cycles_per_formations + 0.5
            x_line = [x_line_i, x_line_i]
            x_text = cycles_per_formations / 4
            y_text = y_max - 0.1 * y_max
            plt.text(x_text, y_text, str(C_rates_list[0]), style='italic', fontsize=20, fontweight='semibold')
            plt.plot(x_line, y_line, color='black', linewidth=2, linestyle='dashed')
            for i in range(1, len(C_rates_list)):
                y_line = [-10, y_max + 0.1 * y_max]
                x_line_i = x_line_i + cycles_per_Crates
                x_line = [x_line_i, x_line_i]
                x_text = x_line_i - (cycles_per_Crates / 2) - 0.5
                y_text = y_max - 0.1 * y_max
                plt.text(x_text, y_text, str(C_rates_list[i]), style='italic', fontsize=20, fontweight='semibold')
                if i < len(C_rates_list) - 1:
                    plt.plot(x_line, y_line, color='black', linewidth=2, linestyle='dashed')
            plt.legend(fontsize=22, loc=3)
            plot_save = directory_res + "\\" + titlef
            plt.savefig(plot_save)

        if formation == 0:
            x_line_i = 0.5
            for i in range(len(C_rates_list)):
                y_line = [-10, y_max + 0.1 * y_max]
                x_line_i = x_line_i + cycles_per_Crates
                x_line = [x_line_i, x_line_i]
                x_text = x_line_i - (cycles_per_Crates / 2) - 0.5
                y_text = y_max - 0.1 * y_max
                plt.text(x_text, y_text, str(C_rates_list[i]), style='italic', fontsize=20, fontweight='semibold')
                if i < len(C_rates_list) - 1:
                    plt.plot(x_line, y_line, color='black', linewidth=2, linestyle='dashed')
            plt.legend(fontsize=22, loc=3)
            plot_save = directory_res + "\\" + titlef
            plt.savefig(plot_save)

# Cycle life test(s)
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    # Calculating Coulombic efficiency
    Qeff = {}
    for i in range(len(C_rates_lct)):
        Qeff[C_rates_lct[i]] = []
        name_Q = str(C_rates_lct[i]) + "_Q"
        for q in range(len(ends_CL_d[name_Q])):
            Qeff_i = (ends_CL_d[name_Q][q] / ends_CL_c[name_Q][q]) * 100
            Qeff[C_rates_lct[i]].append(Qeff_i)

    for i in range(len(C_rates_lct)):
        for k in x_to_plot:
            lenght_plot = len(num_cycles_dict[C_rates_lct[i]]) / 10
            # print(lenght_plot)
            if lenght_plot > 10:
                try:
                    fig = plt.figure(figsize=(len(lenght_plot), 9))
                except:
                    fig = plt.figure(figsize=(lenght_plot, 9))
            else:
                fig = plt.figure(figsize=(11, 9))

            ax1 = fig.add_subplot(1, 1, 1)

            title = "End of (dis)charge capacity as a function of the cycle number - "
            titlef = title + str(k)
            plt.title(titlef, fontsize=20)
            name_end = str(C_rates_lct[i]) + "_" + str(k)
            plt.scatter(num_cycles_dict[C_rates_lct[i]], ends_CL_c[name_end], s=150, color='red', label="charge",
                        marker="^")
            plt.scatter(num_cycles_dict[C_rates_lct[i]], ends_CL_d[name_end], s=150, color='blue', label="discharge",
                        marker="o")
            plt.xticks(fontsize=28)
            plt.yticks(fontsize=28)
            plt.xlabel("Cycle number", fontsize=36)  # fontweight="bold"
            if k == "Q":
                ax1.set_ylabel(k + " (mAh)", fontsize=36)  # fontweight="bold"
            if k == "Qn":
                ax1.set_ylabel(k + " (mAh g$^{-1}$)", fontsize=36)  # fontweight="bold"
            plt.xlim([0, len(num_cycles_dict[C_rates_lct[i]]) + 5])
            if k == "Q":
                plt.ylim([0, max(maxs_Q)])
            if k == "Qn":
                plt.ylim([0, max(maxs_Qn)])
            # Line 80%
            line_80_y = ends_CL_d[name_end][1] * 0.8
            line_80_y_l = [line_80_y, line_80_y]
            line_80_x_l = [-0.5, len(num_cycles_dict[C_rates_lct[i]]) + 1]
            label_line = "80% initial " + str(k) + "$_{(d)}$"
            plt.plot(line_80_x_l, line_80_y_l, color='darkturquoise', linewidth=5, linestyle='dashed', label=label_line)
            num_cycle_80 = 0
            for cap in range(len(ends_CL_d[name_end])):
                if ends_CL_d[name_end][cap] <= line_80_y:
                    num_cycle_80 = cap + 1
                    break
            text_80_y = line_80_y_l[0] + 0.05 * line_80_y_l[0]
            text_80_x = len(num_cycles_dict[C_rates_lct[i]]) - 0.1 * len(num_cycles_dict[C_rates_lct[i]])
            try:
                if num_cycle_80 > 0:
                    plt.text(text_80_x, text_80_y, str(num_cycle_80), color='darkturquoise', style='italic',
                             fontsize=24, fontweight='semibold')
            except:
                pass

            plt.legend(fontsize=24, loc=3)
            ax2 = ax1.twinx()
            ax2.set_ylabel('Coulombic efficiency (%)', color='purple', fontsize=32, rotation=90)
            ax2.scatter(num_cycles_dict[C_rates_lct[i]], Qeff[C_rates_lct[i]], s=150, color='purple',
                        label="Coulombic efficiency", marker="D")
            # ax2.set_ylim(0,110)
            plt.yticks(fontsize=28)
            plot_save = directory_res_CL + "\\" + titlef
            plt.savefig(plot_save)

if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    if len(CV_Qd) > 0 or len(CV_Qc) > 0:
        for i in range(len(C_rates_lct)):
            lenght_plot = len(num_cycles_dict[C_rates_lct[i]]) / 10
            # print(lenght_plot)
            if lenght_plot > 10:
                try:
                    fig = plt.figure(figsize=(len(lenght_plot), 9))
                except:
                    fig = plt.figure(figsize=(lenght_plot, 9))

            title = "Percentage of the end of (dis)charge capacity due to CV step as a function of the cycle number"
            plt.title(title, fontsize=20)
            if len(CV_Qc) > 0:
                plt.scatter(num_cycles_dict[C_rates_lct[i]], CV_Qc, s=150, color='red', label="charge", marker="^")
            if len(CV_Qd) > 0:
                plt.scatter(num_cycles_dict[C_rates_lct[i]], CV_Qd, s=150, color='blue', label="discharge", marker="o")

            plt.xlabel("Cycle number", fontsize=36)
            plt.ylabel("Q contribution of the CV step (%)", fontsize=32)
            plt.xticks(fontsize=28)
            plt.yticks(fontsize=28)
            plt.xlim([0, len(num_cycles_dict[C_rates_lct[i]]) + 5])
            # plt.ylim([0, max()])
            plt.legend(fontsize=24)

            plot_save = directory_res_CL + "\\" + title
            plt.savefig(plot_save)

# Plotting average voltage as a function of the cycle number
# Rate capability test
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    lines = {}
    for l in y_to_plot:
        if len(list_n_cycles) > 10:
            fig = plt.figure(figsize=(len(list_n_cycles), 9))
        else:
            fig = plt.figure(figsize=(11, 9))
        title = "End of (dis)charge average voltage as a function of the cycle number - "
        titlef = title + str(l)
        plt.title(titlef)
        plt.errorbar(list_n_cycles, ave_V_charge, yerr=std_V_charge, color='red', label="charge", capsize=3.5,
                     capthick=1.5, marker="^")
        plt.errorbar(list_n_cycles, ave_V_discharge, yerr=std_V_discharge, color='blue', label="discharge", capsize=3.5,
                     capthick=1.5, marker="o")
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        plt.xlabel("Cycle number", fontsize=28)  # fontweight="bold"
        plt.ylabel(l + " (V)", fontsize=28)  # fontweight="bold"
        plt.xlim([0, len(list_n_cycles) + 0.5])
        ymin_c = limits_plotting["E_down_charge"]
        ymin_d = limits_plotting["E_down_discharge"]
        ymin = min(ymin_d, ymin_c)[0]
        ymax_c = limits_plotting["E_up_charge"]
        ymax_d = limits_plotting["E_up_discharge"]
        ymax = max(ymax_c, ymax_d)[0]
        plt.ylim([ymin, ymax])
        x_line_i = 0
        if formation == 1:
            y_line = [-10, ymax + 0.1 * ymax]
            x_line_i = cycles_per_formations + 0.5
            x_line = [x_line_i, x_line_i]
            x_text = cycles_per_formations / 4
            y_text = ymax - 0.1 * ymax
            plt.text(x_text, y_text, str(C_rates_list[0]), style='italic', fontsize=20, fontweight='semibold')
            plt.plot(x_line, y_line, color='black', linewidth=2, linestyle='dashed')
            for i in range(1, len(C_rates_list)):
                y_line = [-10, ymax + 0.1 * ymax]
                x_line_i = x_line_i + cycles_per_Crates
                x_line = [x_line_i, x_line_i]
                x_text = x_line_i - (cycles_per_Crates / 2) - 0.5
                y_text = ymax - 0.1 * ymax
                plt.text(x_text, y_text, str(C_rates_list[i]), style='italic', fontsize=20, fontweight='semibold')
                if i < len(C_rates_list) - 1:
                    plt.plot(x_line, y_line, color='black', linewidth=2, linestyle='dashed')
            plt.legend(fontsize=22, loc=3)
            plot_save = directory_res + "\\" + titlef
            plt.savefig(plot_save)

        if formation == 0:
            x_line_i = 0.5
            for i in range(len(C_rates_list)):
                y_line = [-10, ymax + 0.1 * ymax]
                x_line_i = x_line_i + cycles_per_Crates
                x_line = [x_line_i, x_line_i]
                x_text = x_line_i - (cycles_per_Crates / 2) - 0.5
                y_text = ymax - 0.1 * ymax
                plt.text(x_text, y_text, str(C_rates_list[i]), style='italic', fontsize=20, fontweight='semibold')
                if i < len(C_rates_list) - 1:
                    plt.plot(x_line, y_line, color='black', linewidth=2, linestyle='dashed')
            plt.legend(fontsize=22, loc=3)
            plot_save = directory_res + "\\" + titlef
            plt.savefig(plot_save)

# Cycle life test(s) - not error bars
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for i in range(len(C_rates_lct)):
        for l in y_to_plot:
            lenght_plot = len(num_cycles_dict[C_rates_lct[i]]) / 10
            # print(lenght_plot)
            if lenght_plot > 10:
                try:
                    fig = plt.figure(figsize=(len(lenght_plot), 9))
                except:
                    fig = plt.figure(figsize=(lenght_plot, 9))
            else:
                fig = plt.figure(figsize=(11, 9))
            title = "End of (dis)charge average voltage as a function of the cycle number - "
            titlef = title + str(l)
            plt.title(titlef, fontsize=20)
            plt.scatter(num_cycles_dict[C_rates_lct[i]], ave_V_CL_c, s=150, color='red', label="charge", marker="^")
            plt.scatter(num_cycles_dict[C_rates_lct[i]], ave_V_CL_d, s=150, color='blue', label="discharge", marker="o")
            plt.xticks(fontsize=28)
            plt.yticks(fontsize=28)
            plt.xlabel("Cycle number", fontsize=36)  # fontweight="bold"

            plt.ylabel(l + " (V)", fontsize=36)  # fontweight="bold"
            plt.xlim([0, len(num_cycles_dict[C_rates_lct[i]]) + 0.5])
            # min_E = min_E - 0.1*min_E
            # max_E = max_E + 0.1*max_E
            plt.ylim([min_E, max_E])

            plt.legend(fontsize=24, loc=3)
            plot_save = directory_res_CL + "\\" + titlef
            plt.savefig(plot_save)

# Cycle life test(s) - error bars
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for i in range(len(C_rates_lct)):
        for l in y_to_plot:
            lenght_plot = len(num_cycles_dict[C_rates_lct[i]]) / 10
            # print(lenght_plot)
            if lenght_plot > 10:
                try:
                    fig = plt.figure(figsize=(len(lenght_plot), 9))
                except:
                    fig = plt.figure(figsize=(lenght_plot, 9))
            else:
                fig = plt.figure(figsize=(11, 9))
            title = "End of (dis)charge average voltage + std as a function of the cycle number - "
            titlef = title + str(l)
            plt.title(titlef, fontsize=20)
            plt.errorbar(num_cycles_dict[C_rates_lct[i]], ave_V_CL_c, yerr=std_V_CL_c, color='red', label="charge",
                         capsize=3.5, capthick=1.5, marker="^")
            plt.errorbar(num_cycles_dict[C_rates_lct[i]], ave_V_CL_d, yerr=std_V_CL_d, color='blue', label="discharge",
                         capsize=3.5, capthick=1.5, marker="o")
            plt.xticks(fontsize=28)
            plt.yticks(fontsize=28)
            plt.xlabel("Cycle number", fontsize=36)  # fontweight="bold"

            plt.ylabel(l + " (V)", fontsize=36)  # fontweight="bold"
            plt.xlim([0, len(num_cycles_dict[C_rates_lct[i]]) + 0.5])
            # min_E = min_E - 0.1*min_E
            # max_E = max_E + 0.1*max_E
            plt.ylim([min_E, max_E])

            plt.legend(fontsize=24, loc=3)
            plot_save = directory_res_CL + "\\" + titlef
            plt.savefig(plot_save)

# Plotting capacity variance as a function of the cycle number
# Cycle life test(s) - not error bars
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for d in range(2):
        if d == 0:
            dictionary = interpolated_CL_discharges_common
            var_listc = var_Q_discharge_CL1
            var_listd = var_Q_charge_CL1
            c_ref = 0
        if d == 1:
            dictionary = interpolated_CL_charges_common
            var_listc = var_Q_discharge_CL2
            var_listd = var_Q_charge_CL2
            c_ref = 9

        for i in range(len(C_rates_lct)):
            lenght_plot = len(num_cycles_dict[C_rates_lct[i]]) / 10
            # print(lenght_plot)
            if lenght_plot > 10:
                try:
                    fig = plt.figure(figsize=(len(lenght_plot), 9))
                except:
                    fig = plt.figure(figsize=(lenght_plot, 9))
            else:
                fig = plt.figure(figsize=(11, 9))
            title = "(dis)charge capacity variance with respect to cycle number " + str(c_ref + 1)
            # titlef = title + str(l)
            plt.title(title, fontsize=20)
            plt.scatter(num_cycles_dict[C_rates_lct[i]][c_ref:], var_listc, s=150, color='red', label="charge",
                        marker="^")
            plt.scatter(num_cycles_dict[C_rates_lct[i]][c_ref:], var_listd, s=150, color='blue', label="discharge",
                        marker="o")
            plt.xticks(fontsize=28)
            plt.yticks(fontsize=28)
            plt.xlabel("Cycle number", fontsize=36)  # fontweight="bold"

            plt.ylabel(u'Δ' + "Q (mAh$^{2}$)", fontsize=36)  # fontweight="bold"
            plt.xlim([0, len(num_cycles_dict[C_rates_lct[i]]) + 0.5])
            # min_E = min_E - 0.1*min_E
            # max_E = max_E + 0.1*max_E
            # plt.ylim([min_E, max_E])

            plt.legend(fontsize=24, loc=3)
            plot_save = directory_res_CL + "\\" + title
            plt.savefig(plot_save)

# Plotting average capacity (+- std) for each C-rate
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    if len(C_rates_list) > 1:
        for k in x_to_plot:
            size_x = (len(C_rates_list) + 1) * 2
            fig = plt.figure(figsize=(size_x, 9))
            title = "Average (dis)charge capacity as a function of the C-rate - "
            titlef = title + str(k)
            plt.title(titlef)
            plt.errorbar(C_rates_list, aves_end_c[k], yerr=stds_end_c[k], linewidth=2, color='red', label="charge",
                         capsize=3.5, capthick=1.5, linestyle='dashed', marker="o")
            plt.errorbar(C_rates_list, aves_end_d[k], yerr=stds_end_d[k], linewidth=2, color='blue', label="discharge",
                         capsize=3.5, capthick=1.5, linestyle='dashed', marker="o")
            plt.xticks(fontsize=22)
            plt.yticks(fontsize=22)
            plt.xlabel("C-rate", fontsize=28)
            if k == "Q":
                plt.ylabel(k + " (mAh)", fontsize=28)  # fontweight="bold"
            if k == "Qn":
                plt.ylabel(k + " (mAh g$^{-1}$)", fontsize=28)  # fontweight="bold"
            y_max_d = max(aves_end_d[k])
            y_max_c = max(aves_end_c[k])
            y_max = max(y_max_d, y_max_c)
            y_max = y_max + (0.05 * y_max)
            plt.ylim([0, y_max])
            plt.legend(fontsize=22)

            plot_save = directory_res + "\\" + titlef
            plt.savefig(plot_save)

# Predicting cycle life
# Beginning from cycle 1
if pred_cycles > 0:
    import warnings

    warnings.filterwarnings("ignore")

    directory_res_CL_pred = directory_res_CL + "\\Future"

    try:
        os.makedirs(directory_res_CL_pred)
    except:
        pass


    def normaliser(x):
        """Create data normalisation function"""
        x_mean, x_std = np.mean(x), np.std(x)
        x_scaled = (x - x_mean) / x_std
        return x_scaled, x_mean, x_std
        # return x / x_max, x_mean, x_std


    def deNorm(x_scaled, x_mean, x_sd):
        x_original = (x_scaled * x_sd) + x_mean
        return x_original


    pkl_filename = "Cycles_predictor.pkl"
    with open(pkl_filename, 'rb') as file:
        best_estimator = pickle.load(file)

    cycles_per_predictions = 30

    c_d = [ends_CL_c, ends_CL_d]

    predictions_print = {}

    for i in range(len(C_rates_lct)):
        for k in x_to_plot:
            for c in c_d:
                fig = plt.figure(figsize=(13, 11))
                plt.xticks(fontsize=22)
                plt.yticks(fontsize=22)
                if c == ends_CL_c:
                    title = "Predicted_charge_" + str(k)
                if c == ends_CL_d:
                    title = "Predicted_discharge_" + str(k)
                dict_b = {}
                dict_e = {}
                name_end = str(C_rates_lct[i]) + "_" + str(k)
                plt.scatter(num_cycles_dict[C_rates_lct[i]], c[name_end], s=45, color='#E95A1A', label="Measured",
                            marker="o")
                name_print_c = title + "_measured_cycles"
                name_print_Q = title + "_measured_Q"
                predictions_print[name_print_c] = num_cycles_dict[C_rates_lct[i]]
                predictions_print[name_print_Q] = c[name_end]
                beg_list = c[name_end][:cycles_per_predictions + 1]
                del beg_list[0]
                beginning, mean_b, std_b = normaliser(beg_list)
                end, mean_e, std_e = normaliser(c[name_end][-cycles_per_predictions:])
                for b in range(len(beginning)):
                    c_name = "Cycle_" + str(b + 1)
                    dict_b[c_name] = []
                    dict_b[c_name].append(beginning[b])
                    dict_e[c_name] = []
                    dict_e[c_name].append(end[b])

                max_pred = pred_cycles + num_cycles_dict[C_rates_lct[i]][-1] - cycles_per_predictions
                Q_pred_n = []
                extra_cycles = []
                for delta in range(1, max_pred + 1):
                    dict_b["Delta future"] = []
                    dict_b["Delta future"].append(delta)
                    extra = cycles_per_predictions + delta + 1
                    extra_cycles.append(extra)
                    pred = pd.DataFrame(dict_b)
                    pred_b = best_estimator.predict(pred)
                    Q_pred_n.append(float(pred_b[0]))
                Q_pred = deNorm(np.array(Q_pred_n), mean_b, std_b)
                plt.scatter(extra_cycles, Q_pred, s=45, color='#6ED2E6', label="Predicted (beginning)", marker="^")
                name_print_c = title + "_pred_beg_cycles"
                name_print_Q = title + "_pred_beg_Q"
                predictions_print[name_print_c] = extra_cycles
                predictions_print[name_print_Q] = Q_pred.tolist()

                Q_pred_n = []
                extra_cycles = []
                for delta in range(1, pred_cycles + 1):
                    dict_e["Delta future"] = []
                    dict_e["Delta future"].append(delta)
                    extra = num_cycles_dict[C_rates_lct[i]][-1] + delta
                    extra_cycles.append(extra)
                    pred = pd.DataFrame(dict_e)
                    pred_e = best_estimator.predict(pred)
                    Q_pred_n.append(float(pred_e[0]))
                Q_pred = deNorm(np.array(Q_pred_n), mean_e, std_e)
                plt.scatter(extra_cycles, Q_pred, s=45, color='#4467B3', label="Predicted (end)", marker="o")
                name_print_c = title + "_pred_end_cycles"
                name_print_Q = title + "_pred_end_Q"
                predictions_print[name_print_c] = extra_cycles
                predictions_print[name_print_Q] = Q_pred.tolist()

                x_line = [cycles_per_predictions + 1, cycles_per_predictions + 1]
                y_line = [-0.5, 1000]
                plt.plot(x_line, y_line, color='#6ED2E6', linewidth=2, linestyle='dashed')  # label="Inputs/outputs"
                x_line = [1, 1]
                y_line = [-0.5, 1000]
                plt.plot(x_line, y_line, color='#6ED2E6', linewidth=2, linestyle='dashed')  # label="Inputs/outputs"
                end_l_limit = num_cycles_dict[C_rates_lct[i]][-1] - cycles_per_predictions
                end_u_limit = num_cycles_dict[C_rates_lct[i]][-1]
                x_line = [end_l_limit, end_l_limit]
                y_line = [-0.5, 1000]
                plt.plot(x_line, y_line, color='#4467B3', linewidth=2, linestyle='dashed')
                x_line = [end_u_limit, end_u_limit]
                y_line = [-0.5, 1000]
                plt.plot(x_line, y_line, color='#4467B3', linewidth=2, linestyle='dashed')

                plt.xlim([-1, len(num_cycles_dict[C_rates_lct[i]]) + pred_cycles + 0.5])
                if k == "Q":
                    plt.ylim([0, max(maxs_Q)])
                if k == "Qn":
                    plt.ylim([0, max(maxs_Qn)])

                plt.xlabel("Cycle number", fontsize=36)
                if k == "Q":
                    plt.ylabel(k + " (mAh)", fontsize=36)  # fontweight="bold"
                if k == "Qn":
                    plt.ylabel(k + " (mAh g$^{-1}$)", fontsize=36)  # fontweight="bold"

                plt.title(title)
                plt.legend(fontsize=28)
                plot_save = directory_res_CL_pred + "\\" + title
                plt.savefig(plot_save)

                # Save predictions in csv
    max_len = 0
    for key in predictions_print:
        if len(predictions_print[key]) > max_len:
            max_len = len(predictions_print[key])

    for key in predictions_print:
        lenght = len(predictions_print[key])
        if lenght < max_len:
            for i in range(lenght, max_len):
                predictions_print[key].append("")

    predictions_print_df = pd.DataFrame(predictions_print)

    path_save = directory_res_CL_pred + '\\' + "Predictions.csv"
    predictions_print_df.to_csv(path_save, index=False)

    for key in predictions_print:
        count_remove = 0
        for i in range(len(predictions_print[key])):
            if predictions_print[key][i] == "":
                try:
                    del predictions_print[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

# Report all results in DataFames and print

# Averages
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    path_save = directory_res + '\\' + "full_averaged_charge.csv"
    ave_results_charge.to_csv(path_save, index=False)
    path_save = directory_res + '\\' + "full_averaged_discharge.csv"
    ave_results_discharge.to_csv(path_save, index=False)

# Raw (dis)charge data
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    discharges_print = {}
    charges_print = {}
    discharges_results = pd.DataFrame()
    charges_results = pd.DataFrame()
    for key in discharges:
        splitted_key = key.split("list_discharge_")
        if "C_rate" not in key:
            discharges_print[splitted_key[1]] = discharges[key]

    for key in charges:
        if "C_rate" not in key:
            splitted_key = key.split("list_charge_")
            charges_print[splitted_key[1]] = charges[key]

    max_lenght = max([len(discharges_print[key]) for key in discharges_print])
    for key in discharges_print:
        lenght = len(discharges_print[key])
        for i in range(lenght, max_lenght):
            discharges_print[key].append("")

    max_lenght = max([len(charges_print[key]) for key in charges_print])
    for key in charges_print:
        lenght = len(charges_print[key])
        for i in range(lenght, max_lenght):
            charges_print[key].append("")

    for key in discharges_results:
        print(len(discharges_results[key]))
    discharges_results = pd.DataFrame(discharges_print)
    charges_results = pd.DataFrame(charges_print)

    path_save = directory_res + '\\' + "raw_discharges.csv"
    discharges_results.to_csv(path_save, index=False)

    path_save = directory_res + '\\' + "raw_charges.csv"
    charges_results.to_csv(path_save, index=False)

# Raw (dis)charge data - CL
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    max_lenght = max([len(cycle_life_discharges[key]) for key in cycle_life_discharges])
    for key in cycle_life_discharges:
        lenght = len(cycle_life_discharges[key])
        for i in range(lenght, max_lenght):
            cycle_life_discharges[key].append("")

    max_lenght = max([len(cycle_life_charges[key]) for key in cycle_life_charges])
    for key in cycle_life_charges:
        lenght = len(cycle_life_charges[key])
        for i in range(lenght, max_lenght):
            cycle_life_charges[key].append("")

    cycle_life_discharges_print = pd.DataFrame(cycle_life_discharges)
    cycle_life_charges_print = pd.DataFrame(cycle_life_charges)

    path_save = directory_res_CL + '\\' + "raw_discharges.csv"
    cycle_life_discharges_print.to_csv(path_save, index=False)

    path_save = directory_res_CL + '\\' + "raw_charges.csv"
    cycle_life_charges_print.to_csv(path_save, index=False)

# dQ/dV
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    dQdV_dict_d_results = pd.DataFrame()
    dQdV_dict_c_results = pd.DataFrame()
    max_lenght = max([len(dQdV_dict_d[key]) for key in dQdV_dict_d])
    for key in dQdV_dict_d:
        lenght = len(dQdV_dict_d[key])
        for i in range(lenght, max_lenght):
            try:
                dQdV_dict_d[key].tolist().append("")
            except:
                dQdV_dict_d[key].append("")

    max_lenght = max([len(dQdV_dict_c[key]) for key in dQdV_dict_c])
    for key in dQdV_dict_c:
        lenght = len(dQdV_dict_c[key])
        for i in range(lenght, max_lenght):
            try:
                dQdV_dict_c[key].tolist().append("")
            except:
                dQdV_dict_c[key].append("")

    dQdV_dict_d_results = pd.DataFrame(dQdV_dict_d)
    dQdV_dict_c_results = pd.DataFrame(dQdV_dict_c)

    path_save = directory_res + '\\' + "dQ_dV_discharges.csv"
    dQdV_dict_d_results.to_csv(path_save, index=False)

    path_save = directory_res + '\\' + "dQ_dV_charges.csv"
    dQdV_dict_c_results.to_csv(path_save, index=False)

# dQ/dV - Cycle life
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    max_lenght = max([len(dQdv_CL_d[key]) for key in dQdv_CL_d])
    for key in dQdv_CL_d:
        lenght = len(dQdv_CL_d[key])
        for i in range(lenght, max_lenght):
            try:
                dQdv_CL_d[key].tolist().append("")
            except:
                dQdv_CL_d[key].append("")

    max_lenght = max([len(dQdv_CL_c[key]) for key in dQdv_CL_c])
    for key in dQdv_CL_c:
        lenght = len(dQdv_CL_c[key])
        for i in range(lenght, max_lenght):
            try:
                dQdv_CL_c[key].tolist().append("")
            except:
                dQdv_CL_c[key].append("")

    dQdV_dict_d_results = pd.DataFrame(dQdv_CL_d)
    dQdV_dict_c_results = pd.DataFrame(dQdv_CL_c)

    path_save = directory_res_CL + '\\' + "dQ_dV_discharges.csv"
    dQdV_dict_d_results.to_csv(path_save, index=False)

    path_save = directory_res_CL + '\\' + "dQ_dV_charges.csv"
    dQdV_dict_c_results.to_csv(path_save, index=False)

# Removing ""
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    for key in discharges:
        count_remove = 0
        for i in range(len(discharges[key])):
            if discharges[key][i] == "":
                try:
                    del discharges[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in charges:
        count_remove = 0
        for i in range(len(charges[key])):
            if charges[key][i] == "":
                try:
                    del charges[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in dQdV_dict_d:
        count_remove = 0
        for i in range(len(dQdV_dict_d[key])):
            if dQdV_dict_d[key][i] == "":
                try:
                    del dQdV_dict_d[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in dQdV_dict_c:
        count_remove = 0
        for i in range(len(dQdV_dict_c[key])):
            if dQdV_dict_c[key][i] == "":
                try:
                    del dQdV_dict_c[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

# Removing "" - CL
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for key in cycle_life_discharges:
        count_remove = 0
        for i in range(len(cycle_life_discharges[key])):
            if cycle_life_discharges[key][i] == "":
                try:
                    del cycle_life_discharges[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in cycle_life_charges:
        count_remove = 0
        for i in range(len(cycle_life_charges[key])):
            if cycle_life_charges[key][i] == "":
                try:
                    del cycle_life_charges[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in dQdv_CL_d:
        count_remove = 0
        for i in range(len(dQdv_CL_d[key])):
            if dQdv_CL_d[key][i] == "":
                try:
                    del dQdv_CL_d[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in dQdv_CL_c:
        count_remove = 0
        for i in range(len(dQdv_CL_c[key])):
            if dQdv_CL_c[key][i] == "":
                try:
                    del dQdv_CL_c[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

# dV/dQ
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    dVdQ_dict_d_results = pd.DataFrame()
    dVdQ_dict_c_results = pd.DataFrame()
    max_lenght = max([len(dVdQ_dict_d[key]) for key in dVdQ_dict_d])
    for key in dVdQ_dict_d:
        lenght = len(dVdQ_dict_d[key])
        for i in range(lenght, max_lenght):
            try:
                dVdQ_dict_d[key].tolist().append("")
            except:
                dVdQ_dict_d[key].append("")

    max_lenght = max([len(dVdQ_dict_c[key]) for key in dVdQ_dict_c])
    for key in dVdQ_dict_c:
        lenght = len(dVdQ_dict_c[key])
        for i in range(lenght, max_lenght):
            try:
                dVdQ_dict_c[key].tolist().append("")
            except:
                dVdQ_dict_c[key].append("")

    dVdQ_dict_d_results = pd.DataFrame(dVdQ_dict_d)
    dVdQ_dict_c_results = pd.DataFrame(dVdQ_dict_c)

    path_save = directory_res + '\\' + "dV_dQ_discharges.csv"
    dVdQ_dict_d_results.to_csv(path_save, index=False)

    path_save = directory_res + '\\' + "dV_dQ_charges.csv"
    dVdQ_dict_c_results.to_csv(path_save, index=False)

# dV/dQ - Cycle life
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    max_lenght = max([len(dVdQ_CL_d[key]) for key in dVdQ_CL_d])
    for key in dVdQ_CL_d:
        lenght = len(dVdQ_CL_d[key])
        for i in range(lenght, max_lenght):
            try:
                dVdQ_CL_d[key].tolist().append("")
            except:
                dVdQ_CL_d[key].append("")

    max_lenght = max([len(dVdQ_CL_c[key]) for key in dVdQ_CL_c])
    for key in dVdQ_CL_c:
        lenght = len(dVdQ_CL_c[key])
        for i in range(lenght, max_lenght):
            try:
                dVdQ_CL_c[key].tolist().append("")
            except:
                dVdQ_CL_c[key].append("")

    dVdQ_dict_d_results = pd.DataFrame(dVdQ_CL_d)
    dVdQ_dict_c_results = pd.DataFrame(dVdQ_CL_c)

    path_save = directory_res_CL + '\\' + "dV_dQ_discharges.csv"
    dVdQ_dict_d_results.to_csv(path_save, index=False)

    path_save = directory_res_CL + '\\' + "dV_dQ_charges.csv"
    dVdQ_dict_c_results.to_csv(path_save, index=False)

# Removing ""
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    for key in discharges:
        count_remove = 0
        for i in range(len(discharges[key])):
            if discharges[key][i] == "":
                try:
                    del discharges[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in charges:
        count_remove = 0
        for i in range(len(charges[key])):
            if charges[key][i] == "":
                try:
                    del charges[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in dVdQ_dict_d:
        count_remove = 0
        for i in range(len(dVdQ_dict_d[key])):
            if dVdQ_dict_d[key][i] == "":
                try:
                    del dVdQ_dict_d[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in dVdQ_dict_c:
        count_remove = 0
        for i in range(len(dVdQ_dict_c[key])):
            if dVdQ_dict_c[key][i] == "":
                try:
                    del dVdQ_dict_c[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

# Removing "" - CL
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for key in cycle_life_discharges:
        count_remove = 0
        for i in range(len(cycle_life_discharges[key])):
            if cycle_life_discharges[key][i] == "":
                try:
                    del cycle_life_discharges[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in cycle_life_charges:
        count_remove = 0
        for i in range(len(cycle_life_charges[key])):
            if cycle_life_charges[key][i] == "":
                try:
                    del cycle_life_charges[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in dVdQ_CL_d:
        count_remove = 0
        for i in range(len(dVdQ_CL_d[key])):
            if dVdQ_CL_d[key][i] == "":
                try:
                    del dVdQ_CL_d[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

    for key in dVdQ_CL_c:
        count_remove = 0
        for i in range(len(dVdQ_CL_c[key])):
            if dVdQ_CL_c[key][i] == "":
                try:
                    del dVdQ_CL_c[key][count_remove:]
                    break
                except:
                    pass
            count_remove += 1

        # End of (dis)charges
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    end_of_discharges_print = {}
    end_of_charges_print = {}
    end_of_discharges_results = pd.DataFrame()
    end_of_charges_results = pd.DataFrame()
    end_of_discharges_print["Cycle number"] = list_n_cycles
    end_of_charges_print["Cycle number"] = list_n_cycles
    for k in x_to_plot:
        end_of_discharges_print[k] = ends_d[k]
        end_of_charges_print[k] = ends_c[k]

    end_of_discharges_results = pd.DataFrame(end_of_discharges_print)
    end_of_charges_results = pd.DataFrame(end_of_charges_print)
    path_save = directory_res + '\\' + "end_of_discharges_capacity.csv"
    end_of_discharges_results.to_csv(path_save, index=False)

    path_save = directory_res + '\\' + "end_of_charges_capacity.csv"
    end_of_charges_results.to_csv(path_save, index=False)

# End of (dis)charges - CL
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for i in range(len(C_rates_lct)):
        end_of_discharges_CL_print = pd.DataFrame()
        end_of_charges_CL_print = pd.DataFrame()
        end_of_discharges_CL_print["Cycle number"] = num_cycles_dict[C_rates_lct[i]]
        end_of_charges_CL_print["Cycle number"] = num_cycles_dict[C_rates_lct[i]]
        for k in x_to_plot:
            key = str(C_rates_lct[i]) + "_" + str(k)
            end_of_discharges_CL_print[k] = ends_CL_d[key]
            end_of_charges_CL_print[k] = ends_CL_c[key]

        path_save = directory_res_CL + '\\' + "end_of_discharges_capacity" + str(C_rates_lct[i]) + ".csv"
        end_of_discharges_CL_print.to_csv(path_save, index=False)

        path_save = directory_res_CL + '\\' + "end_of_charges_capacity" + str(C_rates_lct[i]) + ".csv"
        end_of_charges_CL_print.to_csv(path_save, index=False)

        if len(CV_Qc) > 0 or len(CV_Qd) > 0:
            CV_effect = pd.DataFrame()
            CV_effect["Cycle number"] = num_cycles_dict[C_rates_lct[i]]
            if len(CV_Qc) > 0:
                CV_effect["Qc from CV step (%)"] = CV_Qc
            if len(CV_Qd) > 0:
                CV_effect["Qd from CV step (%)"] = CV_Qd
            path_save = directory_res_CL + '\\' + "Q_from_CV_" + str(C_rates_lct[i]) + ".csv"
            CV_effect.to_csv(path_save, index=False)

# Coulombic efficiency - CL
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for i in range(len(C_rates_lct)):
        coul_eff = pd.DataFrame()
        coul_eff["Cycle number"] = num_cycles_dict[C_rates_lct[i]]
        coul_eff["Coulombic efficiency (%)"] = Qeff[C_rates_lct[i]]

        path_save = directory_res_CL + '\\' + "Coulombic_efficiency_" + str(C_rates_lct[i]) + ".csv"
        coul_eff.to_csv(path_save, index=False)

# Averaged end of (dis)charge
if len(C_rates_list) > 1 and cycles_per_Crates > 0 or len(C_rates_list) > 1 and cycles_per_formations > 0:
    ave_end_of_discharges_print = {}
    ave_end_of_charges_print = {}
    ave_end_of_discharges_results = pd.DataFrame()
    ave_end_of_charges_results = pd.DataFrame()
    ave_end_of_discharges_print["C_rates"] = C_rates_list
    ave_end_of_charges_print["C_rates"] = C_rates_list
    for k in x_to_plot:
        name_std = "std_" + k
        ave_end_of_discharges_print[k] = aves_end_d[k]
        ave_end_of_discharges_print[name_std] = stds_end_d[k]
        ave_end_of_charges_print[k] = aves_end_c[k]
        ave_end_of_charges_print[name_std] = stds_end_c[k]

    ave_end_of_discharges_results = pd.DataFrame(ave_end_of_discharges_print)
    ave_end_of_charges_results = pd.DataFrame(ave_end_of_charges_print)
    path_save = directory_res + '\\' + "ave_end_discharges_capacity.csv"
    ave_end_of_discharges_results.to_csv(path_save, index=False)

    path_save = directory_res + '\\' + "ave_end_charges_capacity.csv"
    ave_end_of_charges_results.to_csv(path_save, index=False)

# Average_voltage (dis)charges
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    average_voltage_RC_charge = pd.DataFrame()
    average_voltage_RC_discharge = pd.DataFrame()
    average_voltage_RC_charge["Cycle number"] = list_n_cycles
    average_voltage_RC_charge["Average voltage"] = ave_V_charge
    average_voltage_RC_charge["Std voltage"] = std_V_charge
    average_voltage_RC_discharge["Cycle number"] = list_n_cycles
    average_voltage_RC_discharge["Average voltage"] = ave_V_discharge
    average_voltage_RC_discharge["Std voltage"] = std_V_discharge

    path_save = directory_res + '\\' + "Average_voltage_charge.csv"
    average_voltage_RC_charge.to_csv(path_save, index=False)

    path_save = directory_res + '\\' + "Average_voltage_discharge.csv"
    average_voltage_RC_discharge.to_csv(path_save, index=False)

# Average_voltage (dis)charges - CL
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for i in range(len(C_rates_lct)):
        average_voltage_RC_charge_CL = pd.DataFrame()
        average_voltage_RC_discharge_CL = pd.DataFrame()
        average_voltage_RC_charge_CL["Cycle number"] = num_cycles_dict[C_rates_lct[i]]
        average_voltage_RC_charge_CL["Average voltage"] = ave_V_CL_c
        average_voltage_RC_charge_CL["Std voltage"] = std_V_CL_c
        average_voltage_RC_discharge_CL["Cycle number"] = num_cycles_dict[C_rates_lct[i]]
        average_voltage_RC_discharge_CL["Average voltage"] = ave_V_CL_d
        average_voltage_RC_discharge_CL["Std voltage"] = std_V_CL_d

        path_save = directory_res_CL + '\\' + "Average_voltage_charge.csv"
        average_voltage_RC_charge_CL.to_csv(path_save, index=False)

        path_save = directory_res_CL + '\\' + "Average_voltage_discharge.csv"
        average_voltage_RC_discharge_CL.to_csv(path_save, index=False)

# Q(V) variance (dis)charges - CL
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    for i in range(len(C_rates_lct)):
        for d in range(2):
            if d == 0:
                var_list1 = var_Q_discharge_CL1
                var_list2 = var_Q_discharge_CL2
                name_file = "Variance_capacity_discharge.csv"
            if d == 1:
                var_list1 = var_Q_charge_CL1
                var_list2 = var_Q_charge_CL2
                name_file = "Variance_capacity_charge.csv"
            c_ref1 = 0
            c_ref2 = 9
            varianceQ_CL = pd.DataFrame()
            name1 = "Cycle number " + str(c_ref1)
            varianceQ_CL[name1] = num_cycles_dict[C_rates_lct[i]][c_ref1:]
            varianceQ_CL["Variance Q(V) - 1"] = var_list1
            name2 = "Cycle number " + str(c_ref2)
            cycles_adjusted = num_cycles_dict[C_rates_lct[i]][c_ref2:]
            if len(cycles_adjusted) < len(num_cycles_dict[C_rates_lct[i]][c_ref1:]):
                for j in range(c_ref2):
                    cycles_adjusted.append("")
            varianceQ_CL[name2] = cycles_adjusted
            variance_adjusted = var_list2
            if len(variance_adjusted) < len(num_cycles_dict[C_rates_lct[i]][c_ref1:]):
                for j in range(c_ref2):
                    variance_adjusted.append("")
            varianceQ_CL["Variance Q(V) - 2"] = var_list2

            path_save = directory_res_CL + '\\' + name_file
            varianceQ_CL.to_csv(path_save, index=False)

# Inputs.log file
logs = []
if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
    path_log_0 = directory_res + '\\' + "Inputs.log"
    log_0 = open(path_log_0, 'w')
    logs.append(log_0)
if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
    path_log_1 = directory_res_CL + '\\' + "Inputs.log"
    log_1 = open(path_log_1, 'w')
    logs.append(log_1)

for log in logs:
    log.write("Software used and file location:" + '\n')
    log.write("Software = " + software + '\n')
    log.write("Folder path  = " + folder_path + '\n')
    log.write("Filename = " + filename + '\n')
    log.write('\n')

    log.write("Rate capability test:" + '\n')
    log.write("C_rate_list = ")
    for i in range(len(C_rates_list)):
        log.write(str(C_rates_list[i]) + ", ")
    log.write('\n')
    log.write("Cycle(s) per formation = " + str(cycles_per_formations) + '\n')
    log.write("Cycle(s) per C_rate = " + str(cycles_per_Crates) + '\n')
    log.write('\n')

    log.write("Cycle life test:" + '\n')
    log.write("C_rate = ")
    for i in range(len(C_rates_lct)):
        log.write(str(C_rates_lct[i]) + ", ")
    log.write('\n')
    log.write("Cycles per cycle life test = ")
    for i in range(len(Cycles_lct)):
        log.write(str(Cycles_lct[i]) + ", ")
    log.write('\n')
    log.write("Cycle(s) to be predicted = " + str(pred_cycles))
    log.write('\n')
    log.write('\n')

    log.write("Smoothing:" + '\n')
    log.write("Points per linear interpolation = " + str(points) + '\n')
    log.write("Savgol filter - windows length = " + str(sv_windows_lenght) + '\n')
    log.write("Savgol filter - polynomial order = " + str(sv_pol_ord) + '\n')
    log.write('\n')

    log.write("Normalization:" + '\n')
    log.write("Normalization mass (if any) = " + str(mass_active_mg) + " (mg)" + '\n')
    log.write('\n')

    log.write("ELab:" + '\n')
    log.write("Experiment name = " + str(Elab_exp) + '\n')
    log.write("User name = " + str(User_name) + '\n')
    log.write("Token = " + str(token_ELab) + '\n')
    log.close()

    try:
        if log == log_0:
            with open(path_log_0, 'r') as f:
                inputs = f.readlines()
                f.close()
    except:
        continue
    try:
        if log == log_1:
            with open(path_log_1, 'r') as f:
                inputs = f.readlines()
                f.close()
    except:
        continue

# Picks' characteristics - dQdV
try:
    if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
        max_peaks = 0
        for key in peaks_der_c:
            if len(peaks_der_c[key]) > max_peaks:
                max_peaks = len(peaks_der_c[key])
            if len(peaks_der_d[key]) > max_peaks:
                max_peaks = len(peaks_der_d[key])

        for d_c_direction in range(2):
            if d_c_direction == 0:
                path_peaks = directory_res + '\\' + 'Identified_Peaks_dQdV_charge.csv'
            if d_c_direction == 1:
                path_peaks = directory_res + '\\' + 'Identified_Peaks_dQdV_discharge.csv'
            header = []
            header.append("Cycle")
            for i in range(max_peaks):
                pos = "peak_pos_" + str(i + 1)
                header.append(pos)
                hieight = "peak_height_" + str(i + 1)
                header.append(hieight)
                width = "peak_width_" + str(i + 1)
                header.append(width)

            with open(path_peaks, 'w') as f:
                write = csv.writer(f)
                write.writerow(header)
                for i in range(len(C_rates_list)):
                    if formation == 1 and C_rates_list[i] == C_rates_list[0]:
                        cycles_to_count = cycles_per_formations
                    elif finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                        cycles_to_count = cycle_per_last
                    else:
                        cycles_to_count = cycles_per_Crates
                    for j in range(cycles_to_count):
                        list_to_write = []
                        name_dict_y = str(C_rates_list[i]) + "_dQdV_" + str(j + 1)
                        name_dict_x = str(C_rates_list[i]) + "_E_" + str(j + 1)
                        name_cond = str(C_rates_list[i]) + " - cycle: " + str(j)
                        list_to_write.append(name_cond)
                        if d_c_direction == 0:
                            try:
                                num_peaks_i = len(peaks_der_c[name_dict_x])
                            except:
                                num_peaks_i = 0
                            for p in range(num_peaks_i):
                                list_to_write.append(peaks_der_c[name_dict_x][p])
                                list_to_write.append(peaks_der_c[name_dict_y][p])
                                list_to_write.append(peaks_der_c_w[name_dict_y][p])
                            write.writerow(list_to_write)
                        if d_c_direction == 1:
                            try:
                                num_peaks_i = len(peaks_der_d[name_dict_x])
                            except:
                                num_peaks_i = 0
                            for p in range(num_peaks_i):
                                list_to_write.append(peaks_der_d[name_dict_x][p])
                                list_to_write.append(peaks_der_d[name_dict_y][p])
                                list_to_write.append(peaks_der_d_w[name_dict_y][p])
                            write.writerow(list_to_write)
except:
    pass

# Picks' characteristics - CL
try:
    if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
        max_peaks = 0
        for key in peaks_derCL_c:
            if len(peaks_derCL_c[key]) > max_peaks:
                max_peaks = len(peaks_derCL_c[key])
            if len(peaks_derCL_d[key]) > max_peaks:
                max_peaks = len(peaks_derCL_d[key])

        for d_c_direction in range(2):
            if d_c_direction == 0:
                path_peaks = directory_res_CL + '\\' + 'Identified_Peaks_dQdV_charge.csv'
            if d_c_direction == 1:
                path_peaks = directory_res_CL + '\\' + 'Identified_Peaks_dQdV_discharge.csv'
            header = []
            header.append("Cycle")
            for i in range(max_peaks):
                pos = "peak_pos_" + str(i + 1)
                header.append(pos)
                hieight = "peak_height_" + str(i + 1)
                header.append(hieight)
                width = "peak_width_" + str(i + 1)
                header.append(width)

            with open(path_peaks, 'w') as f:
                write = csv.writer(f)
                write.writerow(header)
                for i in range(len(C_rates_lct)):
                    range_cycles = len(num_cycles_dict[C_rates_lct[i]])
                    for j in range(range_cycles):
                        list_to_write = []
                        name_dict_y = str(C_rates_lct[i]) + "_dQdV_" + str(j)
                        name_dict_x = str(C_rates_lct[i]) + "_E_" + str(j)
                        name_cond = str(C_rates_lct[i]) + " - cycle: " + str(j)
                        list_to_write.append(name_cond)
                        if d_c_direction == 0:
                            try:
                                num_peaks_i = len(peaks_derCL_c[name_dict_x])
                            except:
                                num_peaks_i = 0
                            for p in range(num_peaks_i):
                                list_to_write.append(peaks_derCL_c[name_dict_x][p])
                                list_to_write.append(peaks_derCL_c[name_dict_y][p])
                                list_to_write.append(peaks_derCL_c_w[name_dict_y][p])
                            write.writerow(list_to_write)
                        if d_c_direction == 1:
                            try:
                                num_peaks_i = len(peaks_derCL_d[name_dict_x])
                            except:
                                num_peaks_i = 0
                            for p in range(num_peaks_i):
                                list_to_write.append(peaks_derCL_d[name_dict_x][p])
                                list_to_write.append(peaks_derCL_d[name_dict_y][p])
                                list_to_write.append(peaks_derCL_d_w[name_dict_y][p])
                            write.writerow(list_to_write)
except:
    pass

# Picks' characteristics - dVdQ
try:
    if len(C_rates_list) > 0 and cycles_per_Crates > 0 or len(C_rates_list) > 0 and cycles_per_formations > 0:
        max_peaks = 0
        for key in peaks_der_c1:
            if len(peaks_der_c1[key]) > max_peaks:
                max_peaks = len(peaks_der_c1[key])
            if len(peaks_der_d1[key]) > max_peaks:
                max_peaks = len(peaks_der_d1[key])

        for d_c_direction in range(2):
            try:
                if d_c_direction == 0:
                    path_peaks = directory_res + '\\' + 'Identified_Peaks_dVdQ_charge.csv'
                if d_c_direction == 1:
                    path_peaks = directory_res + '\\' + 'Identified_Peaks_dVdQ_discharge.csv'
                header = []
                header.append("Cycle")
                for i in range(max_peaks):
                    pos = "peak_pos_" + str(i + 1)
                    header.append(pos)
                    hieight = "peak_height_" + str(i + 1)
                    header.append(hieight)
                    width = "peak_width_" + str(i + 1)
                    header.append(width)

                with open(path_peaks, 'w') as f:
                    write = csv.writer(f)
                    write.writerow(header)
                    for i in range(len(C_rates_list)):
                        if formation == 1 and C_rates_list[i] == C_rates_list[0]:
                            cycles_to_count = cycles_per_formations
                        elif finish_check == 0 and C_rates_list[i] == C_rates_list[-1]:
                            cycles_to_count = cycle_per_last
                        else:
                            cycles_to_count = cycles_per_Crates
                        for j in range(cycles_to_count):
                            list_to_write = []
                            name_dict_y = str(C_rates_list[i]) + "_dVdQ_" + str(j + 1)
                            name_dict_x = str(C_rates_list[i]) + "_Q_" + str(j + 1)
                            name_cond = str(C_rates_list[i]) + " - cycle: " + str(j)
                            list_to_write.append(name_cond)
                            if d_c_direction == 0:
                                try:
                                    num_peaks_i = len(peaks_der_c1[name_dict_x])
                                except:
                                    num_peaks_i = 0
                                for p in range(num_peaks_i):
                                    list_to_write.append(peaks_der_c1[name_dict_x][p])
                                    list_to_write.append(peaks_der_c1[name_dict_y][p])
                                    list_to_write.append(peaks_der_c1_w[name_dict_y][p])
                                write.writerow(list_to_write)
                            if d_c_direction == 1:
                                try:
                                    num_peaks_i = len(peaks_der_d1[name_dict_x])
                                except:
                                    num_peaks_i
                                for p in range(num_peaks_i):
                                    list_to_write.append(peaks_der_d1[name_dict_x][p])
                                    list_to_write.append(peaks_der_d1[name_dict_y][p])
                                    list_to_write.append(peaks_der_d1_w[name_dict_y][p])
                                write.writerow(list_to_write)

            except:
                continue
except:
    pass

# Picks' characteristics - CL
try:
    if len(C_rates_lct) > 0 and len(Cycles_lct) > 0:
        max_peaks = 0
        for key in peaks_derCL1_c:
            if len(peaks_derCL1_c[key]) > max_peaks:
                max_peaks = len(peaks_derCL1_c[key])
            if len(peaks_derCL1_d[key]) > max_peaks:
                max_peaks = len(peaks_derCL1_d[key])

        for d_c_direction in range(2):
            if d_c_direction == 0:
                path_peaks = directory_res_CL + '\\' + 'Identified_Peaks_dVdQ_charge.csv'
            if d_c_direction == 1:
                path_peaks = directory_res_CL + '\\' + 'Identified_Peaks_dVdQ_discharge.csv'
            header = []
            header.append("Cycle")
            for i in range(max_peaks):
                pos = "peak_pos_" + str(i + 1)
                header.append(pos)
                hieight = "peak_height_" + str(i + 1)
                header.append(hieight)
                width = "peak_width_" + str(i + 1)
                header.append(width)

            with open(path_peaks, 'w') as f:
                write = csv.writer(f)
                write.writerow(header)
                for i in range(len(C_rates_lct)):
                    range_cycles = len(num_cycles_dict[C_rates_lct[i]])
                    for j in range(range_cycles):
                        list_to_write = []
                        name_dict_y = str(C_rates_lct[i]) + "_dVdQ_" + str(j)
                        name_dict_x = str(C_rates_lct[i]) + "_Q_" + str(j)
                        name_cond = str(C_rates_lct[i]) + " - cycle: " + str(j)
                        list_to_write.append(name_cond)
                        if d_c_direction == 0:
                            try:
                                num_peaks_i = len(peaks_derCL1_c[name_dict_x])
                            except:
                                num_peaks_i = 0
                            for p in range(num_peaks_i):
                                list_to_write.append(peaks_derCL1_c[name_dict_x][p])
                                list_to_write.append(peaks_derCL1_c[name_dict_y][p])
                                list_to_write.append(peaks_derCL1_c_w[name_dict_y][p])
                            write.writerow(list_to_write)
                        if d_c_direction == 1:
                            try:
                                num_peaks_i = len(peaks_derCL1_d[name_dict_x])
                            except:
                                num_peaks_i = 0
                            for p in range(num_peaks_i):
                                list_to_write.append(peaks_derCL1_d[name_dict_x][p])
                                list_to_write.append(peaks_derCL1_d[name_dict_y][p])
                                list_to_write.append(peaks_derCL1_d_w[name_dict_y][p])
                            write.writerow(list_to_write)
except:
    continue



# Send data to experiment in ELab
if len(token_ELab) > 3:
    import elabapy
    import json
    import zipfile


    def zip_directory(folder_path, zip_path):
        with zipfile.ZipFile(zip_path, mode='w') as zipf:
            len_dir_path = len(folder_path)
            for root, _, files in os.walk(folder_path):
                for file in files:
                    file_path = os.path.join(root, file)
                    zipf.write(file_path, file_path[len_dir_path:])


    manager = elabapy.Manager(endpoint="https://fb08-eln.phys.chemie.uni-giessen.de/api/v1/", token=token_ELab,
                              verify=False)

    params_get_exp = {'limit': 100000}

    all_exp = manager.get_all_experiments(params_get_exp)
    all_exp_df = pd.DataFrame.from_records(all_exp)

    for i in range(len(all_exp_df["title"])):
        if all_exp_df["title"][i] == Elab_exp:
            exp_id = all_exp_df["id"][i]

    folders_to_zip = [directory_res, directory_res_CL]

    for f in folders_to_zip:
        path_for_zip = f.split("\\Results_")[0] + "\\Results_" + f.split("\\Results_")[1].split(".txt")[0] + ".zip"
        zip_directory(f, path_for_zip)
        with open(path_for_zip, 'r+b') as myfold:
            params = {'file': myfold}
            manager.upload_to_experiment(exp_id, params)
        os.remove(path_for_zip)

else:
    print("No Token entered")

# Messaging
from tkinter import *
root = Tk()
MyLable = Label(root, text="Data recovered, analysed, and plotted correctly. Enjoy!")
MyLable.pack()
root.mainloop()

##################################################################################### END OF THE PROGRAM #####################################################################################