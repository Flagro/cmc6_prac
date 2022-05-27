from unittest import result
import pandas as pd


with open("./output.txt") as file:
    lines = file.readlines()

lines = [line.replace("\n", "") for line in lines]
lines = [line.split() for line in lines]

result_dict_list = []
for line in lines:
    result_dict_list.append({
        "N_x": line[0],
        "N_y": line[1],
        "N_z": line[2],
        "Test Id": line[3],
        "Threads Cnt": line[4],
        "Expected Threads Cnt": line[5],
        "CG Max Iterations": line[6],
        "CG Epsilon": line[7],
        "Total Time": line[8],
        "Dot Product Time": line[9],
        "SpMV Time": line[10],
        "axpby Time": line[11],
        "Dot Product BW": line[12],
        "SpMV BW": line[13],
        "axpby BW": line[14],
        "CG Iterations Used": line[15],
        "CG Residual": line[16],
        "Residual": line[17],
        "Error": line[18],
        "Mode": line[19].replace("_", " "),
        "Threading Mode": str(line[4]) + " (" + str(line[19].replace("_", " ")) + ")"
    })

result_df = pd.DataFrame(result_dict_list)
result_df = result_df.groupby([
    "N_x", "N_y", "N_z", "Test Id", "Threads Cnt", "Expected Threads Cnt", "Mode", "Threading Mode"])[[
        "CG Max Iterations", "CG Epsilon", "Total Time", "Dot Product Time", "SpMV Time", "axpby Time", 
        "Dot Product BW", "SpMV BW", "axpby BW", "CG Iterations Used", "CG Residual", "Residual", "Error"]].mean().reset_index()
result_df = result_df.loc[result_df.groupby(["N_x", "N_y", "N_z", "Threads Cnt"])["Total Time"].idxmin()].reset_index(drop=True)
del result_df["Error"]
del result_df["CG Max Iterations"]
del result_df["CG Epsilon"]
del result_df["Test Id"]
del result_df["Expected Threads Cnt"]
del result_df["Dot Product BW"]
del result_df["SpMV BW"]
del result_df["axpby BW"]
del result_df["CG Residual"]
del result_df["Mode"]
result_df["N_x"] = result_df["N_x"].astype(int)
result_df["N_y"] = result_df["N_y"].astype(int)
result_df["N_z"] = result_df["N_z"].astype(int)
result_df["N"] = result_df["N_x"] * result_df["N_y"] * result_df["N_z"]
del result_df["N_x"]
del result_df["N_y"]
del result_df["N_z"]
print(result_df)
result_df.to_csv("./results.csv")
