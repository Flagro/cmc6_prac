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
pd.DataFrame(result_dict_list).to_csv("./results.csv")
