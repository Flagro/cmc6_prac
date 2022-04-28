import pandas as pd


with open("./output.txt") as file:
    lines = file.readlines()

lines = [line.replace("\n", "") for line in lines]
lines = [line.split() for line in lines]

result_dict_list = []
for line in lines:
    result_dict_list.append({
        "N": line[0],
        "Test Id": line[1],
        "Threads Cnt": line[2],
        "Expected Threads Cnt": line[9],
        "Triangulation Time": line[3],
        "Backward Move Time": line[4],
        "Total Time": line[5],
        "Residual": line[6],
        "Error": line[7],
        "Mode": line[8].replace("_", " ")
    })
pd.DataFrame(result_dict_list).to_csv("./results.csv")
