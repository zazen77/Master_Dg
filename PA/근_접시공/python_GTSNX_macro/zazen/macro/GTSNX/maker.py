import pandas as pd
import numpy as np
from openpyxl import load_workbook
import time

start = time.time()

wb1 = load_workbook("강제변위node.xlsx") # sample.xlsx 파일에서 wb 을 불러옴
ws1 = wb1.active # 활성화된 Sheet

#노드넘버 파일
Nodenum_xl = pd.read_excel('Node_number.xlsx', index_col=0)
Nodenum_list = np.transpose(Nodenum_xl.to_numpy())
Nodenum = Nodenum_list.reshape(Nodenum_list.size, 1)



#좌표파일 열기
xycoord_xl = pd.read_excel('Volumeloss2xycoord.xlsx', index_col=0)
xycoord_list = xycoord_xl.to_numpy()

p = 0

for k in range(150, -51, -10):
    N = -k + 150
    for i in range(0, 4040):
        # load set name
        ws1["A" + str(i+3)] = '강제변위' + str(k) + "m"
        # numbering
        ws1["B" + str(i+3)] = '강제변위-' + str(i)
        # node number
        ws1["C" + str(i+3)] = Nodenum[i, 0]

        ws1["D" + str(i+3)] = 'Total'
        ws1["E" + str(i+3)] = 'None'
        ws1["F" + str(i+3)] = "011000"

        # x and y coordinate
        row = i%40
        col = (i//40)*2 + (N*2)
        ws1["H" + str(i+3)] = xycoord_list[row, col]
        ws1["I" + str(i+3)] = xycoord_list[row, col + 1]

        # print(i)
    file_name = "강제변위" + str(k) + "m.xlsx"
    wb1.save(file_name)
    print(k)
    p = p + 1

print("time :", time.time() - start)