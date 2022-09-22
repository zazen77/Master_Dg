import pandas as pd
import numpy as np
from openpyxl import load_workbook
import time

start = time.time()
wb1 = load_workbook("강제변위node.xlsx") # sample.xlsx 파일에서 wb 을 불러옴
ws1 = wb1.active # 활성화된 Sheet

#노드넘버 파일
Nodenum_xl = pd.read_excel('Node_number.xlsx', index_col=0)
Nodenum_list = Nodenum_xl.to_numpy()
Nodenum = Nodenum_list.reshape(1, Nodenum_list.size)



#좌표파일 열기
xycoord_xl = pd.read_excel('Volumeloss2xycoord.xlsx', index_col=0)
xycoord_list = xycoord_xl.to_numpy()

p = 0

for k in range(-40, -35, 1):

    for i in range(0, 1):
        # load set name

        row = i%40
        col = (i//40)*2

        test = np.array(['강제변위' + str(k) + "m", '강제변위-' + str(i), Nodenum[0, i], 'Total', 'None', "110000", xycoord_list[row, col], xycoord_list[row, col+1]])
        # ws1.append(test)
        # ws1["A" + str(i+3)] = '강제변위' + str(k) + "m"
        # # numbering
        # ws1["B" + str(i+3)] = '강제변위-' + str(i)
        # # node number
        # ws1["C" + str(i+3)] = Nodenum[0, i]
        #
        # ws1["D" + str(i+3)] = 'Total'
        # ws1["E" + str(i+3)] = 'None'
        # ws1["F" + str(i+3)] = "110000"

        # x and y coordinate

        # ws1["G" + str(i+3)] = xycoord_list[row, col]
        # ws1["H" + str(i+3)] = xycoord_list[row, col+1]
        # #
        file_name = "강제변위" + str(k) + "m.xlsx"
        T1 = np.concatenate([T1, T2], axis=0)
        wb1.save(file_name)
        print(i)

    p = p + 1

print("time :", time.time() - start)