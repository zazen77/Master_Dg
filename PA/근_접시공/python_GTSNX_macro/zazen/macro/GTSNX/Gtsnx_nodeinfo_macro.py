import pyautogui
import pyperclip
import numpy as np
import math
import pandas as pd

w = pyautogui.getWindowsWithTitle("GTS NX")[0]
print(w)
if not w.isActive:  # 현재 활성화가 되지 않았다면
    w.activate()  # 활성화 (맨 앞으로 가져오기)

if not w.isMaximized:  # 현재 최대화가 되지 않았다면
    w.maximize()  # 최대화

# 신규터널의 절점정보 정리하기 (절점 40개 원형이면 좌표 안변함)
point_1 = 959, 230
point_11 = 1301, 572

r1 = abs(point_11[0] - point_1[0])
r2 = abs(point_11[1] - point_1[1])

if abs(r1 - r2) > 3:
    pyautogui.alert("1번 좌표와 11번 좌표를 다시 설정해주세요.", "경고")  # 확인 버튼만 있는 팝업

# for i in range(40):
#     t = math.cos(2 * math.pi / 40* i)

# 절점 드래그 해서 정보 얻기 (좌측상단좌표 기준)


# 처음 절점 정보 0번

i = 0
square_x = r1 * math.sin(2 * math.pi / 40 * i)
square_y = r1 * math.cos(2 * math.pi / 40 * i)

# print(square_x)
# print(square_y)
square_Lconer = square_x + point_1[0] - 5, point_11[1] - square_y - 5

# print(square_1_Lconer,square_1_Rconer)
# pyautogui.PAUSE = 1
pyautogui.hotkey("ctrl", "shift", "c", duration=0.25)

pyautogui.moveTo(square_Lconer)
pyautogui.drag(10, 10)  # 너무 빠른 동작으로 drag 수행이 안될때는 duration 값 설정

pyautogui.click(201, 700)  # 오름차순 체크
pyautogui.click(83, 699)  # x 정렬
pyautogui.click(129, 671)  # 절점 정보들
pyautogui.hotkey("ctrl", "a")
pyautogui.hotkey("ctrl", "c")
pyautogui.click(244, 819)  # 닫기클릭

# line_1 = np.array(pyperclip.paste().split())
line_1 = pyperclip.paste().split()  # 공백으로 분리
line_node_1 = np.reshape(list(map(int, line_1)), (101, 1))  # map 함수로 정수로 바꾸고 리스트로 묶기

# 절점정보 반복문 1~39
for i in range(1, 40):
    square_x = r1 * math.sin(2 * math.pi / 40 * i)
    square_y = r1 * math.cos(2 * math.pi / 40 * i)

    # print(square_x)
    # print(square_y)
    square_Lconer = square_x + point_1[0] - 5, point_11[1] - square_y - 5

    # print(square_1_Lconer,square_1_Rconer)
    # pyautogui.PAUSE = 1
    pyautogui.hotkey("ctrl", "shift", "c", duration=0.25)

    pyautogui.moveTo(square_Lconer)
    pyautogui.drag(10, 10)  # 너무 빠른 동작으로 drag 수행이 안될때는 duration 값 설정

    pyautogui.click(201, 700)  # 오름차순 체크
    pyautogui.click(83, 699)  # x 정렬
    pyautogui.click(129, 671)  # 절점 정보들
    pyautogui.hotkey("ctrl", "a", )
    pyautogui.hotkey("ctrl", "c", )
    pyautogui.click(244, 819)  # 닫기클릭

    line_2 = pyperclip.paste().split()  # 공백으로 분리
    line_node_2 = np.reshape(list(map(int, line_2)), (101, 1))  # map 함수로 정수로 바꾸고 리스트로 묶기

    line_node_1 = np.concatenate((line_node_1, line_node_2), axis=1)

print(line_node_1)

# # 엑셀파일 쓰기
# wb = Workbook() # 새 워크북 생성
# ws = wb.active # 현재 활성화된 sheet 가져옴
# ws.title = "Node_number" # sheet 의 이름을 변경
# ws["A1"] = line_node_1
# wb.save("node_number.xlsx")
# wb.close()


df = pd.DataFrame(line_node_1).T
df.to_excel(excel_writer="D:/OneDrive - 한양대학교/조재은/학교수업 석사3기/paper/기존터널_근접시공/python/zazen/Node_number.xlsx")
