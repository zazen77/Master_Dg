from openpyxl import load_workbook # 파일 불러오기
import pyautogui

pyautogui.click(767, 769, duration=0.1) #Y정렬 버튼 누르기
pyautogui.click(822,889, duration=0.1) #테이블 버튼 누르기
pyautogui.sleep(2) #테이블 생성시간 슬립
pyautogui.click(437,219, duration=0.1) #절점1번 누르기
pyautogui.sleep(0.5)

pyautogui.hotkey("ctrl", "c")
pyautogui.click(531, 707, duration=0.1) #테이블 닫기

pyautogui.click(372, 1063, duration=0.1) #엑셀 창 누르기

pyautogui.click(-1767, 246, duration=0.5) # 엑셀 a2셀 누르기
pyautogui.sleep(0.1)


pyautogui.sleep(0.1)
pyautogui.hotkey("ctrl", "v")

pyautogui.click(-31, 1004, duration=0.1)#우측버튼 1번 누르기

pyautogui.click(1054, 1061, duration=0.1) # gtsnx 창 누르기
pyautogui.click(842, 203, duration=0.1) # 전체선택해제 누르기 


