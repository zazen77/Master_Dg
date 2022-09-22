from openpyxl import load_workbook
wb = load_workbook("nodecopy.xlsx")
ws = wb.active



for x in range(1, 66):
    ws.move_range("B1:B40", rows=40*x, cols=-1)
    ws.delete_cols(2)

wb.save("sample_korean.xlsx")



