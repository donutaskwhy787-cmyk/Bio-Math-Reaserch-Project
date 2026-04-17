from openpyxl import load_workbook

filePath = "twentyFourAndUp.xls"
wb = load_workbook(filePath)
sheet = wb.active

