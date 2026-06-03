import xlrd as xl
def get_label(label_name, filedata_root, label_sheet, station_sheet, label_col, station_col):
    label = []
    station_label = []
    wb = xl.open_workbook(filedata_root + label_name)
    sheet_1 = wb.sheet_by_index(label_sheet)
    sheet_2 = wb.sheet_by_index(station_sheet)
    label_data = sheet_1.col_values(label_col)
    station_data = sheet_2.col_values(label_col)
    label = list(map(int, label_data))
    station_label = list(map(int, station_data))
    return label, station_label