"""Adds the Column_definitions sheet to Supplementary Data 1, 3, 4, 5, 6, 7 (definitions in supp_data_definitions.py).
usage: python3 supp_data_add_column_definitions.py <folder with the workbooks from step 15> <output folder>"""
import openpyxl, sys, os
from openpyxl.styles import Font,Alignment,PatternFill
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from supp_data_definitions import C
IN, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
SRC = IN + '/Supplementary_Data_%d.xlsx'
for n in (1,3,4,5,6,7):
    wb=openpyxl.load_workbook(SRC%n); rows=[]
    for (k,sh),d in C.items():
        if k!=n: continue
        ws=wb[sh]
        hdr=[c.value for c in next(ws.iter_rows(min_row=1,max_row=1))]
        missing=[h for h in hdr if h not in d and h is not None]
        extra=[h for h in d if h not in hdr]
        if extra: print('EXTRA',n,sh,extra)
        for h in hdr:
            if h in d: rows.append((sh,h,d[h]))
            elif h is not None: rows.append((sh,h,'(column name is a full description)'))
        if missing: print('full-name columns',n,sh,len(missing))
    ws=wb.create_sheet('Column_definitions'); ws.append(['Sheet','Column','Definition'])
    for r in rows: ws.append(list(r))
    for c in ws[1]: c.font=Font(bold=True); c.fill=PatternFill('solid',fgColor='EAF0F6')
    ws.column_dimensions['A'].width=30; ws.column_dimensions['B'].width=36; ws.column_dimensions['C'].width=110
    for row in ws.iter_rows(min_row=2):
        for c in row: c.alignment=Alignment(wrap_text=True,vertical='top')
    ws.freeze_panes='A2'
    rd=wb['README']; rd.cell(rd.max_row+1,1,'Column_definitions: meaning and unit of every column in every sheet, including abbreviations (GE, RPKG, RPKM, pLDDT, TM-score, null_RPKG and others).')
    wb.save(f'{OUT}/Supplementary_Data_{n}.xlsx'); print(n,len(rows),'rows')
