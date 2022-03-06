import json
from test_data import testData
from calcrate import FrameCalc
from dataInput import dataInput
import pandas as pd

js = testData()

# 入力データを読み込む
inp = dataInput()

error = inp.setJSON(js)
if error != None:
    pass

for i in range(0, 5151, 101):
    for j in range(51):
        n = i+j
        lo = inp.load["1"]
        load_node = lo["load_node"]
        load = load_node[0]
        load["n"] = str(n)
        calc = FrameCalc(inp)

        # 計算開始
        error, result = calc.calcrate()
        if error != None:
            pass
        
        res = result["1"]
        disg = res["disg"]
        reac = res["reac"]
        fsec = res["fsec"]

        df_disg = pd.DataFrame.from_dict(disg)
        df_reac = pd.DataFrame.from_dict(reac)
        df_fsec = pd.DataFrame.from_dict(fsec)

        df_disg.to_csv("./csv/disg/{}.csv".format(n))
        df_reac.to_csv("./csv/reac/{}.csv".format(n))
        df_fsec.to_csv("./csv/fsec/{}.csv".format(n))

