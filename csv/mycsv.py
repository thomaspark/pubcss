from os.path import join, split
from glob import glob
import pandas as pd
import numpy as np

def insert_csv():

    dz_path = './dz.csv'
    rx_path = './rx.csv'
    ry_path = './ry.csv'

    csv_files=glob('./disg/*.csv')

    # 抽出する節点番号
    apoint = []
    for i in range(0, 10101, 1010):
        for j in range(0, 101, 10):
            apoint.append(i+j)

    # 荷重を載荷した節点番号
    lpoint = []
    for i in range(0, 5151, 101):
        for j in range(51):
            lpoint.append(i+j)

    # 結果を格納する変数
    dz = pd.DataFrame(np.zeros((len(lpoint),len(apoint))), index=lpoint, columns=apoint) # z軸方向の変位量
    rx = pd.DataFrame(np.zeros((len(lpoint),len(apoint))), index=lpoint, columns=apoint) # y軸周りの回転角
    ry = pd.DataFrame(np.zeros((len(lpoint),len(apoint))), index=lpoint, columns=apoint) # y軸周りの回転角

    # すべての csvファイル
    for file in csv_files:

        r = split(file)[1] # テーブル名作成
        r = int(r[0:-4]) # テーブル名から拡張子を削除

        l = pd.read_csv(file, index_col=0)

        cdz = l[2:3] # z軸方向の変位量
        crx = l[3:4] # x軸周りの回転角
        cry = l[4:5] # y軸周りの回転角

        for c in apoint:
            col = str(c)
            dz.at[r, c] = cdz[col]
            rx.at[r, c] = crx[col]
            ry.at[r, c] = cry[col]

    dz.to_csv(dz_path)
    rx.to_csv(rx_path)
    ry.to_csv(ry_path)

if __name__ == '__main__':
    insert_csv()