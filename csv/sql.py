# -*- coding: utf-8 -*-
import sqlite3
from os.path import join, split
from glob import glob
import pandas as pd

def insert_csv():

    db_path = './disg.db'

    csv_files=glob('./disg/*.csv')

    for csv in csv_files:

        table_name = split(csv)[1] # テーブル名作成

        table_name = table_name[0:-4] # テーブル名から拡張子を削除

        df = pd.read_csv(csv, dtype=object) # CSV 読込

        with sqlite3.connect(db_path) as conn:
            df.T.to_sql(table_name, con=conn) # SQLiteにCSVをインポート

if __name__ == '__main__':
    insert_csv()