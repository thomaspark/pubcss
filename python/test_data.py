# -*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import axes3d    
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def testData():
  span = [3, 4, 3]
  pich = 0.1

  min = 0
  max = sum(span) + pich

  gpoint = np.arange(min, max, pich)
  gpoint

  """# 入力データの作成"""

  body = dict()

  """## 格点 を定義する"""

  line = []
  z = 0
  for x in gpoint:
    for y in gpoint:
      line.append([x, y, z])

  df = pd.DataFrame(line, columns=["x", "y", "z"])
  body["node"] = df.T.to_dict('dict')


  """## 部材 を定義する

  骨組構造における部材の番号と下記の情報を定義する

  *   ni: 部材の左端の点
  *   nj: 部材の右端の点
  *   e: 後述する材料番号

  """

  # i行当りの分割数
  n = int(sum(span) / pich)

  # 縦列の節点番号
  index2, j = [0], 0
  for i in span:
    j += int(i / pich)
    index2.append(j)

  # 横列の節点番号
  index1 = []
  for i in range(len(index2)):
    j = index2[i] * n + index2[i]
    index1.append(j)

  member = dict()
  key = 0

  # 縦桁
  for i in index1:
    for ni in range(i, i+n):
      member[str(key)] = {
          "ni":ni, 
          "nj":ni+1, 
          "e":"1"
      }
      key += 1

  # 横桁
  for ni in index2:
    nj = ni+101
    for i in range(n):
      member[str(key)] = {
          "ni":ni, 
          "nj":nj, 
          "e":"1"
      }
      ni = nj
      nj = ni+101
      key += 1

  body["member"] = member


  """
  ## シェル要素を定義する

  3次元構造解析における4節点シェル要素の番号と下記の情報を定義する

  *    nodes: シェルを構成する点
  *    e: 後述する材料番号
  """

  shell = dict()

  key = 0
  for i in range(len(gpoint)-1):
    j = i * n + i
    for k in range(n):
      p1 = j + k
      p2 = p1 + 1
      p3 = p2 + n
      p4 = p3 + 1
      shell[str(key)] = {
          "nodes": [p1, p2, p4, p3],
          "e":"2"
      }
      key+=1


  body["shell"] = shell


  """## 材料の番号と下記の情報を定義する

  ### 棒要素の場合

  *  E: 弾性係数
  *  G: せん断弾性係数
  *  Xp: 膨張係数
  *  A: 断面積
  *  J: ねじり定数
  *  Iy: y軸周りの断面二次モーメント
  *  Iz: z軸周りの断面二次モーメント

  ### シェル要素の場合

  *  E: 弾性係数
  *  G: せん断弾性係数
  *  A: 板厚
  """

  element = {
      "1":{ "E":26500000, "G":11000000, "Xp":0.00001, "A":1.3, "J":0.2288, "Iy":0.183083, "Iz":0.108333 },
      "2":{ "E":26500000, "G":11000000, "A":0.3 }
    }

  body["element"]  = {
      "1": element
      }


  """## 支点の情報を定義する

  *   n: 固定する格点番号
  *   tx: x方向の変位を固定する場合は1, 固定しない場合は0
  *   ty: y方向の変位を固定する場合は1, 固定しない場合は0
  *   tz: z方向の変位を固定する場合は1, 固定しない場合は0
  *   rx: x軸周りの回転を固定する場合は1, 固定しない場合は0
  *   ry: y軸周りの回転を固定する場合は1, 固定しない場合は0
  *   rz: z軸周りの回転を固定する場合は1, 固定しない場合は0
  """

  fix_node = []

  for i in index1: # [0, 3030, 7070, 10100]
    for j in index2: # [0, 30, 70, 100]
      fix_node.append({
        "n":str(i + j),"tx":1,"ty":1,"tz":1,"rx":0,"ry":0,"rz":0
      })

  body["fix_node"] = { 
      "1": fix_node
      }


  """## 荷重の情報を定義する

  *   fix_node: 支点の情報の番号
  *   element: 材料情報の番号
  *   load_node
      *   n: 荷重を載荷する格点番号
      *   tx: x方向に作用する荷重
      *   ty: y方向に作用する荷重
      *   tz: z方向に作用する荷重
      *   rx: x軸周りに作用するモーメント
      *   ry: y軸周りに作用するモーメント
      *   rz: z軸周りに作用するモーメント
  """

  load = []

  load.append(
      { "n":"1","tx":0,"ty":0,"tz":100,"rx":0,"ry":0,"rz":0 }
  )

  body["load"] = {
        "1":{
          "fix_node": "1",
          "element": "1",
          "load_node": load
          }
      }

  return body

