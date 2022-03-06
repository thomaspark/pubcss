
import numpy as np
import copy
from collections import OrderedDict

from dataManager import dataManager
from matrixBar import barMatrix
from matrixShell import shellMatrix


class gMatrix:

    def __init__(self, _inp: dataManager, _nodeIndex: dict):

        self.nfree = 6 # 支点の自由度

        self.inp = _inp

        # 支点条件を整理し解析用変数にセットする
        self.mpfix = self.set_fix(self.inp)

        # 要素剛性マトリックスの計算
        self.kmat = self.set_kMatrix(self.inp)

        # 全体剛性マトリックスの合計サイズを計測する
        self.Gindex, n = self.setMatrixIndex(self.inp, self.kmat, _nodeIndex)


        # 全体荷重マトリックス
        # xternal force vector
        self.fg = np.zeros(n, dtype=np.float64) # External force vector
        for ln in self.inp.load_node:
            lp = self.Gindex[ln['n']]  # 剛性マトリックスのインデックス（先頭）
            self.fg[lp+0] += ln['tx']  # load in x-direction
            self.fg[lp+1] += ln['ty']  # load in y-direction
            self.fg[lp+2] += ln['tz']  # load in z-direction
            self.fg[lp+3] += ln['rx']  # moment around x-axis
            self.fg[lp+4] += ln['ry']  # moment around y-axis
            self.fg[lp+5] += ln['rz']  # moment around z-axis


        # 全体剛性マトリックス
        # Global stiffness matrix
        self.gk = np.zeros([n, n], dtype=np.float64)
        for k in self.kmat:
            if k.label == 'Bar':
                self.set_barMatrix(k) # Bar部材の剛性マトリックス
            elif k.label == 'Shell':
                self.set_shellMatrix(k) # Shell部材の剛性マトリックス

        # 強制変位も考慮した荷重
        self.dg = np.zeros(n, dtype=np.float64) # External displacement vector
        for ln in self.inp.load_node:
            lp = self.Gindex[ln['n']]  # 剛性マトリックスのインデックス（先頭）
            self.dg[lp+0] += ln['dx']  # displacement in x-direction
            self.dg[lp+1] += ln['dy']  # displacement in y-direction
            self.dg[lp+2] += ln['dz']  # displacement in z-direction
            self.dg[lp+3] += ln['ax']  # rotate around x-axis
            self.dg[lp+4] += ln['ay']  # rotate around y-axis
            self.dg[lp+5] += ln['az']  # rotate around z-axis

        self.org_dp = copy.deepcopy(self.dg)

        fg = np.dot(self.gk, self.dg)
        self.fp  = self.fg - fg

        self.org_fp = copy.deepcopy(self.fp)

        # 強制変位が入力されている成分に対する処理
        for npos in range(len(self.dg)):
            if self.dg[npos] != 0:
                # 変位拘束が存在する行の成分を0にする
                self.gk[npos, :] = 0
                self.gk[:, npos] = 0

                # 対角成分を1にする
                self.gk[npos, npos] = 1

                # 強制変位成分に対応する荷重ベクトルを修正
                self.fp[npos] = self.dg[npos]


        # 変位が拘束されている成分に対する処理
        for ID in self.mpfix:
            fix = self.mpfix[ID]
            for i in range(len(fix)):
                if fix[i] == 0:
                    continue

                value1 = fix[i]
                npos = self.Gindex[ID] + i

                if value1 == 1:
                    # 変位拘束が存在する行の成分を0にする
                    self.gk[npos, :] = 0
                    self.gk[:, npos] = 0

                    # 対角成分を1にする
                    self.gk[npos, npos] = 1

                    # 変位拘束成分に対応する荷重ベクトルを修正
                    self.fp[npos] = 0.0

                else: # value != 0:
                    # バネ値がある場合対角成分にバネ剛性を加算する
                    self.gk[npos, npos] += value1


        # 対角成分がゼロとなる場合（ピン結合などの場合）1 とする
        for npos in range(len(self.gk)):
            if self.gk[npos, npos] == 0:
                self.gk[npos, :] = 0
                self.gk[:, npos] = 0
                self.gk[npos, npos] = 1
                self.fp[npos] = 0.0

        return


    # 要素剛性マトリックスの計算
    def set_kMatrix(self, inp: dataManager):
        kmat = []

        # Bar部材の剛性マトリックス
        for ID in inp.member:
            k = barMatrix(inp, ID)
            kmat.append(k)

        # Shell部材の剛性マトリックス
        for ID in inp.shell:
            k = shellMatrix(inp, ID)
            kmat.append(k)

        return kmat

    # Bar部材の剛性マトリックス を 全体剛性マトリックスに登録する
    def set_barMatrix(self, k: barMatrix):
        # ne = Array's Index - 1 = Member ID

        # 全体剛性マトリックスでの１インデックスを決定する
        index = []
        for ID in k.nodes:
            for i in range(k.nfree):
                index.append(self.Gindex[ID] + i)

        for i in range(len(index)):
            it = index[i]
            for j in range(len(index)):
                jt = index[j]
                self.gk[it, jt] += k.ck[i, j]

        for i in range(len(index)):
            it = index[i]
            self.fg[it] += k.wfe[i] + k.tfe[i]

        return

    # shell部材の剛性マトリックス を 全体剛性マトリックスに登録する
    def set_shellMatrix(self, k: shellMatrix):
        # ne = Array's Index - 1 = Shell ID

        index = []
        for ID in k.nodes:
            for i in range(k.nfree):
                index.append(self.Gindex[ID] + i)

        for i in range(len(index)):
            it = index[i]
            for j in range(len(index)):
                jt = index[j]
                self.gk[it, jt] += k.ck[i, j]

        return

    # 支点条件を整理し解析用変数にセットする
    def set_fix(self, inp: dataManager):

        mpfix = dict()
        key = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz' ]

        for fn in inp.fix_node:

            fix = np.zeros([self.nfree], dtype = np.float64)
            for i in range(self.nfree):
                if fix[i] == 1:
                    continue # すでに固定支点となっていたら何もしない
                if fn[key[i]] == 1:
                    fix[i] = 1 # 固定支点
                else:
                    fix[i] += fn[key[i]] # バネ支点
            n = fn['n']
            mpfix[n] = fix

        return mpfix

    # 全体剛性マトリックスの合計サイズを計測する
    ## この関数で Gindex が決まる
    def setMatrixIndex(self, inp: dataManager ,kmat: list, _nodeIndex: dict):

        nodeDict = dict() # 各節点の自由度
        index = dict()

        for k in kmat:
            for ID in k.nodes:
                nodeDict[ID] = k.nfree
        n = 0
        for ID in _nodeIndex:
            index[ID] = n
            n += nodeDict[ID]

        return index, n

