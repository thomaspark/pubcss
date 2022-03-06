# ======================
# 3D Frame Analysis
# ======================
import numpy as np
from collections import OrderedDict

from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix


from matrixBar import barMatrix
from matrixShell import shellMatrix
from matrixG import gMatrix

from dataInput import dataInput
from dataManager import dataManager


class FrameCalc:

    def __init__(self, _inp):
        self.inp = dataManager(_inp)

        # 全ての荷重ケースの解析を開始する
    def calcrate(self):

        result = dict() # 結果を格納するdictionary

        # calcration
        for id in self.inp.loadList:
            error, res = self.calcLinear3D(id)
            if error != None:
                return error, id
            result[id] = res
        return None, result


    # ターゲットケースの情報をセットする
    def calcLinear3D(self, caseID: str, debug=False):

        # ターゲットケースの基本情報をセットする
        error = self.inp.setTargetData(caseID)
        if error != None:
            return error, None

        # 節点に剛性マトリックスのインデックスを割り当てる
        self.nodeIndex = self.reNumbering()

        # 全体剛性マトリックス
        gmat = gMatrix(self.inp, self.nodeIndex)

        # X, LU = torch.solve(fp, gk)
        d = spsolve(gmat.gk, gmat.fp, use_umfpack=True)
        # result = np.linalg.solve(gmat.gk, fmat.fp)

        # 解析結果を集計する
        disg, reac, fsec, shell_fsec = self.getResultData(gmat, d)

        return None, self.inp.getResultData(gmat, disg, reac, fsec, shell_fsec)


    # 節点に剛性マトリックスのインデックスを割り当てる
    def reNumbering(self):

        # 節点に剛性マトリックスのインデックスを割り当てる
        nodeIndex = list(self.inp.node.keys())
        nd = dict(zip(nodeIndex, list(range(len(nodeIndex)))))

        members = self.inp.member.values()

        # i端節点No ＞ ｊ端節点No の場合は i端 と j端をひっくり返す・・・この処理重いいらねんじゃね？
        flg = True
        while flg:
            flg = False
            for m in members:
                i = nd[m['ni']]
                j = nd[m['nj']]
                if i > j:
                    nodeIndex[j] = m['ni']
                    nodeIndex[i] = m['nj']
                    nd[m['ni']] = j
                    nd[m['nj']] = i
                    flg = True
                    break
        return nodeIndex


    def getResultData(self, gmat: gMatrix, disg):

        reac = np.zeros_like(disg)  # Section force vector
        fsec = dict()  # Section force vector
        shell_fsec = dict()  # Section force vector

        # calculation of section force
        for k in gmat.kmat:
            if k.label == 'Shell':
                list_fsec, fseg = k.getFsec(disg, gmat.Gindex) # Shell部材の剛性マトリックス
                for ID in list_fsec:
                    shell_fsec[ID] = list_fsec[ID]
                nodes = k.nodes

            elif k.label == 'Bar':
                fsec[k.ID], fseg = k.getFsec(disg, gmat.Gindex) # Bar部材の剛性マトリックス
                nodes = [k.IDi, k.IDj]

            # 支点の反力
            i = 0
            for nd in nodes:
                if nd in gmat.mpfix:
                    fix = gmat.mpfix[nd]
                    for j in range(k.nfree):
                        if fix[j] != 0:
                            # 固定支点がある場合
                            iz = gmat.Gindex[nd] + j
                            reac[iz] -= fseg[j + i]
                i += k.nfree



        return disg, reac, fsec, shell_fsec
