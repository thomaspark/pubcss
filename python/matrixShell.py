import numpy as np
from collections import OrderedDict
from matrixBar import tMatrix

from dataManager import dataManager

class shellMatrix:

    def __init__(self, inp: dataManager, ID: str):

        #解析用変数をセット
        self.label = 'Shell'
        self.nod = 4                                # 要素の節点数
        self.nfree = 6                              # 節点自由度
        self.dof = self.nod * self.nfree

        self.ID = ID

        self.shell = inp.shell[ID]

        self.nodes = {}
        for id in self.shell['nodes']:
            self.nodes[id] = inp.node[id]
        self.p = []
        for id in self.nodes:
            n = self.nodes[id]
            self.p.append([n['x'], n['y'], n['z']])

        self.eNo   = self.shell['e']

        element = inp.element[self.eNo]
        self.ee = element['E']    # elastic modulus
        self.gg = element['G']    # self.ee/2/(1+Poisson's ratio)
        self.t  = element['A']    # section area
        self.nu = (0.5 * self.ee  / self.gg) - 1 # ポアソン比

        # 四角形1次要素の積分点のξ,η座標,重み係数
        self.xg = np.array([
            [-1 / np.sqrt(3), -1 / np.sqrt(3)],
            [ 1 / np.sqrt(3), -1 / np.sqrt(3)],
            [-1 / np.sqrt(3),  1 / np.sqrt(3)],
            [ 1 / np.sqrt(3),  1 / np.sqrt(3)]])

        self.msh = self.dMatrix()

        self.ck = self.stiffnessMatrix()

        return

    # Dマトリックス
    ## ee: ヤング率 (縦弾性係数)
    ## gg: 横弾性係数
    def dMatrix(self):
        # Dマトリックス
        KS_RECT = 5 / 6
        coef = self.ee / (1 - self.nu * self.nu)
        s2 = coef * self.nu
        msh = np.array([[coef, s2, 0, 0, 0],
                        [s2, coef, 0, 0, 0],
                        [0, 0, self.gg, 0, 0],
                        [0, 0, 0, KS_RECT * self.gg, 0],
                        [0, 0, 0, 0, KS_RECT * self.gg]])

        return msh


    # 形状関数行列 [ Ni dNi/dξ dNi/dη ] を返す
    ## xsi,eta - 要素内部ξ,η座標
    def shapeFunction(self, xsi, eta):
        return np.array((
            [0.25 * (1 - xsi) * (1 - eta), -0.25 * (1 - eta), -0.25 * (1 - xsi)],
            [0.25 * (1 + xsi) * (1 - eta), 0.25 * (1 - eta), -0.25 * (1 + xsi)],
            [0.25 * (1 + xsi) * (1 + eta), 0.25 * (1 + eta), 0.25 * (1 + xsi)],
            [0.25 * (1 - xsi) * (1 + eta), -0.25 * (1 + eta), 0.25 * (1 - xsi)]
            ))


    # 法線ベクトルを返す
    # p - 頂点座標
    def normalVector(self, p):
        a = np.subtract(p[2], p[0])
        b = np.subtract(p[3], p[1])
        c = np.cross(a, b)
        d = c / np.linalg.norm(c)
        return d


    # 方向余弦マトリックスを返す
    # p - 頂点座標
    def dirVectors(self, p):
        v3 = self.normalVector(p)
        v2 = np.subtract(p[1], p[0])
        v2 = np.cross(v3, v2)
        v2 = v2 / np.linalg.norm(v2)
        v1 = np.cross(v2, v3)
        return v1, v2, v3


    # 方向余弦マトリックスを返す
    # p - 頂点座標
    # axis - 断面基準方向ベクトル
    def dirMatrix(self, p):
        v1, v2, v3 = self.dirVectors(p)
        return np.array((
            [v1[0], v2[0], v3[0]],
            [v1[1], v2[1], v3[1]],
            [v1[2], v2[2], v3[2]]
        ))


    # ヤコビ行列を返す
    ## p - 要素節点
    ## sf - 形状関数行列
    ## n - 法線ベクトル
    ## t - 要素厚さ
    def jacobianMatrix(self, p, sf, n, t):
        jac = np.zeros((3, 3))
        for i in range(self.nod):
            sfi = sf[i]
            ppi = p[i]
            pix = ppi[0]
            piy = ppi[1]
            piz = ppi[2]
            for j in range(2):
                sfij = sfi[j + 1]
                jac[0, j] += sfij * pix
                jac[1, j] += sfij * piy
                jac[2, j] += sfij * piz
        jac[0, 2] = 0.5 * t * n[0]
        jac[1, 2] = 0.5 * t * n[1]
        jac[2, 2] = 0.5 * t * n[2]

        return jac


    # 逆ヤコビ行列を返す
    ## ja - ヤコビ行列
    ## d - 方向余弦マトリックス
    def jacobInv(self, e1, d):
        jd = np.array((
            [e1[0,0] * d[0,0] + e1[1,0] * d[1,0] + e1[2,0] * d[2,0],
            e1[0,0] * d[0,1] + e1[1,0] * d[1,1] + e1[2,0] * d[2,1],
            e1[0,0] * d[0,2] + e1[1,0] * d[1,2] + e1[2,0] * d[2,2]],
            [e1[0,1] * d[0,0] + e1[1,1] * d[1,0] + e1[2,1] * d[2,0],
            e1[0,1] * d[0,1] + e1[1,1] * d[1,1] + e1[2,1] * d[2,1],
            e1[0,1] * d[0,2] + e1[1,1] * d[1,2] + e1[2,1] * d[2,2]],
            [0,
            0,
            e1[0,2] * d[0,2] + e1[1,2] * d[1,2] + e1[2,2] * d[2,2]]
            ))
        return np.linalg.inv(jd)


    # 歪 - 変位マトリックスの転置行列を返す
    ## ja - ヤコビ行列
    ## sf - 形状関数行列
    ## d - 方向余弦マトリックス
    def strainMatrix1(self, ja, sf, d):
        m = np.zeros((self.nod, 4))
        ji = self.jacobInv(ja, d)
        for i in range(self.nod):
            mi = m[i]
            sfi = sf[i]
            for j in range(3):
                mi[j] = ji[0, j] * sfi[1] + ji[1, j] * sfi[2]
            mi[3] = ji[2, 2] * sfi[0]
        return m


    # 剛性マトリックスの方向を修正する
    # d - 方向余弦マトリックス
    # k - 剛性マトリックス
    def toDir(self, d, k):
        a = np.dot(d, k)
        for i in range(len(k)):
            ki = k[i]
            ai = a[i]
            for j in range(len(ki)):
                # ki[j] = self.dotVV(ai, d[j])
                ki[j] = np.dot(ai, d[j])

    # 積分点の剛性マトリックスを返す
    # p - 要素節点
    # d1 - 応力 - 歪マトリックス
    # n - 法線ベクトル
    # xsi,eta - 要素内部ξ,η座標
    # t - 要素厚さ
    def quadstiffPart(self, p, d1, n, xsi, eta, t):
        d = self.dirMatrix(p)
        sf = self.shapeFunction(xsi, eta)
        ja = self.jacobianMatrix(p, sf, n, t)
        bc0 = self.strainMatrix1(ja, sf, d)
        sf1 = self.shapeFunction(xsi, 0)
        ja1 = self.jacobianMatrix(p, sf1, n, t)
        sf2 = self.shapeFunction(0, eta)
        ja2 = self.jacobianMatrix(p, sf2, n, t)
        bc = np.array((self.strainMatrix1(ja1, sf1, d), self.strainMatrix1(ja2, sf2, d)))
        kk = np.zeros([self.dof, self.dof])
        dete = np.linalg.det(ja)
        jacob = abs(dete)

        tt6 = t * t / 6.0
        ce1 = 1e-3 * t * t * d1[3,3]
        ce2 = -ce1 / (self.nod - 1)
        k1 = np.zeros((3,3))
        k2 = np.zeros((3,3))
        k3 = np.zeros((3,3))
        k4 = np.zeros((3,3))
        for i in range(self.nod):
            for j in range(self.nod):
                for j1 in range(3):
                    for j2 in range(3):
                        k1[j1][j2] = 0
                        k2[j1][j2] = 0
                        k3[j1][j2] = 0
                        k4[j1][j2] = 0
                for j1 in range(2):
                    for j2 in range(2):
                        k1[j1][j2] = bc0[i][j1] * d1[j1][j2] * bc0[j][j2] + bc0[i][1 - j1] * d1[2][2] * bc0[j][1 - j2]
                    dd = d1[4 - j1][4 - j1]
                    k1[j1][j1] += bc[1 - j1][i][2] * dd * bc[1 - j1][j][2]
                    k1[j1][2] = bc[1 - j1][i][2] * dd * bc[j1][j][j1]
                    k1[2][j1] = bc[j1][i][j1] * dd * bc[1 - j1][j][2]
                    k2[j1][j1] = bc[1 - j1][i][2] * dd * bc[1 - j1][j][3]
                    k2[2][j1] = bc[1 - j1][i][j1] * dd * bc[1 - j1][j][3]
                    k3[j1][j1] = bc[1 - j1][i][3] * dd * bc[1 - j1][j][2]
                    k3[j1][2] = bc[1 - j1][i][3] * dd * bc[1 - j1][j][j1]

                k1[2][2] = bc[0][i][1] * d1[3][3] * bc[0][j][1] + bc[1][i][0] * d1[4][4] * bc[1][j][0]
                k4[0][0] = k1[1][1] + 3 * bc[0][i][3] * d1[3][3] * bc[0][j][3]
                k4[0][1] = -k1[1][0]
                k4[1][0] = -k1[0][1]
                k4[1][1] = k1[0][0] + 3 * bc[1][i][3] * d1[4][4] * bc[1][j][3]
                for j1 in range(3):
                    kt = k2[j1][0]
                    k2[j1][0] = -k2[j1][1]
                    k2[j1][1] = kt
                    kt = k3[0][j1]
                    k3[0][j1] = -k3[1][j1]
                    k3[1][j1] = kt

                if (i == j):
                    k4[2][2] = ce1
                else:
                    k4[2][2] = ce2
                self.toDir(d, k1)
                self.toDir(d, k2)
                self.toDir(d, k3)
                self.toDir(d, k4)
                i0 = 6 * i
                j0 = 6 * j
                for j1 in range(3):
                    for j2 in range(3):
                        kk[i0 + j1][j0 + j2] = 2 * jacob * k1[j1][j2]
                        kk[i0 + j1][j0 + 3 + j2] = t * jacob * k2[j1][j2]
                        kk[i0 + 3 + j1][j0 + j2] = t * jacob * k3[j1][j2]
                        kk[i0 + 3 + j1][j0 + 3 + j2] = tt6 * jacob * k4[j1][j2]

        return kk


    # 行列の和を計算する
    ## a - 基準行列
    ## da - 加える行列
    def addMatrix(self, a, da):
        for i in range(len(a)):
            for j in range(len(a[i])):
                a[i, j] += da[i, j]


    # 剛性マトリックスを返す
    # p - 要素節点
    # d1 - 応力 - 歪マトリックス
    # sp - シェルパラメータ
    def stiffnessMatrix(self):

        kk = np.zeros((self.dof, self.dof))
        n = self.normalVector(self.p)

        for intP in self.xg:
            ks = self.quadstiffPart(self.p, self.msh, n, intP[0], intP[1], self.t)
            self.addMatrix(kk, ks)

        return kk

    # 変位量から断面力を計算する
    def getFsec(self, disg, Gindex):

        work = []
        for ID in self.nodes:
            npos = Gindex[ID]
            work.extend([
                disg[npos],
                disg[npos + 1],
                disg[npos + 2],
                disg[npos + 3],
                disg[npos + 4],
                disg[npos + 5]
            ])


        #全体座標系の断面力 --------------------------------------------
        fseg = np.dot(self.ck, work)

        #要素座標系の断面力 ---------------------------------------------
        list_fsec = dict()
        ## 各節点を結ぶ線を仮想Bar部材とし、断面力を求める
        list_node = self.shell['nodes']
        matrix_nodes = [fseg[0:6], fseg[6:12], fseg[12:18], fseg[18:24]]
        IDi = list_node[-1]
        matrixI = matrix_nodes[-1]
        for i in range(len(list_node)):
            matrixJ = matrix_nodes[i]
            IDj = list_node[i]
            ID = IDi + '-' + IDj
            ip = self.nodes[IDi]
            jp = self.nodes[IDj]
            # Work vector for matrix assembly
            ts = self.normalVector(self.p) 
            tm = tMatrix(ip['x'], ip['y'], ip['z'], jp['x'], jp['y'], jp['z'], 0)
            # ２つのベクトルの成す角を求める。
            theta = self.tangent_angle(ts, tm.t3[2])
            if theta != 0:
                tm = tMatrix(ip['x'], ip['y'], ip['z'], jp['x'], jp['y'], jp['z'], theta)

            fsegl = np.concatenate([matrixI, matrixJ])
            fsec = np.dot(tm.tt, fsegl)

            # 構造力学の定義
            for i in [0, 3, 5, 7, 8, 10]:
                fsec[i] *= -1

            IDi = IDj
            matrixI = matrixJ
            list_fsec[ID] = fsec

        return list_fsec, fseg

    # ２つのベクトルの成す角を求める。
    def tangent_angle(self, ts, tz):
        x = np.inner(ts, tz) # 内積を計算。
        theta = np.arccos(x/(ts*tz))
        return theta[2]