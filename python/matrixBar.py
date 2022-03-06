import numpy as np
from collections import OrderedDict
import copy

class barMatrix:

    def __init__(self, inp, ID: str):

        #解析用変数をセット
        # info
        self.label = 'Bar'
        self.nod  = 2   # Number of nodes per member
        self.nfree= 6   # Degree of freedom per node
        self.ID = ID

        member = inp.member[ID]

        # node
        self.IDi   = member['ni']
        self.IDj   = member['nj']
        ip = inp.node[self.IDi]
        jp = inp.node[self.IDj]

        self.nodes = {
            self.IDi: ip,
            self.IDj: jp
        }


        # 部材長さ
        self.el = inp.GetDistance(ip, jp) #入力データに格納されている部材長さ'L'は 1000倍しているのであらためて計算

        # element
        self.eNo   = member['e']
        element    = inp.element[self.eNo]
        self.ee    = element['E']    # elastic modulus
        self.jj    = element['G']    # self.ee/2/(1+Poisson's ratio)
        self.aa    = element['A']    # section area
        self.aix   = element['J']    # tortional constant
        self.aiy   = element['Iy']   # moment of inertia around y-axis
        self.aiz   = element['Iz']   # moment of inertia around z-axis
        self.alpha = element['Xp']   # alpha : thermal expansion coefficient

        self.GJ  = self.jj * self.aix
        self.EA  = self.ee * self.aa
        self.EIy = self.ee * self.aiy
        self.EIz = self.ee * self.aiz

        # joint
        self.k1x   = member['xi']
        self.k1y   = member['yi']
        self.k1z   = member['zi']
        self.k2x   = member['xj']
        self.k2y   = member['yj']
        self.k2z   = member['zj']
        self.theta = member['cg'] # theta : chord angle

        # MemberSupport
        self.tx  = 0
        self.ty  = 0
        self.tz  = 0
        self.tr  = 0
        for fm in inp.fix_member:
            if ID  ==  fm['m']:
                self.tx = fm['tx']  # fixed in x-direction
                self.ty = fm['ty']  # fixed in x-direction
                self.tz = fm['tz']  # fixed in x-direction
                self.tr = fm['tr']  # fixed in x-direction
                break

        # Work vector for matrix assembly
        self.tm = tMatrix(ip['x'], ip['y'], ip['z'], jp['x'], jp['y'], jp['z'], self.theta)
        self.ek = self.SM_3DFRM()                                   # 要素座標系
        self.t3 = self.tm.t3                                        # 座標変換行列
        self.tt = self.tm.tt                                        # 座標変換行列
        self.gk = np.dot( np.dot( self.tt.T, self.ek ), self.tt )   # 全体座標系
        self.ck = self.S2_3DFRM( copy.copy( self.gk ) )             # ピン結合を考慮した 全体座標系
        self.ek = np.dot( np.dot( self.tt, self.ck ), self.tt.T )   # ピン結合を考慮した 要素座標系

        # loaded member
        fe = self.set_load(inp)
        self.wxi = fe[0]
        self.wxj = fe[1]
        self.wyi = fe[2]
        self.wyj = fe[3]
        self.wzi = fe[4]
        self.wzj = fe[5]
        self.wti = fe[6]
        self.wtj = fe[7]
        self.wte = fe[8] # Temperature change of element


        # 分布荷重を節点荷重に変換する
        self.tfe, self.wfe = self.Convert_memberload_to_nodeload()

        return


    def SM_3DFRM(self):
        ek = np.zeros([12,12],dtype=np.float64) # local stiffness matrix

        # 節点に作用する軸方向力
        #fxi =
        ek[0, 0] = self.ka(self.tx)
        ek[6, 0] = -self.ka(-self.tx)

        #fxj =
        ek[0, 6] = -self.ka(-self.tx)
        ek[6, 6] = self.ka(self.tx)


        #節点に作用するy方向のせん断力
        #fyi = 12EIz/l^3 × uyi + 6EIz/l^2 × θzi - 12EIz/l^3 × uyj + 6EIz/l^2 × θzj
        ek[1, 1] = self.kb(self.ty)
        ek[5, 1] = self.kc(self.ty)
        ek[7, 1] = -self.kb(-self.ty)
        ek[11, 1] = self.kc(-self.ty)

        #fyj = -12EIz/l^3 × uyi - 6EIz/l^2 × θzi + 12EIz/l^3 × uyj - 6EIz/l^2 × θzj
        ek[1, 7] = -self.kb(-self.ty)
        ek[5, 7] = -self.kc(-self.ty)
        ek[7, 7] = self.kb(self.ty)
        ek[11, 7] = -self.kc(self.ty)


        #節点に作用するz方向のせん断力
        #fzi = 12EIy/l^3 × uzi + 6EIy/l^2 × θyi - 12EIy/l^3 × uzj + 6EIy/l^2 × θyj
        ek[2, 2] = self.kd(self.tz)
        ek[4, 2] = -self.ke(self.tz)
        ek[8, 2] = -self.kd(-self.tz)
        ek[10, 2] = -self.ke(-self.tz)

        #fzj = -12EIy/l^3 × uzi - 6EIy/l^2 × θyi + 12EIy/l^3 × uzj - 6EIy/l^2 × θyj
        ek[2, 8] = -self.kd(-self.tz)
        ek[4, 8] = self.ke(-self.tz)
        ek[8, 8] = self.kd(self.tz)
        ek[10, 8] = self.ke(self.tz)


        #節点に作用するx軸回りのねじり
        #Mxi = GJ/l × θxi - GJ/l × θxj
        ek[3, 3] = self.kf(self.tr)
        ek[9, 3] = -self.kf(-self.tr)

        #Mxj = -GJ/l × θxi + GJ/l × θxj
        ek[3, 9] = -self.kf(-self.tr)
        ek[9, 9] = self.kf(self.tr)


        #節点に作用するy軸回りの曲げモーメント
        #Myi = -6EIy/l^2 × uzi + 4EIy/l × θyi + 6EIy/l^2 × uzj + 2EIy/l × θyj
        ek[2, 4] = -self.ke(self.tz)
        ek[4, 4] = self.kg(self.tz)
        ek[8, 4] = self.ke(-self.tz)
        ek[10, 4] = self.kh(self.tz)

        #Myj = -6EIy/l^2 × uzi + 2EIy/l × θyi + 6EIy/l^2 × uzj + 4EIy/l × θyj
        ek[2, 10] = -self.ke(-self.tz)
        ek[4, 10] = self.kh(self.tz)
        ek[8, 10] = self.ke(self.tz)
        ek[10, 10] = self.kg(self.tz)


        #節点に作用するz軸回りの曲げモーメント
        #Mzi = 6EIz/l^2 × uyi + 4EIz/l × θzi - 6EIz/l^2 × uyj + 2EIz/l × θzj
        ek[1, 5] = self.kc(self.ty)
        ek[5, 5] = self.ki(self.ty)
        ek[7, 5] = -self.kc(-self.ty)
        ek[11, 5] = self.kj(self.ty)

        #Mzj = 6EIz/l^2 × uyi + 2EIz/l × θzi - 6EIz/l^2 × uyj + 4EIz/l × θzj
        ek[1, 11] = self.kc(-self.ty)
        ek[5, 11] = self.kj(self.ty)
        ek[7, 11] = -self.kc(self.ty)
        ek[11, 11] = self.ki(self.ty)

        return ek

    def ka(self, Ku = 0.0):
        if Ku == 0:
            result = self.EA / self.el
        else:
            L0 = np.sqrt(abs(Ku) / (self.EA))
            E0 = np.exp(L0 * self.el)
            C8 = (E0 + 1 / E0) / 2
            S8 = (E0 - 1 / E0) / 2
            if Ku > 0:
                result = self.EA * L0 * C8 / S8
            if Ku < 0:
                result = self.EA * L0 / S8

        return result

    def kb(self, Kv = 0.0):
        if Kv == 0:
            result = 12 * self.EIz / (self.el**3)
        else:
            L0 = np.sqrt(np.sqrt(4 * self.EIz / abs(Kv)))
            L9 = self.el / L0
            E0 = np.exp(L9)
            C8 = (E0 + 1 / E0) / 2
            S8 = (E0 - 1 / E0) / 2
            C9 = np.cos(L9)
            S9 = np.sin(L9)
            if Kv > 0:
                result = (4 * self.EIz / L0**3) * ((S8 * C8 + S9 * C9) / (S8**2 - S9**2))
            if Kv < 0:
                result = (4 * self.EIz / L0**3) * ((C8 * S9 + S8 * C9) / (S8**2 - S9**2))

        return result

    def kc(self, Kv = 0.0):

        if Kv == 0:
            result = 6 * self.EIz / (self.el**2)
        else:
            L0 = np.sqrt(np.sqrt(4 * self.EIz / abs(Kv)))
            L9 = self.el / L0
            E0 = np.exp(L9)
            C8 = (E0 + 1 / E0) / 2
            S8 = (E0 - 1 / E0) / 2
            C9 = np.cos(L9)
            S9 = np.sin(L9)
            if Kv > 0:
                result = (2 * self.EIz / L0**2) * ((S8**2 + S9**2) / (S8**2 - S9**2))
            if Kv < 0:
                result = (4 * self.EIz / L0**2) * ((S8 * S9) / (S8**2 - S9**2))

        return result

    def kd(self, Kv = 0.0):

        if Kv == 0:
            result = 12 * self.EIy / (self.el**3)
        else:
            L0 = np.sqrt(np.sqrt(4 * self.EIy / abs(Kv)))
            L9 = self.el / L0
            E0 = np.exp(L9)
            C8 = (E0 + 1 / E0) / 2
            S8 = (E0 - 1 / E0) / 2
            C9 = np.cos(L9)
            S9 = np.sin(L9)
            if Kv > 0:
                result = (4 * self.EIy / L0**3) * ((S8 * C8 + S9 * C9) / (S8**2 - S9**2))
            if Kv < 0:
                result = (4 * self.EIy / L0**3) * ((C8 * S9 + S8 * C9) / (S8**2 - S9**2))

        return result

    def ke(self, Kv = 0.0):

        if Kv == 0:
            result = 6 * self.EIy / (self.el**2)
        else:
            L0 = np.sqrt(np.sqrt(4 * self.EIy / abs(Kv)))
            L9 = self.el / L0
            E0 = np.exp(L9)
            C8 = (E0 + 1 / E0) / 2
            S8 = (E0 - 1 / E0) / 2
            C9 = np.cos(L9)
            S9 = np.sin(L9)
            if Kv > 0:
                result = (2 * self.EIy / L0**2) * ((S8**2 + S9**2) / (S8**2 - S9**2))
            if Kv < 0:
                result = (4 * self.EIy / L0**2) * ((S8 * S9) / (S8**2 - S9**2))

        return result

    def kf(self, Kt = 0.0):

        if Kt == 0:
            result = self.GJ / self.el
        else:
            w = np.sqrt(abs(Kt) / self.GJ)
            if Kt > 0:
                result = self.GJ * w * np.cosh(w * self.el) / np.sinh(w * self.el)
            if Kt < 0:
                result = self.GJ * w / np.sinh(w * self.el)

        return result

    def kg(self, Kv = 0.0):

        if Kv == 0:
            result = 4 * self.EIy / self.el
        else:
            L0 = np.sqrt(np.sqrt(4 * self.EIy / abs(Kv)))
            L9 = self.el / L0
            E0 = np.exp(L9)
            C8 = (E0 + 1 / E0) / 2
            S8 = (E0 - 1 / E0) / 2
            C9 = np.cos(L9)
            S9 = np.sin(L9)
            result = (2 * self.EIy / L0) * ((S8 * C8 - S9 * C9) / (S8**2 - S9**2))

        return result

    def kh(self, Kv = 0.0):

        if Kv == 0:
            result = 2 * self.EIy / self.el
        else:
            L0 = np.sqrt(np.sqrt(4 * self.EIy / abs(Kv)))
            L9 = self.el / L0
            E0 = np.exp(L9)
            C8 = (E0 + 1 / E0) / 2
            S8 = (E0 - 1 / E0) / 2
            C9 = np.cos(L9)
            S9 = np.sin(L9)
            result = (2 * self.EIy / L0) * ((C8 * S9 - S8 * C9) / (S8**2 - S9**2))

        return result

    def ki(self, Kv = 0.0):

        if Kv == 0:
            result = 4 * self.EIz / self.el
        else:
            L0 = np.sqrt(np.sqrt(4 * self.EIz / abs(Kv)))
            L9 = self.el / L0
            E0 = np.exp(L9)
            C8 = (E0 + 1 / E0) / 2
            S8 = (E0 - 1 / E0) / 2
            C9 = np.cos(L9)
            S9 = np.sin(L9)
            result = (2 * self.EIz / L0) * ((S8 * C8 - S9 * C9) / (S8**2 - S9**2))

        return result

    def kj(self, Kv = 0.0):

        if Kv == 0:
            result = 2 * self.EIz / self.el
        else:
            L0 = np.sqrt(np.sqrt(4 * self.EIz / abs(Kv)))
            L9 = self.el / L0
            E0 = np.exp(L9)
            C8 = (E0 + 1 / E0) / 2
            S8 = (E0 - 1 / E0) / 2
            C9 = np.cos(L9)
            S9 = np.sin(L9)
            result = (2 * self.EIz / L0) * ((C8 * S9 - S8 * C9) / (S8**2 - S9**2))

        return result

    def S2_3DFRM(self, ek):

        #1端1軸まわりのみピンをもつ要素の場合
        if (self.k1x == 0 and self.k1y == 1 and self.k1z == 1) and (self.k2x == 1 and self.k2y == 1 and self.k2z == 1):
            ek = self.elka1(ek,3)
        elif (self.k1x == 1 and self.k1y == 0 and self.k1z == 1) and (self.k2x == 1 and self.k2y == 1 and self.k2z == 1):
            ek =self.elka1(ek,4)
        elif (self.k1x == 1 and self.k1y == 1 and self.k1z == 0) and (self.k2x == 1 and self.k2y == 1 and self.k2z == 1):
            ek =self.elka1(ek,5)
        elif (self.k1x == 1 and self.k1y == 1 and self.k1z == 1) and (self.k2x == 0 and self.k2y == 1 and self.k2z == 1):
            ek =self.elka1(ek,9)
        elif (self.k1x == 1 and self.k1y == 1 and self.k1z == 1) and (self.k2x == 1 and self.k2y == 0 and self.k2z == 1):
            ek =self.elka1(ek,10)
        elif (self.k1x == 1 and self.k1y == 1 and self.k1z == 1) and (self.k2x == 1 and self.k2y == 1 and self.k2z == 0):
            ek =self.elka1(ek,11)

        #両端1軸まわりのみピンをもつ要素の場合
        elif (self.k1x == 0 and self.k1y == 1 and self.k1z == 1) and (self.k2x == 0 and self.k2y == 1 and self.k2z == 1):
            ek =self.elka2(ek,3)
        elif (self.k1x == 1 and self.k1y == 0 and self.k1z == 1) and (self.k2x == 1 and self.k2y == 0 and self.k2z == 1):
            ek =self.elka2(ek,4)
        elif (self.k1x == 1 and self.k1y == 1 and self.k1z == 0) and (self.k2x == 1 and self.k2y == 1 and self.k2z == 0):
            ek =self.elka2(ek,5)

        #1端2軸まわりのみピンをもつ要素の場合
        elif (self.k1x == 0 and self.k1y == 0 and self.k1z == 1) and (self.k2x == 1 and self.k2y == 1 and self.k2z == 1):
            ek =self.elka3(ek,5)
        elif (self.k1x == 0 and self.k1y == 1 and self.k1z == 0) and (self.k2x == 1 and self.k2y == 1 and self.k2z == 1):
            ek =self.elka3(ek,4)
        elif (self.k1x == 1 and self.k1y == 0 and self.k1z == 0) and (self.k2x == 1 and self.k2y == 1 and self.k2z == 1):
            ek =self.elka3(ek,3)
        elif (self.k1x == 1 and self.k1y == 1 and self.k1z == 1) and (self.k2x == 0 and self.k2y == 0 and self.k2z == 1):
            ek =self.elka3(ek,11)
        elif (self.k1x == 1 and self.k1y == 1 and self.k1z == 1) and (self.k2x == 0 and self.k2y == 1 and self.k2z == 0):
            ek =self.elka3(ek,10)
        elif (self.k1x == 1 and self.k1y == 1 and self.k1z == 1) and (self.k2x == 1 and self.k2y == 0 and self.k2z == 0):
            ek =self.elka3(ek,9)

        #1端3軸まわりのみピンをもつ要素の場合
        elif (self.k1x == 0 and self.k1y == 0 and self.k1z == 0) and (self.k2x == 1 and self.k2y == 1 and self.k2z == 1):
            ek =self.elka4(ek,3)
        elif (self.k1x == 1 and self.k1y == 1 and self.k1z == 1) and (self.k2x == 0 and self.k2y == 0 and self.k2z == 0):
            ek =self.elka4(ek,9)

        #両端2軸まわりのみピンをもつ要素の場合
        elif (self.k1x == 0 and self.k1y == 0 and self.k1z == 1) and (self.k2x == 0 and self.k2y == 0 and self.k2z == 1):
            ek =self.elka5(ek,0)
        elif (self.k1x == 0 and self.k1y == 1 and self.k1z == 0) and (self.k2x == 0 and self.k2y == 1 and self.k2z == 0):
            ek =self.elka5(ek,1)
        elif (self.k1x == 1 and self.k1y == 0 and self.k1z == 0) and (self.k2x == 1 and self.k2y == 0 and self.k2z == 0):
            ek =self.elka5(ek,2)

        return ek

    def elka1(self, ek, jj):
        ####################################################
        ##     トラス（０）　ラーメン（１）の結合状態     ##
        ##        self.kTR1[i] - self.kTR2[j)             ##
        ##                    x,y,z - x,y,z               ##
        ##             jj=4  (0,1,1 - 1,1,1)              ##
        ##             jj=5  (1,0,1 - 1,1,1)              ##
        ##     ３次元  jj=6  (1,1,0 - 1,1,1)              ##
        ##             jj=10 (1,1,1 - 0,1,1)              ##
        ##             jj=11 (1,1,1 - 1,0,1)              ##
        ##             jj=12 (1,1,1 - 1,1,0)              ##
        ##                ** ELkA1.for **                 ##
        ####################################################
        b1  = ek[jj,0]  / ek[jj,jj]
        b2  = ek[jj,1]  / ek[jj,jj]
        b3  = ek[jj,2]  / ek[jj,jj]
        b4  = ek[jj,3]  / ek[jj,jj]
        b5  = ek[jj,4]  / ek[jj,jj]
        b6  = ek[jj,5]  / ek[jj,jj]
        b7  = ek[jj,6]  / ek[jj,jj]
        b8  = ek[jj,7]  / ek[jj,jj]
        b9  = ek[jj,8]  / ek[jj,jj]
        b10 = ek[jj,9]  / ek[jj,jj]
        b11 = ek[jj,10] / ek[jj,jj]
        b12 = ek[jj,11] / ek[jj,jj]

        for i in range(12):
            ekk = ek[i,jj]
            ek[i,0]  = ek[i,0]  - ekk * b1
            ek[i,1]  = ek[i,1]  - ekk * b2
            ek[i,2]  = ek[i,2]  - ekk * b3
            ek[i,3]  = ek[i,3]  - ekk * b4
            ek[i,4]  = ek[i,4]  - ekk * b5
            ek[i,5]  = ek[i,5]  - ekk * b6
            ek[i,6]  = ek[i,6]  - ekk * b7
            ek[i,7]  = ek[i,7]  - ekk * b8
            ek[i,8]  = ek[i,8]  - ekk * b9
            ek[i,9]  = ek[i,9]  - ekk * b10
            ek[i,10] = ek[i,10] - ekk * b11
            ek[i,11] = ek[i,11] - ekk * b12

        return ek

    def elka2(self, ek, j1):
        ####################################################
        ##     トラス（０）　ラーメン（１）の結合状態     ##
        ##        self.kTR1[i] - self.kTR2[j)             ##
        ##                    x,y,z - x,y,z               ##
        ##             j1=4  (0,1,1 - 0,1,1)              ##
        ##     ３次元  j1=5  (1,0,1 - 1,0,1)              ##
        ##             j1=6  (1,1,0 - 1,1,0)              ##
        ##                ** ELkA2.for **                 ##
        ####################################################
        j2 = j1 + 6
        bunbo1 = ek[j1,j1] * ek[j2,j2] - ek[j2,j1] * ek[j1,j2]
        bunbo2 = -bunbo1

        b11  = (ek[j1,j2] * ek[j2,0]  - ek[j2,j2] * ek[j1,0])  / bunbo1
        b12  = (ek[j1,j2] * ek[j2,1]  - ek[j2,j2] * ek[j1,1])  / bunbo1
        b13  = (ek[j1,j2] * ek[j2,2]  - ek[j2,j2] * ek[j1,2])  / bunbo1
        b14  = (ek[j1,j2] * ek[j2,3]  - ek[j2,j2] * ek[j1,3])  / bunbo1
        b15  = (ek[j1,j2] * ek[j2,4]  - ek[j2,j2] * ek[j1,4])  / bunbo1
        b16  = (ek[j1,j2] * ek[j2,5]  - ek[j2,j2] * ek[j1,5])  / bunbo1
        b17  = (ek[j1,j2] * ek[j2,6]  - ek[j2,j2] * ek[j1,6])  / bunbo1
        b18  = (ek[j1,j2] * ek[j2,7]  - ek[j2,j2] * ek[j1,7])  / bunbo1
        b19  = (ek[j1,j2] * ek[j2,8]  - ek[j2,j2] * ek[j1,8])  / bunbo1
        b110 = (ek[j1,j2] * ek[j2,9]  - ek[j2,j2] * ek[j1,9])  / bunbo1
        b111 = (ek[j1,j2] * ek[j2,10] - ek[j2,j2] * ek[j1,10]) / bunbo1
        b112 = (ek[j1,j2] * ek[j2,11] - ek[j2,j2] * ek[j1,11]) / bunbo1

        b21  = (ek[j1,j1] * ek[j2,0]  - ek[j2,j1] * ek[j1,0])  / bunbo2
        b22  = (ek[j1,j1] * ek[j2,1]  - ek[j2,j1] * ek[j1,1])  / bunbo2
        b23  = (ek[j1,j1] * ek[j2,2]  - ek[j2,j1] * ek[j1,2])  / bunbo2
        b24  = (ek[j1,j1] * ek[j2,3]  - ek[j2,j1] * ek[j1,3])  / bunbo2
        b25  = (ek[j1,j1] * ek[j2,4]  - ek[j2,j1] * ek[j1,4])  / bunbo2
        b26  = (ek[j1,j1] * ek[j2,5]  - ek[j2,j1] * ek[j1,5])  / bunbo2
        b27  = (ek[j1,j1] * ek[j2,6]  - ek[j2,j1] * ek[j1,6])  / bunbo2
        b28  = (ek[j1,j1] * ek[j2,7]  - ek[j2,j1] * ek[j1,7])  / bunbo2
        b29  = (ek[j1,j1] * ek[j2,8]  - ek[j2,j1] * ek[j1,8])  / bunbo2
        b210 = (ek[j1,j1] * ek[j2,9]  - ek[j2,j1] * ek[j1,9])  / bunbo2
        b211 = (ek[j1,j1] * ek[j2,10] - ek[j2,j1] * ek[j1,10]) / bunbo2
        b212 = (ek[j1,j1] * ek[j2,11] - ek[j2,j1] * ek[j1,11]) / bunbo2

        for i in range(12):
            ek1 = ek[i,j1]
            ek2 = ek[i,j2]
            ek[i,0]  = ek[i,0]  + ek1 * b11  + ek2 * b21
            ek[i,1]  = ek[i,1]  + ek1 * b12  + ek2 * b22
            ek[i,2]  = ek[i,2]  + ek1 * b13  + ek2 * b23
            ek[i,3]  = ek[i,3]  + ek1 * b14  + ek2 * b24
            ek[i,4]  = ek[i,4]  + ek1 * b15  + ek2 * b25
            ek[i,5]  = ek[i,5]  + ek1 * b16  + ek2 * b26
            ek[i,6]  = ek[i,6]  + ek1 * b17  + ek2 * b27
            ek[i,7]  = ek[i,7]  + ek1 * b18  + ek2 * b28
            ek[i,8]  = ek[i,8]  + ek1 * b19  + ek2 * b29
            ek[i,9]  = ek[i,9]  + ek1 * b110 + ek2 * b210
            ek[i,10] = ek[i,10] + ek1 * b111 + ek2 * b211
            ek[i,11] = ek[i,11] + ek1 * b112 + ek2 * b212

        return ek

    def elka3(self, ek, jj):
        ####################################################
        ##     トラス（０）　ラーメン（１）の結合状態     ##
        ##        self.kTR1[i] - self.kTR2[j)             ##
        ##                    x,y,z - x,y,z               ##
        ##             jj=6  (0,0,1 - 1,1,1)              ##
        ##             jj=5  (0,1,0 - 1,1,1)              ##
        ##             jj=4  (1,0,0 - 1,1,1)              ##
        ##     ３次元  jj=12 (1,1,1 - 0,0,1)              ##
        ##             jj=11 (1,1,1 - 0,1,0)              ##
        ##             jj=10 (1,1,1 - 1,0,0)              ##
        ##                ** ELkA3.for **                 ##
        ####################################################

        if jj == 5:
            j1 = 3
            j2 = 4
        elif jj == 4:
            j1 = 3
            j2 = 5
        elif jj == 3:
            j1 = 4
            j2 = 5
        elif jj == 11:
            j1 = 9
            j2 = 10
        elif jj == 10:
            j1 = 9
            j2 = 11
        elif jj == 9:
            j1 = 10
            j2 = 11
        else:
            return

        bunbo1 = ek[j1,j1] * ek[j2,j2] - ek[j2,j1] * ek[j1,j2]
        bunbo2 = -bunbo1
        b11  = (ek[j1,j2] * ek[j2,0]  - ek[j2,j2] * ek[j1,0])  / bunbo1
        b12  = (ek[j1,j2] * ek[j2,1]  - ek[j2,j2] * ek[j1,1])  / bunbo1
        b13  = (ek[j1,j2] * ek[j2,2]  - ek[j2,j2] * ek[j1,2])  / bunbo1
        b14  = (ek[j1,j2] * ek[j2,3]  - ek[j2,j2] * ek[j1,3])  / bunbo1
        b15  = (ek[j1,j2] * ek[j2,4]  - ek[j2,j2] * ek[j1,4])  / bunbo1
        b16  = (ek[j1,j2] * ek[j2,5]  - ek[j2,j2] * ek[j1,5])  / bunbo1
        b17  = (ek[j1,j2] * ek[j2,6]  - ek[j2,j2] * ek[j1,6])  / bunbo1
        b18  = (ek[j1,j2] * ek[j2,7]  - ek[j2,j2] * ek[j1,7])  / bunbo1
        b19  = (ek[j1,j2] * ek[j2,8]  - ek[j2,j2] * ek[j1,8])  / bunbo1
        b110 = (ek[j1,j2] * ek[j2,9]  - ek[j2,j2] * ek[j1,9])  / bunbo1
        b111 = (ek[j1,j2] * ek[j2,10] - ek[j2,j2] * ek[j1,10]) / bunbo1
        b112 = (ek[j1,j2] * ek[j2,11] - ek[j2,j2] * ek[j1,11]) / bunbo1

        b21  = (ek[j1,j1] * ek[j2,0]  - ek[j2,j1] * ek[j1,0])  / bunbo2
        b22  = (ek[j1,j1] * ek[j2,1]  - ek[j2,j1] * ek[j1,1])  / bunbo2
        b23  = (ek[j1,j1] * ek[j2,2]  - ek[j2,j1] * ek[j1,2])  / bunbo2
        b24  = (ek[j1,j1] * ek[j2,3]  - ek[j2,j1] * ek[j1,3])  / bunbo2
        b25  = (ek[j1,j1] * ek[j2,4]  - ek[j2,j1] * ek[j1,4])  / bunbo2
        b26  = (ek[j1,j1] * ek[j2,5]  - ek[j2,j1] * ek[j1,5])  / bunbo2
        b27  = (ek[j1,j1] * ek[j2,6]  - ek[j2,j1] * ek[j1,6])  / bunbo2
        b28  = (ek[j1,j1] * ek[j2,7]  - ek[j2,j1] * ek[j1,7])  / bunbo2
        b29  = (ek[j1,j1] * ek[j2,8]  - ek[j2,j1] * ek[j1,8])  / bunbo2
        b210 = (ek[j1,j1] * ek[j2,9]  - ek[j2,j1] * ek[j1,9])  / bunbo2
        b211 = (ek[j1,j1] * ek[j2,10] - ek[j2,j1] * ek[j1,10]) / bunbo2
        b212 = (ek[j1,j1] * ek[j2,11] - ek[j2,j1] * ek[j1,11]) / bunbo2

        for i in range(12):
            ek1 = ek[i,j1]
            ek2 = ek[i,j2]
            ek[i,0]  = ek[i,0]  + ek1 * b11  + ek2 * b21
            ek[i,1]  = ek[i,1]  + ek1 * b12  + ek2 * b22
            ek[i,2]  = ek[i,2]  + ek1 * b13  + ek2 * b23
            ek[i,3]  = ek[i,3]  + ek1 * b14  + ek2 * b24
            ek[i,4]  = ek[i,4]  + ek1 * b15  + ek2 * b25
            ek[i,5]  = ek[i,5]  + ek1 * b16  + ek2 * b26
            ek[i,6]  = ek[i,6]  + ek1 * b17  + ek2 * b27
            ek[i,7]  = ek[i,7]  + ek1 * b18  + ek2 * b28
            ek[i,8]  = ek[i,8]  + ek1 * b19  + ek2 * b29
            ek[i,9]  = ek[i,9]  + ek1 * b110 + ek2 * b210
            ek[i,10] = ek[i,10] + ek1 * b111 + ek2 * b211
            ek[i,11] = ek[i,11] + ek1 * b112 + ek2 * b212

        return ek

    def elka4(self, ek, j):
        ####################################################
        ##     トラス（０）　ラーメン（１）の結合状態     ##
        ##        self.kTR1[i] - self.kTR2[j)             ##
        ##                    x,y,z - x,y,z               ##
        ##     ３次元  j=4   (0,0,0 - 1,1,1)              ##
        ##             j=10  (1,1,1 - 0,0,0)              ##
        ##                ** ELkA4.for **                 ##
        ####################################################
        j1  = j
        j2  = j + 1
        j3  = j + 2
        b1  = ek[j1,j1]
        b2  = ek[j1,j2]
        b3  = ek[j1,j3]
        b4  = ek[j2,j1]
        b5  = ek[j2,j2]
        b6  = ek[j2,j3]
        b7  = ek[j3,j1]
        b8  = ek[j3,j2]
        b9  = ek[j3,j3]

        bb1 = 1 / b9
        am1 = b1 - b3 * b7 * bb1
        am2 = b2 - b3 * b8 * bb1
        am3 = b4 - b6 * b7 * bb1
        am4 = b5 - b6 * b8 * bb1
        am  = 1 / (am1 * am4 - am3 * am2)
        s1  = -am4 * am
        s2  = am2 * am
        s3  = b3 * am4 * bb1 * am - b6 * am2 * bb1 * am

        am  = 1 / (am2 * am3 - am4 * am1)
        t1  = -am3 * am
        T2  = am1 * am
        t3  = b3 * am3 * bb1 * am - b6 * am1 * bb1 * am

        am  = 1 / (am2 * am3 - am4 * am1)
        bm  = -am
        r1  = b8 * am3 * bb1 * am + b7 * am4 * bb1 * bm
        r2  = -b7 * am2 * bb1 * bm - b8 * am1 * bb1 * am
        bb2 = bb1 / b9
        rr1 = bb1 + b3 * b7 * bb2 * am4 * bm
        rr2 = b6 * b7 * bb2 * am2 * bm
        rr3 = b3 * b8 * bb2 * am3 * am
        rr4 = b6 * b8 * bb2 * am1 * am
        r3  = -rr1 + rr2 - rr3 + rr4

        bx1  = s1 * ek[j1,0]  + s2 * ek[j2,0]  + s3 * ek[j3,0]
        bx2  = s1 * ek[j1,1]  + s2 * ek[j2,1]  + s3 * ek[j3,1]
        bx3  = s1 * ek[j1,2]  + s2 * ek[j2,2]  + s3 * ek[j3,2]
        bx4  = s1 * ek[j1,3]  + s2 * ek[j2,3]  + s3 * ek[j3,3]
        bx5  = s1 * ek[j1,4]  + s2 * ek[j2,4]  + s3 * ek[j3,4]
        bx6  = s1 * ek[j1,5]  + s2 * ek[j2,5]  + s3 * ek[j3,5]
        bx7  = s1 * ek[j1,6]  + s2 * ek[j2,6]  + s3 * ek[j3,6]
        bx8  = s1 * ek[j1,7]  + s2 * ek[j2,7]  + s3 * ek[j3,7]
        bx9  = s1 * ek[j1,8]  + s2 * ek[j2,8]  + s3 * ek[j3,8]
        bx10 = s1 * ek[j1,9]  + s2 * ek[j2,9]  + s3 * ek[j3,9]
        bx11 = s1 * ek[j1,10] + s2 * ek[j2,10] + s3 * ek[j3,10]
        bx12 = s1 * ek[j1,11] + s2 * ek[j2,11] + s3 * ek[j3,11]

        by1  = t1 * ek[j1,0]  + T2 * ek[j2,0]  + t3 * ek[j3,0]
        by2  = t1 * ek[j1,1]  + T2 * ek[j2,1]  + t3 * ek[j3,1]
        by3  = t1 * ek[j1,2]  + T2 * ek[j2,2]  + t3 * ek[j3,2]
        by4  = t1 * ek[j1,3]  + T2 * ek[j2,3]  + t3 * ek[j3,3]
        by5  = t1 * ek[j1,4]  + T2 * ek[j2,4]  + t3 * ek[j3,4]
        by6  = t1 * ek[j1,5]  + T2 * ek[j2,5]  + t3 * ek[j3,5]
        by7  = t1 * ek[j1,6]  + T2 * ek[j2,6]  + t3 * ek[j3,6]
        by8  = t1 * ek[j1,7]  + T2 * ek[j2,7]  + t3 * ek[j3,7]
        by9  = t1 * ek[j1,8]  + T2 * ek[j2,8]  + t3 * ek[j3,8]
        by10 = t1 * ek[j1,9]  + T2 * ek[j2,9]  + t3 * ek[j3,9]
        by11 = t1 * ek[j1,10] + T2 * ek[j2,10] + t3 * ek[j3,10]
        by12 = t1 * ek[j1,11] + T2 * ek[j2,11] + t3 * ek[j3,11]

        bz1  = r1 * ek[j1,0]  + r2 * ek[j2,0]  + r3 * ek[j3,0]
        bz2  = r1 * ek[j1,1]  + r2 * ek[j2,1]  + r3 * ek[j3,1]
        bz3  = r1 * ek[j1,2]  + r2 * ek[j2,2]  + r3 * ek[j3,2]
        bz4  = r1 * ek[j1,3]  + r2 * ek[j2,3]  + r3 * ek[j3,3]
        bz5  = r1 * ek[j1,4]  + r2 * ek[j2,4]  + r3 * ek[j3,4]
        bz6  = r1 * ek[j1,5]  + r2 * ek[j2,5]  + r3 * ek[j3,5]
        bz7  = r1 * ek[j1,6]  + r2 * ek[j2,6]  + r3 * ek[j3,6]
        bz8  = r1 * ek[j1,7]  + r2 * ek[j2,7]  + r3 * ek[j3,7]
        bz9  = r1 * ek[j1,8]  + r2 * ek[j2,8]  + r3 * ek[j3,8]
        bz10 = r1 * ek[j1,9]  + r2 * ek[j2,9]  + r3 * ek[j3,9]
        bz11 = r1 * ek[j1,10] + r2 * ek[j2,10] + r3 * ek[j3,10]
        bz12 = r1 * ek[j1,11] + r2 * ek[j2,11] + r3 * ek[j3,11]

        for i in range(12):
            ek[i,0]  = ek[i,0] \
                        + ek[i,j1] * bx1  + ek[i,j2] * by1  + ek[i,j3] * bz1
            ek[i,1]  = ek[i,1] \
                        + ek[i,j1] * bx2  + ek[i,j2] * by2  + ek[i,j3] * bz2
            ek[i,2]  = ek[i,2] \
                        + ek[i,j1] * bx3  + ek[i,j2] * by3  + ek[i,j3] * bz3
            ek[i,3]  = ek[i,3] \
                        + ek[i,j1] * bx4  + ek[i,j2] * by4  + ek[i,j3] * bz4
            ek[i,4]  = ek[i,4] \
                        + ek[i,j1] * bx5  + ek[i,j2] * by5  + ek[i,j3] * bz5
            ek[i,5]  = ek[i,5] \
                        + ek[i,j1] * bx6  + ek[i,j2] * by6  + ek[i,j3] * bz6
            ek[i,6]  = ek[i,6] \
                        + ek[i,j1] * bx7  + ek[i,j2] * by7  + ek[i,j3] * bz7
            ek[i,7]  = ek[i,7] \
                        + ek[i,j1] * bx8  + ek[i,j2] * by8  + ek[i,j3] * bz8
            ek[i,8]  = ek[i,8] \
                        + ek[i,j1] * bx9  + ek[i,j2] * by9  + ek[i,j3] * bz9
            ek[i,9]  = ek[i,9] \
                        + ek[i,j1] * bx10 + ek[i,j2] * by10 + ek[i,j3] * bz10
            ek[i,10] = ek[i,10] \
                        + ek[i,j1] * bx11 + ek[i,j2] * by11 + ek[i,j3] * bz11
            ek[i,11] = ek[i,11] \
                        + ek[i,j1] * bx12 + ek[i,j2] * by12 + ek[i,j3] * bz12

        return ek

    def elka5(self, ek, j):
        ####################################################
        ##     トラス（０）　ラーメン（１）の結合状態     ##
        ##        self.kTR1[i] - self.kTR2[j)             ##
        ##                    x,y,z - x,y,z               ##
        ##     ３次元   j=1  (0,0,1 - 0,0,1)              ##
        ##              j=2  (0,1,0 - 0,1,0)              ##
        ##              j=3  (1,0,0 - 1,0,0)              ##
        ##                ** ELkA5.for **                 ##
        ####################################################
        if j == 0:
            j1 = 3
            j2 = 4
            j3 = 9
            j4 = 10
        elif j == 1:
            j1 = 3
            j2 = 5
            j3 = 9
            j4 = 11
        elif j == 2:
            j1 = 4
            j2 = 5
            j3 = 10
            j4 = 11
        else:
            return


        s1 = 1 / (ek[j1,j1] * ek[j2,j2] - ek[j2,j1] * ek[j1,j2])

        t5               = (ek[j1,j2] * ek[j2,0]  - ek[j2,j2] * ek[j1,0])  * s1
        t6               = (ek[j1,j2] * ek[j2,1]  - ek[j2,j2] * ek[j1,1])  * s1
        t7               = (ek[j1,j2] * ek[j2,2]  - ek[j2,j2] * ek[j1,2])  * s1
        if j == 0: t8    = (ek[j1,j2] * ek[j2,5]  - ek[j2,j2] * ek[j1,5])  * s1
        if j == 1: t8    = (ek[j1,j2] * ek[j2,4]  - ek[j2,j2] * ek[j1,4])  * s1
        if j == 2: t8    = (ek[j1,j2] * ek[j2,3]  - ek[j2,j2] * ek[j1,3])  * s1
        t9               = (ek[j1,j2] * ek[j2,6]  - ek[j2,j2] * ek[j1,6])  * s1
        t10              = (ek[j1,j2] * ek[j2,7]  - ek[j2,j2] * ek[j1,7])  * s1
        t11              = (ek[j1,j2] * ek[j2,8]  - ek[j2,j2] * ek[j1,8])  * s1
        if j == 0: t12   = (ek[j1,j2] * ek[j2,11] - ek[j2,j2] * ek[j1,11]) * s1
        if j == 1: t12   = (ek[j1,j2] * ek[j2,10] - ek[j2,j2] * ek[j1,10]) * s1
        if j == 2: t12   = (ek[j1,j2] * ek[j2,9]  - ek[j2,j2] * ek[j1,9])  * s1

        v5               = (ek[j2,j1] * ek[j1,0]  - ek[j1,j1] * ek[j2,0])  * s1
        v6               = (ek[j2,j1] * ek[j1,1]  - ek[j1,j1] * ek[j2,1])  * s1
        v7               = (ek[j2,j1] * ek[j1,2]  - ek[j1,j1] * ek[j2,2])  * s1
        if j == 0: v8    = (ek[j2,j1] * ek[j1,5]  - ek[j1,j1] * ek[j2,5])  * s1
        if j == 1: v8    = (ek[j2,j1] * ek[j1,4]  - ek[j1,j1] * ek[j2,4])  * s1
        if j == 2: v8    = (ek[j2,j1] * ek[j1,3]  - ek[j1,j1] * ek[j2,3])  * s1
        v9               = (ek[j2,j1] * ek[j1,6]  - ek[j1,j1] * ek[j2,6])  * s1
        v10              = (ek[j2,j1] * ek[j1,7]  - ek[j1,j1] * ek[j2,7])  * s1
        v11              = (ek[j2,j1] * ek[j1,8]  - ek[j1,j1] * ek[j2,8])  * s1
        if j == 0: v12   = (ek[j2,j1] * ek[j1,11] - ek[j1,j1] * ek[j2,11]) * s1
        if j == 1: v12   = (ek[j2,j1] * ek[j1,10] - ek[j1,j1] * ek[j2,10]) * s1
        if j == 2: v12   = (ek[j2,j1] * ek[j1,9]  - ek[j1,j1] * ek[j2,9])  * s1


        smx5             = ek[j3,0]  + ek[j3,j1] * t5  + ek[j3,j2] * v5
        smx6             = ek[j3,1]  + ek[j3,j1] * t6  + ek[j3,j2] * v6
        smx7             = ek[j3,2]  + ek[j3,j1] * t7  + ek[j3,j2] * v7
        if j == 0: smx8  = ek[j3,5]  + ek[j3,j1] * t8  + ek[j3,j2] * v8
        if j == 1: smx8  = ek[j3,4]  + ek[j3,j1] * t8  + ek[j3,j2] * v8
        if j == 2: smx8  = ek[j3,3]  + ek[j3,j1] * t8  + ek[j3,j2] * v8
        smx9             = ek[j3,6]  + ek[j3,j1] * t9  + ek[j3,j2] * v9
        smx10            = ek[j3,7]  + ek[j3,j1] * t10 + ek[j3,j2] * v10
        smx11            = ek[j3,8]  + ek[j3,j1] * t11 + ek[j3,j2] * v11
        if j == 0: smx12 = ek[j3,11] + ek[j3,j1] * t12 + ek[j3,j2] * v12
        if j == 1: smx12 = ek[j3,10] + ek[j3,j1] * t12 + ek[j3,j2] * v12
        if j == 2: smx12 = ek[j3,9]  + ek[j3,j1] * t12 + ek[j3,j2] * v12

        smy5             = ek[j4,0]  + ek[j4,j1] * t5  + ek[j4,j2] * v5
        smy6             = ek[j4,1]  + ek[j4,j1] * t6  + ek[j4,j2] * v6
        smy7             = ek[j4,2]  + ek[j4,j1] * t7  + ek[j4,j2] * v7
        if j == 0: smy8  = ek[j4,5]  + ek[j4,j1] * t8  + ek[j4,j2] * v8
        if j == 1: smy8  = ek[j4,4]  + ek[j4,j1] * t8  + ek[j4,j2] * v8
        if j == 2: smy8  = ek[j4,3]  + ek[j4,j1] * t8  + ek[j4,j2] * v8
        smy9             = ek[j4,6]  + ek[j4,j1] * t9  + ek[j4,j2] * v9
        smy10            = ek[j4,7]  + ek[j4,j1] * t10 + ek[j4,j2] * v10
        smy11            = ek[j4,8]  + ek[j4,j1] * t11 + ek[j4,j2] * v11
        if j == 0: smy12 = ek[j4,11] + ek[j4,j1] * t12 + ek[j4,j2] * v12
        if j == 1: smy12 = ek[j4,10] + ek[j4,j1] * t12 + ek[j4,j2] * v12
        if j == 2: smy12 = ek[j4,9]  + ek[j4,j1] * t12 + ek[j4,j2] * v12

        c2   = (-ek[j2,j2] * ek[j1,j3] + ek[j1,j2] * ek[j2,j3]) * s1
        c3   = (-ek[j2,j2] * ek[j1,j4] + ek[j1,j2] * ek[j2,j4]) * s1
        d2   = ( ek[j2,j1] * ek[j1,j3] - ek[j1,j1] * ek[j2,j3]) * s1
        d3   = ( ek[j2,j1] * ek[j1,j4] - ek[j1,j1] * ek[j2,j4]) * s1
        sk11 =   ek[j3,j1] * c2 + ek[j3,j2] * d2 + ek[j3,j3]
        sk12 =   ek[j3,j1] * c3 + ek[j3,j2] * d3 + ek[j3,j4]
        sk21 =   ek[j4,j1] * c2 + ek[j4,j2] * d2 + ek[j4,j3]
        sk22 =   ek[j4,j1] * c3 + ek[j4,j2] * d3 + ek[j4,j4]
        skk  = 1 / (sk11 * sk22 - sk21 * sk12)

        x35  = (sk12 * smy5  - sk22 * smx5)  * skk
        x36  = (sk12 * smy6  - sk22 * smx6)  * skk
        x37  = (sk12 * smy7  - sk22 * smx7)  * skk
        x38  = (sk12 * smy8  - sk22 * smx8)  * skk
        x39  = (sk12 * smy9  - sk22 * smx9)  * skk
        x310 = (sk12 * smy10 - sk22 * smx10) * skk
        x311 = (sk12 * smy11 - sk22 * smx11) * skk
        x312 = (sk12 * smy12 - sk22 * smx12) * skk

        x45  = (sk21 * smx5  - sk11 * smy5)  * skk
        x46  = (sk21 * smx6  - sk11 * smy6)  * skk
        x47  = (sk21 * smx7  - sk11 * smy7)  * skk
        x48  = (sk21 * smx8  - sk11 * smy8)  * skk
        x49  = (sk21 * smx9  - sk11 * smy9)  * skk
        x410 = (sk21 * smx10 - sk11 * smy10) * skk
        x411 = (sk21 * smx11 - sk11 * smy11) * skk
        x412 = (sk21 * smx12 - sk11 * smy12) * skk

        x15  = t5  + c2 * x35  + c3 * x45
        x16  = t6  + c2 * x36  + c3 * x46
        x17  = t7  + c2 * x37  + c3 * x47
        x18  = t8  + c2 * x38  + c3 * x48
        x19  = t9  + c2 * x39  + c3 * x49
        x110 = t10 + c2 * x310 + c3 * x410
        x111 = t11 + c2 * x311 + c3 * x411
        x112 = t12 + c2 * x312 + c3 * x412

        x25  = v5  + d2 * x35  + d3 * x45
        x26  = v6  + d2 * x36  + d3 * x46
        x27  = v7  + d2 * x37  + d3 * x47
        x28  = v8  + d2 * x38  + d3 * x48
        x29  = v9  + d2 * x39  + d3 * x49
        x210 = v10 + d2 * x310 + d3 * x410
        x211 = v11 + d2 * x311 + d3 * x411
        x212 = v12 + d2 * x312 + d3 * x412

        for i in range(12):
            if not (j == 0 and i == 3) and  not (j == 0 and i == 4) and  \
               not (j == 1 and i == 3) and  not (j == 1 and i == 5) and  \
               not (j == 2 and i == 4) and  not (j == 2 and i == 5) and  \
               not (j == 0 and i == 9) and  not (j == 0 and i == 10) and \
               not (j == 1 and i == 9) and  not (j == 1 and i == 11) and \
               not (j == 2 and i == 10) and not (j == 2 and i == 11):

                ek1 = ek[i,0]  + ek[i,j1] * x15 \
                    + ek[i,j2] * x25 + ek[i,j3] * x35 + ek[i,j4] * x45
                ek2 = ek[i,1]  + ek[i,j1] * x16 \
                    + ek[i,j2] * x26 + ek[i,j3] * x36 + ek[i,j4] * x46
                ek3 = ek[i,2]  + ek[i,j1] * x17 \
                    + ek[i,j2] * x27 + ek[i,j3] * x37 + ek[i,j4] * x47

                if j == 0:
                    ek6 = ek[i,5]  + ek[i,j1] * x18 \
                        + ek[i,j2] * x28 + ek[i,j3] * x38 + ek[i,j4] * x48

                if j == 1:
                    ek5 = ek[i,4]  + ek[i,j1] * x18 \
                        + ek[i,j2] * x28 + ek[i,j3] * x38 + ek[i,j4] * x48

                if j == 2:
                    ek4 = ek[i,3]  + ek[i,j1] * x18 \
                        + ek[i,j2] * x28 + ek[i,j3] * x38 + ek[i,j4] * x48

                ek7 = ek[i,6]  + ek[i,j1] * x19 \
                    + ek[i,j2] * x29 + ek[i,j3] * x39 + ek[i,j4] * x49

                ek8 = ek[i,7]  + ek[i,j1] * x110 \
                    + ek[i,j2] * x210 + ek[i,j3] * x310 + ek[i,j4] * x410

                ek9 = ek[i,8]  + ek[i,j1] * x111 \
                    + ek[i,j2] * x211 + ek[i,j3] * x311 + ek[i,j4] * x411

                if j == 0:
                    ek12 = ek[i,11] + ek[i,j1] * x112 \
                         + ek[i,j2] * x212 + ek[i,j3] * x312 + ek[i,j4] * x412

                if j == 1:
                    ek11 = ek[i,10] + ek[i,j1] * x112 \
                         + ek[i,j2] * x212 + ek[i,j3] * x312 + ek[i,j4] * x412

                if j == 2:
                    ek10 = ek[i,9]  + ek[i,j1] * x112 \
                         + ek[i,j2] * x212 + ek[i,j3] * x312 + ek[i,j4] * x412


                ek[i,0]  = ek1
                ek[i,1]  = ek2
                ek[i,2]  = ek3
                ek[i,3]  = ek4
                ek[i,4]  = ek5
                ek[i,5]  = ek6
                ek[i,6]  = ek7
                ek[i,7]  = ek8
                ek[i,8]  = ek9
                ek[i,9]  = ek10
                ek[i,10] = ek11
                ek[i,11] = ek12


        for i in range(12):
            ek[j1,i] = 0
            ek[j2,i] = 0
            ek[j3,i] = 0
            ek[j4,i] = 0
            ek[i,j1] = 0
            ek[i,j2] = 0
            ek[i,j3] = 0
            ek[i,j4] = 0

        return ek

    def get_world_vector(self, vector):
        return self.tm.get_world_vector(vector)

    def get_member_vector(self, vector):
        return self.tm.get_member_vector(vector)

    # member に載荷している分布荷重を節点荷重に変換する
    def Convert_memberload_to_nodeload(self):
        tt    = self.tt

        wfe   = np.dot(tt.T, self.WBUNPU_3DFRM())    # Thermal load vector in global coordinate

        ### ピン結合を考慮した分布荷重を換算 #####
        # x軸回り曲げモーメント（ねじりモーメント）
        if self.k1x ==0: # i端ピン結合
            self.set_pin(wfe, 3)
        if self.k2x ==0: # j端ピン結合
            self.set_pin(wfe, 9)
        # y軸回り曲げモーメント
        if self.k1y ==0: # i端ピン結合
            self.set_pin(wfe, 4)
        if self.k2y ==0: # j端ピン結合
            self.set_pin(wfe, 10)
        # z軸回り曲げモーメント
        if self.k1z ==0: # i端ピン結合
            self.set_pin(wfe, 5)
        if self.k2z ==0: # j端ピン結合
            self.set_pin(wfe, 11)

        # Thermal load vector in global coordinate
        tfe   = np.dot(tt.T, self.TFVEC_3DFRM())

        return tfe, wfe


    def WBUNPU_3DFRM(self):
        wfe_l = np.zeros(12,dtype=np.float64)
        el    = self.el

        # 軸力（x軸方向力）
        if self.tx == 0:
            wfe_l[0]  =  el/6 * (2 * self.wxi + self.wxj)
            wfe_l[6]  =  el/6 * (self.wxi + 2 * self.wxj)
        else:
            wfe_l[0] ,wfe_l[6]  =  self.set_wx(self.tx, self.EA, el, self.wxi, self.wxj)

        # x軸回りモーメント（ねじりモーメント）
        if self.tr == 0:
            wfe_l[3]  =  el/6 * (2 * self.wti + self.wtj)
            wfe_l[9]  =  el/6 * (self.wti + 2 * self.wtj)
        else:
            wfe_l[3] ,wfe_l[9]  =  self.set_wx(self.tr, self.GJ, el, self.wti, self.wtj)

        if self.ty == 0:
            # y軸方向せん断力
            wfe_l[1]  =  el/20 * (7 * self.wyi + 3 * self.wyj)
            wfe_l[7]  =  el/20 * (3 * self.wyi + 7 * self.wyj)
            # z軸回り曲げモーメント
            wfe_l[5]  = (el**2)/60 * (3 * self.wyi + 2 * self.wyj)
            wfe_l[11] = -(el**2)/60 * (2 * self.wyi + 3 * self.wyj)
        else:
            wfe_l[1], wfe_l[7], wfe_l[5], wfe_l[11] =  self.set_wy(self.ty, self.EIz, el, self.wyi, self.wyj)
            wfe_l[11] = -wfe_l[11]

        if self.tz == 0:
            # z軸方向せん断力
            wfe_l[2]  =  el/20 * (7 * self.wzi + 3 * self.wzj)
            wfe_l[8]  =  el/20 * (3 * self.wzi + 7 * self.wzj)
            # y軸回り曲げモーメント
            wfe_l[4]  = -(el**2)/60 * (3 * self.wzi + 2 * self.wzj)
            wfe_l[10] = (el**2)/60 * (2 * self.wzi + 3 * self.wzj)
        else:
            wfe_l[2], wfe_l[8], wfe_l[4], wfe_l[10] =  self.set_wy(self.tz, self.EIy, el, self.wzi, self.wzj)
            wfe_l[4] = -wfe_l[4]

        return wfe_l


    def TFVEC_3DFRM(self):
        # Thermal load vector  in local coordinate system
        tfe_l    = np.zeros(12,dtype=np.float64)
        E        = self.ee    # elastic modulus
        AA       = self.aa    # section area
        alpha    = self.alpha # thermal expansion coefficient
        tempe    = self.wte   # 0.5*(self.inp.deltaT[i]+self.inp.deltaT[j])
        tfe_l[0] = -E * AA * alpha * tempe
        tfe_l[6] =  E * AA * alpha * tempe
        return tfe_l


    def set_wx(self, K2, EA, L1, P8, P9):
        # 軸方向のバネが付いた部材の分布荷重

        D2 = np.sqrt(K2 / EA)
        E0 = np.exp(L1 * D2)
        Q4 = (E0 + 1 / E0) / 2
        Q5 = (E0 - 1 / E0) / 2

        W6 = P8
        W7 = P9

        n = 19
        D = 0
        if L1 > 6:
            n = 19
        elif L1 > 4:
            n = 15
        elif L1 > 3:
            n = 9
        elif L1 > 2:
            n = 7
        elif L1 > 1.2:
            n = 5
        else:
            n = 3
            D = L1 / (n+1)

        for k in range(n):
            i = k + 1
            j = i + 1
            X =  j * D
            P7 = (P9 - P8) / L1 * X + P8
            E0 = np.exp((L1 - X) * D2)
            Q6 = (E0 + 1 / E0) / 2
            Q7 = (E0 - 1 / E0) / 2
            W1 = 1 / Q5 * Q7 * P7
            W2 = (Q6 - Q4 / Q5 * Q7) * P7

            if j % 2 == 0:
                W6 = W6 + 2 * W1
                W7 = W7 + 2 * W2
            else:
                W6 = W6 + 4 * W1
                W7 = W7 + 4 * W2

        W6 = W6 * D / 3
        W7 = W7 * D / 3

        return W6, W7


    def set_wy(self, K1, EI, L1, P8, P9):
        # 軸直角方向のバネが付いた部材の分布荷重
        if P8==0 and P9==0:
            return 0, 0, 0, 0

        D1 = (K1 / (4 * EI))**0.25
        W1 = L1 * D1
        E0 = np.exp(W1)
        C4 = (E0 + 1 / E0) / 2
        S4 = (E0 - 1 / E0) / 2
        C5 = np.cos(W1)
        S5 = np.sin(W1)
        Q = S4 * S4 - S5 * S5

        W6 = P8
        W7 = P9
        W8 = 0
        W9 = 0

        n = 19
        if L1 > 6:
            n = 19
        elif L1 > 4:
            n = 15
        elif L1 > 3:
            n = 9
        elif L1 > 2:
            n = 7
        elif L1 > 1.2:
            n = 5
        else:
            n = 3
        D = L1 / (n+1)

        for k in range(n):
            i = k + 1
            j = i + 1
            X =  i * D
            P7 = (P9 - P8) / L1 * X + P8
            W1 = (L1 - X) * D1
            E0 = np.exp(W1)
            C6 = (E0 + 1 / E0) / 2
            S6 = (E0 - 1 / E0) / 2
            C7 = np.cos(W1)
            S7 = np.sin(W1)
            W1 = X * D1
            E0 = np.exp(W1)
            C8 = (E0 + 1 / E0) / 2
            S8 = (E0 - 1 / E0) / 2
            C9 = np.cos(W1)
            S9 = np.sin(W1)
            W1 = (S4 * (C6 * S9 + S6 * C9) - S5 * (C7 * S8 + S7 * C8)) * P7 / Q
            W2 = (S4 * (C8 * S7 + S8 * C7) - S5 * (C9 * S6 + S9 * C6)) * P7 / Q
            W3 = (S4 * S6 * S9 - S5 * S7 * S8) * P7 / (Q * D1)
            W4 = (S4 * S8 * S7 - S5 * S9 * S6) * P7 / (Q * D1)

            if j % 2 == 0:
                W6 = W6 + 4 * W1
                W7 = W7 + 4 * W2
                W8 = W8 + 4 * W3
                W9 = W9 + 4 * W4
            else:
                W6 = W6 + 2 * W1
                W7 = W7 + 2 * W2
                W8 = W8 + 2 * W3
                W9 = W9 + 2 * W4

        W6 = W6 * D / 3
        W7 = W7 * D / 3
        W8 = W8 * D / 3
        W9 = W9 * D / 3

        return W6, W7, W8, W9


    def set_pin(self, wfe_l, n):
        # ピン結合を考慮した分布荷重の換算関数

        ARW = copy.copy(self.gk)
        ARB1 = ARW[n,:]
        B3 = ARB1[n]
        if B3 == 0:
            # 荷重がない
            return wfe_l

        ARB1[n] = 0

        ARF = copy.copy(wfe_l)
        ARF3 = ARF[n]
        ARF[n] = 0

        ART1 = np.zeros([n+1],dtype=np.float64)
        ART2 = np.zeros([n+1],dtype=np.float64)

        for i in range(n+1):
            ART1[i] = ARB1[i] / B3
            ART2[i] = ART1[i] * ARF3
            wfe_l[i] = ARF[i] - ART2[i]

        return wfe_l


    # 要素荷重を整理し解析用変数にセットする
    def set_load(self, inp):
        fe = np.zeros([9], dtype=np.float64)

        for lm in inp.load_member:
            if lm['m'] != self.ID:
                continue
            if  lm['mark'] == 2:
                self.set_load_mk2(fe, lm)
            elif  lm['mark'] == 9:
                self.set_load_mk9(fe, lm)

        return fe

    # 温度荷重を整理し解析用変数にセットする
    def set_load_mk9(self, fe: np.ndarray, target_load):
        P1 = target_load['P1']
        fe[8] += P1 # Temperature change of element
        return


    # 分布荷重を整理し解析用変数にセットする
    def set_load_mk2(self, fe: np.ndarray, target_load):

        load_direct = target_load['direction']
        Pi = target_load['P1']
        Pj = target_load['P2']

        if load_direct == "x" :
            fe[0] += Pi # wxi
            fe[1] += Pj # wxj
        elif load_direct == "y":
            fe[2] += Pi # wyi
            fe[3] += Pj # wyj
        elif load_direct == "z":
            fe[4] += Pi # wzi
            fe[5] += Pj # wzj
        elif load_direct == "r":
            fe[6] += Pi # wri
            fe[7] += Pj # wrj
        else:
            if load_direct == "gx":
                Pii = self.get_member_vector([Pi,0,0])
                Pjj = self.get_member_vector([Pj,0,0])
            elif load_direct == "gy":
                Pii = self.get_member_vector([0,Pi,0])
                Pjj = self.get_member_vector([0,Pj,0])
            elif load_direct == "gz":
                Pii = self.get_member_vector([0,0,Pi])
                Pjj = self.get_member_vector([0,0,Pj])

            fe[0] += Pii[0] # wxi
            fe[1] += Pjj[0] # wxj
            fe[2] += Pii[1] # wyi
            fe[3] += Pjj[1] # wyj
            fe[4] += Pii[2] # wzi
            fe[5] += Pjj[2] # wzj

        return


    # 変位量から断面力を計算する
    def getFsec(self, disg, Gindex):

        i = Gindex[self.IDi]
        j = Gindex[self.IDj]
        work = np.array([
            disg[i],
            disg[i + 1],
            disg[i + 2],
            disg[i + 3],
            disg[i + 4],
            disg[i + 5],
            disg[j],
            disg[j + 1],
            disg[j + 2],
            disg[j + 3],
            disg[j + 4],
            disg[j + 5]
        ])

        #全体座標系の断面力 --------------------------------------------
        fseg = np.dot(self.ck, work)

        # 分布荷重を換算した分を減じる
        for i in range(self.nod * self.nfree):
            fseg[i] -= self.wfe[i]
            # 温度荷重を換算した分を減じる
            fseg[i] -= self.tfe[i]

        #要素座標系の断面力 ---------------------------------------------
        fsec = np.dot(self.tt, fseg)

        # 構造力学の定義
        for i in [0, 3, 5, 7, 8, 10]:
            fsec[i] *= -1


        return fsec, fseg


class tMatrix:

    def __init__(self, x1, y1, z1, x2, y2, z2, theta):

        xx = x2-x1
        yy = y2-y1
        zz = z2-z1
        el = np.sqrt(xx**2+yy**2+zz**2)

        tt      = np.zeros([12,12],dtype=np.float64) # transformation matrix
        t1      = np.zeros([3,3],dtype=np.float64)
        t2      = np.zeros([3,3],dtype=np.float64)
        theta   = np.radians(theta)             # chord angle
        t1[0,0] = 1
        t1[1,1] = np.cos(theta)
        t1[1,2] = np.sin(theta)
        t1[2,1] =-np.sin(theta)
        t1[2,2] = np.cos(theta)
        ll=(x2-x1)/el
        mm=(y2-y1)/el
        nn=(z2-z1)/el

        if x2-x1==0.0 and y2-y1==0.0:
            t2[0,2]=nn
            t2[1,0]=nn
            t2[2,1]=1.0
        else:
            qq=np.sqrt(ll**2+mm**2)
            t2[0,0]=ll
            t2[0,1]=mm
            t2[0,2]=nn
            t2[1,0]=-mm/qq
            t2[1,1]= ll/qq
            t2[2,0]=-ll*nn/qq
            t2[2,1]=-mm*nn/qq
            t2[2,2]=qq
        t3=np.dot(t1,t2)
        tt[ 0:3,  0:3] = t3[0:3, 0:3]
        tt[ 3:6,  3:6] = t3[0:3, 0:3]
        tt[ 6:9,  6:9] = t3[0:3, 0:3]
        tt[9:12, 9:12] = t3[0:3, 0:3]

        self.t3 = t3
        self.tt = tt


    def get_world_vector(self, vector):
        result =  np.dot(self.t3.T, vector)
        return result


    def get_member_vector(self, vector):
        result =  np.dot(self.t3, vector)
        return result

