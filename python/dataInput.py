import numpy as np
import json
from collections import OrderedDict
import copy

# compress
import gzip
import base64

class dataInput:

    def __init__(self):
        self.isCompress = False

    def setCompressJSON(self, data):
        self.isCompress = True
        js = self.decompress(data)
        self.setJSONString(js)

    def readTestData(self, fnameR='test.json'):
        f = open(fnameR, encoding="utf-8")
        fstr = f.read()  # ファイル終端まで全て読んだデータを返す
        f.close()
        self.setJSONString(fstr)

    def setJSONString(self, fstr):
        js = json.loads(fstr, object_pairs_hook=dict)
        self.setJSON(js)

    def setJSON(self, js: dict):

        self.production = js['production'] if 'production' in js else True
        self.uid = js['uid'] if 'uid' in js else None
        self.username = js['username'] if 'username' in js else ''
        self.password = js['password'] if 'password' in js else ''

        # ケースによって変わらないパラメータを取得

        ## 節点
        self.node = OrderedDict()
        # 型をfloat に 統一する
        a = js['node']
        for ID in a:
            n = a[ID]
            n['x'] = float(n['x'])
            n['y'] = float(n['y'])
            n['z'] = float(n['z'])
            self.node[str(ID)] = n


        ## 材料
        self.element = []
        if 'element' in js:
            self.element = js['element']
            for no in self.element:
                element = self.element[no]
                for en in element:
                    el = element[en]
                    el['E']  = float(el['E'])                       # elastic modulus
                    el['G']  = float(el['G'])                       # self.ee/2/(1+Poisson's ratio)
                    el['A']  = float(el['A'])                       # bar: section area, shell: thickness
                    el['J']  = float(el['J']) if 'J' in el else 0   # tortional constant
                    el['Iy'] = float(el['Iy']) if 'Iy' in el else 0 # moment of inertia around y-axis
                    el['Iz'] = float(el['Iz']) if 'Iz' in el else 0 # moment of inertia around z-axis
                    el['Xp'] = float(el['Xp']) if 'Xp' in el else 0 # alpha : thermal expansion coefficient

        ## 部材
        self.member = []
        if 'member' in js:
            self.member = js['member']
            # 型をstring に 統一する, 部材長さを計算しておく
            for ID in self.member:
                m = self.member[ID]
                m['ni'] = str(m['ni'])
                m['nj'] = str(m['nj'])
                m['e'] = str(m['e'])
                m['xi'] = 1
                m['yi'] = 1
                m['zi'] = 1
                m['xj'] = 1
                m['yj'] = 1
                m['zj'] = 1
                m['cg'] = m['cg'] if 'cg' in m else 0
                m['L'] =  round(self.GetLength(m) * 1000) #部材長さを計算しておく
                # 節点に関する整合性をチェック
                for n in [m['ni'], m['nj'] ]:
                    # 存在しない節点を含むシェルはエラー
                    if not n in self.node:
                        return 'have a shell that contains nodes that do not exist.'


        ## シェル要素
        self.shell = []
        if 'shell' in js:
            self.shell = js['shell']
            # 型をstring に 統一する
            for ID in self.shell:
                sh = self.shell[ID]
                sh['e'] = str(sh['e'])
                # 節点に関する整合性をチェック
                nods = []
                for n in sh['nodes']:
                    n = str(n)
                    # 存在しない節点を含むシェルはエラー
                    if not n in self.node:
                        return 'have a shell that contains nodes that do not exist.'
                    nods.append(n)
                sh['nodes'] = nods



        ## 着目点
        self.notice_points = []
        if 'notice_points' in js:
            self.notice_points = js['notice_points']

            for nc in self.notice_points[:]:
                # 型をstring に 統一する
                nc['m'] = str(nc['m'])
                # 存在しない要素の着目点は削除
                if not nc['m'] in self.member:
                    self.notice_points.remove(nc)
                    continue

            # 着目点を1000倍する
            for nc in self.notice_points:
                points = nc['Points']
                newlist = []
                for pos in points:
                    newlist.append(round(pos * 1000))
                nc['Points'] = newlist

            # 着目点に関する整合性をチェック
            for nc in self.notice_points:
                points = nc['Points']
                for pos in points[:]:
                    if pos == 0:
                        # 距離ゼロの着目点は削除
                        points.remove(pos)
                        continue
                    m = self.member[nc['m']]
                    if pos >= m['L']:
                        # 部材長と同じ距離の着目点は削除
                        points.remove(pos)
                        continue

            # 同じ部材に複数の着目点入力があったら結合する
            Flag = True
            index = 0
            while Flag:
                Flag = False
                for i in range(index, len(self.notice_points)):
                    nc = self.notice_points[i]
                    for j in range(index, len(self.notice_points)):
                        if i == j:
                            continue
                        nc2 = self.notice_points[j]
                        if nc['m'] == nc2['m']:
                            nc['Points'].extend(nc2['Points'])
                            self.notice_points.remove(nc2)
                            Flag = True
                            break
                    if Flag == True:
                        index = i
                        break


        ## 支点
        self.fix_node = []
        if 'fix_node' in js:
            self.fix_node = js['fix_node']

            for fix_node in self.fix_node.values():
                # 型をstring に 統一する
                for fn in fix_node:
                    fn['n'] = str(fn['n'])
                    fn['tx'] = float(fn['tx']) if 'tx' in fn else 0  #fixed in x-direction
                    fn['ty'] = float(fn['ty']) if 'ty' in fn else 0  #fixed in y-direction
                    fn['tz'] = float(fn['tz']) if 'tz' in fn else 0  #fixed in z-direction
                    fn['rx'] = float(fn['rx']) if 'rx' in fn else 0  #fixed in rotation around x-axis
                    fn['ry'] = float(fn['ry']) if 'ry' in fn else 0  #fixed in rotation around y-axis
                    fn['rz'] = float(fn['rz']) if 'rz' in fn else 0  #fixed in rotation around z-axis
                # 存在しない節点の支点は削除
                for fn in fix_node[:]:
                    if not fn['n'] in self.node:
                        fix_node.remove(fn)

        ## バネ
        self.fix_member = []
        if 'fix_member' in js:
            self.fix_member = js['fix_member']
            for fix_member in self.fix_member.values():
                # 型をstring に 統一する
                for fm in fix_member:
                    fm['m'] = str(fm['m'])
                    fm['tx'] = float(fm['tx']) if 'tx' in fm else 0
                    fm['ty'] = float(fm['ty']) if 'ty' in fm else 0
                    fm['tz'] = float(fm['tz']) if 'tz' in fm else 0
                    fm['tr'] = float(fm['tr']) if 'tr' in fm else 0
                # 存在しない要素のバネは削除 ...しなくてもいい使わないから
                # for fm in fix_member[:]:
                #     if not fm['m'] in self.member:
                #         fix_member.remove(fm)

        ## 結合
        self.joint = []
        if 'joint' in js:
            self.joint = js['joint']
            for joint in self.joint.values():
                # 型をstring に 統一する
                for jo in joint:
                    jo['m'] = str(jo['m'])
                    jo['xi'] = int(jo['xi']) if 'xi' in jo else 1
                    jo['yi'] = int(jo['yi']) if 'yi' in jo else 1
                    jo['zi'] = int(jo['zi']) if 'zi' in jo else 1
                    jo['xj'] = int(jo['xj']) if 'xj' in jo else 1
                    jo['yj'] = int(jo['yj']) if 'yj' in jo else 1
                    jo['zj'] = int(jo['zj']) if 'zj' in jo else 1
                # 存在しない要素の結合は削除 ...しなくてもいい使わないから
                # for jo in joint[:]:
                #     if not jo['m'] in self.member:
                #         joint.remove(jo)

        ## 荷重
        self.load = []
        if 'load' in js:
            self.load = js['load']

            ### 型を 統一する
            for load in self.load.values():
                load['fix_node']   = str(load['fix_node'])   if 'fix_node'   in load else '1'
                load['fix_member'] = str(load['fix_member']) if 'fix_member' in load else '1'
                load['element']    = str(load['element'])    if 'element'    in load else '1'
                load['joint']      = str(load['joint'])      if 'joint'      in load else '1'

                if 'load_node' in load:
                    load_node =load['load_node']
                    # 型をstring に 統一する
                    for ln in load_node:
                        ln['n'] = str(ln['n'])
                        # 荷重
                        ln['tx'] = float(ln['tx']) if 'tx' in ln else 0  # load in x-direction
                        ln['ty'] = float(ln['ty']) if 'ty' in ln else 0  # load in y-direction
                        ln['tz'] = float(ln['tz']) if 'tz' in ln else 0  # load in z-direction
                        ln['rx'] = float(ln['rx']) if 'rx' in ln else 0  # moment around x-axis
                        ln['ry'] = float(ln['ry']) if 'ry' in ln else 0  # moment around y-axis
                        ln['rz'] = float(ln['rz']) if 'rz' in ln else 0  # moment around z-axis
                        # 強制変位
                        ln['dx'] = float(ln['dx']) if 'dx' in ln else 0  # deload in x-direction
                        ln['dy'] = float(ln['dy']) if 'dy' in ln else 0  # load in y-direction
                        ln['dz'] = float(ln['dz']) if 'dz' in ln else 0  # load in z-direction
                        ln['ax'] = float(ln['ax']) if 'ax' in ln else 0  # moment around x-axis
                        ln['ay'] = float(ln['ay']) if 'ay' in ln else 0  # moment around y-axis
                        ln['az'] = float(ln['az']) if 'az' in ln else 0  # moment around z-axis

                    # 存在しない節点の荷重は削除
                    for ln in load_node[:]:
                        if not ln['n'] in self.node:
                            load_node.remove(ln)
                else:
                    load['load_node'] = []

                if 'load_member' in load:
                    load_member = load['load_member']
                    # 型を統一する
                    for lm in load_member:
                        lm['m'] = str(lm['m'])
                        lm['mark'] = int(lm['mark'])
                        lm['L1'] = round(lm['L1']*1000) if 'L1' in lm else 0 # 1000倍して int型で管理
                        lm['L2'] = round(lm['L2']*1000) if 'L2' in lm else 0 # 1000倍して int型で管理
                        lm['P1'] = float(lm['P1']) if 'P1' in lm else 0
                        lm['P2'] = float(lm['P2']) if 'P2' in lm else 0

                    # 存在しない要素の荷重は削除
                    for lm in load_member[:]:
                        if not lm['m'] in self.member:
                            load_member.remove(lm)
                else:
                    load['load_member'] = []


        return None

    # 部材の長さを計算する
    def GetLength(self, target: dict):
        IDi   = target['ni']
        IDj   = target['nj']
        pi = self.node[IDi]
        pj = self.node[IDj]
        return self.GetDistance(pi, pj)

    # 2点間の距離を計算する
    def GetDistance(self, pi: dict, pj: dict):
        xx = pj['x']-pi['x']
        yy = pj['y']-pi['y']
        zz = pj['z']-pi['z']
        return np.sqrt(xx**2+yy**2+zz**2)

    def compress(self, js: str):
        # gzip圧縮する
        l = gzip.compress(js.encode())
        # Base64エンコードする
        byteBase64 = base64.b64encode(l)
        # string に変換
        return byteBase64.decode()

    def decompress(self, byteBase64):
        # base64型 を もとに戻す
        b = base64.b64decode(byteBase64)
        # str型に変換し、カンマでばらして int 配列に変換する
        l = eval(b) #[int(n) for n in b.decode().split(',')]
        # gzipを解凍する
        return gzip.decompress(bytes(l))


