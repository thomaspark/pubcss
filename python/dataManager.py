
import numpy as np
from collections import OrderedDict
import math

from dataInput import dataInput
from matrixBar import tMatrix
import copy

class dataManager:

    def __init__(self, _inp):

        self.SplitChar1 = 'n'
        self.SplitChar2 = 'l'
        self.base = _inp
        self.price = len(self.base.node) # 課金のポイントサイズ

        # 要素分割処理を行う
        ## 一時的にローカル変数とする
        self.node = self.base.node                      ### 節点情報をセット
        self.member = self.base.member                  ### 要素情報を編集
        self.notice_points = self.base.notice_points[:] ### 着目点をセット

        ## 分割後の部材名を格納する変数を用意する
        for ID in self.member:
            m = self.member[ID]
            m['child'] = []

        ### 着目点で要素を分割する
        for nc in self.notice_points:
            targetID = nc['m']
            baseMember = self.member[targetID]
            points = sorted(list(set(nc['Points'])))
            currentPos = 0
            for index in range(len(points)):
                pos = points[index] - currentPos
                newID = nc['m'] + self.SplitChar1 + str(index + 1)
                # 荷重以外の分割をする
                targetMember, newMember = self.SplitMember1of3Base(targetID, pos, newID)
                # 要素バネを分割する
                for type in self.base.fix_member:
                    target_fix_member = self.base.fix_member[type]
                    self.SplitMember2of3FixMember(targetID, newID, target_fix_member)
                # 分布荷重を分割する
                for load_ID in self.base.load:
                    target_load = self.base.load[load_ID]
                    self.SplitMember3of3Load(targetMember, newMember, target_load['load_member'], target_load['load_node'])
                # 着目点を更新する
                currentPos += targetMember['L']
                # 親部材に子部材を登録する
                baseMember['child'].append(newID)
                targetID = newID

        return


    # 全荷重IDのリスト
    @property
    def loadList(self):
        return list(self.base.load)


    # 分布荷重など変換処理を行ない解析できる状態にする
    def setTargetData(self, CaseID: str):

        target_load = self.base.load[CaseID]

        ## 節点情報をセット
        self.node = copy.deepcopy(self.base.node)

        ## 材料情報をセット
        self.element = []
        eNo = target_load['element']
        if eNo in self.base.element:
            self.element = copy.deepcopy(self.base.element[eNo])
        else:
            return 'element no not found'

        ## 要素情報をセット
        self.member = copy.deepcopy(self.base.member)
        # シェル情報の生合成をチェック
        for ID in self.member:
            m = self.member[ID]
            # 存在しない材料をもつシェルはエラー
            if not m['e'] in self.element:
                return 'have a member that contains element that do not exist.'

        ## シェル情報をセット
        self.shell = copy.deepcopy(self.base.shell)
        # シェル情報の生合成をチェック
        for ID in self.shell:
            sh = self.shell[ID]
            # 存在しない材料をもつシェルはエラー
            if not sh['e'] in self.element:
                return 'have a shell that contains element that do not exist.'

        ## 支点情報をセット
        self.fix_node = []
        fix_node_no = target_load['fix_node']
        if fix_node_no in self.base.fix_node:
            self.fix_node = copy.deepcopy(self.base.fix_node[fix_node_no])

        ## バネ情報をセット
        self.fix_member = []
        fix_member_no = target_load['fix_member']
        if fix_member_no in self.base.fix_member:
            self.fix_member = copy.deepcopy(self.base.fix_member[fix_member_no])

        ## 結合情報をセット
        joint_no = target_load['joint']
        self.joint = []
        if joint_no in self.base.joint:
            self.joint = copy.deepcopy(self.base.joint[joint_no])
            for jo in self.joint:
                m = self.member[jo['m']]
                m['xi'] = jo['xi']
                m['yi'] = jo['yi']
                m['zi'] = jo['zi']
                m['xj'] = jo['xj']
                m['yj'] = jo['yj']
                m['zj'] = jo['zj']

        # 分布荷重
        self.load_member = []
        if 'load_member' in target_load:
            self.load_member = target_load['load_member']

        # 節点荷重
        self.load_node = []
        if 'load_node' in target_load:
            self.load_node = target_load['load_node']

        # 分布荷重の端点 で部材を分割する
        split_point = dict()
        for lc in self.load_member:

            if lc['mark'] == 9:
                continue

            L1 = lc['L1']
            L2 = lc['L2']
            targetID = lc['m']
            m = self.member[targetID]
            L = m['L']
            if lc['mark'] == 2:
                L2 = L - L2 # L2 を j端からの距離から i端からの距離に変換している

            if 0 < L1 and L1 < L:
                if targetID in split_point:
                    split_point[targetID].append(L1)
                else:
                    split_point[targetID] = [L1]
            if 0 < L2 and L2 < L:
                if targetID in split_point:
                    split_point[targetID].append(L2)
                else:
                    split_point[targetID] = [L2]
            if lc['mark'] == 1 or lc['mark'] == 11:        
                if 0 == L1:
                    ln = {'m': lc['m'],'n': m['ni'], 'mark':lc['mark'], 'Pd': lc['P1'], 'direction': lc['direction'] }
                    self.load_node.append(self.changeNodeLoad(ln))
                    lc['P1'] = 0
                if L1 == L:
                    ln = {'m': lc['m'],'n': m['nj'], 'mark':lc['mark'], 'Pd': lc['P1'], 'direction': lc['direction'] }
                    self.load_node.append(self.changeNodeLoad(ln))
                    lc['P1'] = 0
                if 0 == L2:
                    ln = {'m': lc['m'],'n': m['ni'], 'mark':lc['mark'], 'Pd': lc['P2'], 'direction': lc['direction'] }
                    self.load_node.append(self.changeNodeLoad(ln))
                    lc['P2'] = 0
                if L2 == L:
                    ln = {'m': lc['m'],'n': m['nj'], 'mark':lc['mark'], 'Pd': lc['P2'], 'direction': lc['direction'] }
                    self.load_node.append(self.changeNodeLoad(ln))
                    lc['P2'] = 0

        if len(split_point) > 0:

            # 着目点で要素を分割する
            for ID in split_point:
                targetID = ID
                baseID = targetID.split(self.SplitChar1)[0].split(self.SplitChar2)[0]
                baseMember = self.member[baseID]
                points = sorted(list(set(split_point[targetID])))
                currentPos = 0
                for index in range(len(points)):
                    pos = points[index] - currentPos
                    newID = ID + self.SplitChar2 + str(index + 1)
                    # 荷重以外の分割をする
                    targetMember, newMember = self.SplitMember1of3Base(targetID, pos, newID)
                    # 要素バネを分割する
                    self.SplitMember2of3FixMember(targetID, newID, self.fix_member)
                    # 分布荷重を分割する
                    self.SplitMember3of3Load(targetMember, newMember, self.load_member, self.load_node)
                    # 着目点を更新する
                    currentPos += targetMember['L']
                    # 親部材に子部材を登録する
                    baseMember['child'].append(newID)
                    targetID = newID


        return None


    # 着目点で荷重を分割する
    def SplitMember3of3Load(self, targetMember, newMember, load_member, load_node):

        targetID = targetMember['ID']
        newID = newMember['ID']

        Li = targetMember['L']  # 分割位置
        Lj = newMember['L']     # 残りの長さ

        delete_list = []
        append_list = []

        for lc in [lm for lm in load_member if targetID == lm['m']]:

            # 既存の部材荷重の入力を直す
            L1 = lc['L1']
            L2 = lc['L2']
            P1 = lc['P1']
            P2 = lc['P2']
            mk = lc['mark']
            if mk == 2:
                if L1 >= Li:
                    ## new_member にのみ載荷する荷重
                    lc['m'] = newID ## 荷重を付け替える
                    lc['L1'] = L1 - Li
                elif L2 >= Lj:
                    ## target_member にのみ載荷する荷重
                    lc['L2'] = L2 - Lj
                else:
                    ## target_memberとnew_member にまたがって載荷する荷重
                    new_load = copy.deepcopy(lc)
                    # 既存の分布荷重の中間値を計算する

                    LL = Li + Lj
                    P3 = P1 + ( P2 - P1 ) / ( LL - L1 - L2 ) * ( Li - L1 )
                    ## target_member の荷重を修正する
                    lc['L2'] = 0
                    lc['P2'] = P3
                    ## new_member の荷重を追加する
                    new_load['m'] = newID
                    new_load['L1'] = 0
                    new_load['P1'] = P3
                    append_list.append(new_load)# 追加登録


            elif mk == 9:
                # 温度荷重
                new_load = copy.deepcopy(lc)
                new_load['m'] = newID
                append_list.append(new_load)# 追加登録
            else:
                # 集中荷重
                for i in ['1', '2']:
                    Px = lc['P'+i]
                    Lx = lc['L'+i]
                    if Px == 0:
                        continue
                    if math.isclose(Li, Lx):
                        ## もし節点と同じ位置に荷重があったら節点荷重として登録する
                        ln = {'m': newID,'n': newID, 'mark':mk, 'Pd': Px, 'direction': lc['direction'] }
                        load_node.append(self.changeNodeLoad(ln))
                        lc['P'+i] = 0
                    elif Li <= Lx:
                        lc['L'+i ] = max(Lx - Li, 0)

                if lc['P1'] == 0 and lc['P2'] == 0:
                    # P1, P2 両方とも節点荷重に変換してしまったら要素荷重から削除
                    delete_list.append(lc)

                elif Li < L1 and Li < L2:
                    ## new_member にのみ載荷する荷重
                    lc['m'] = newID

                else:
                    ## target_memberとnew_member にまたがって載荷する荷重
                    new_load = copy.deepcopy(lc)
                    if Li < L1:
                        ## P1はtarget_member, P2はnew_member
                        new_load['P1'] = 0
                        new_load['L1'] = 0
                        lc['m'] = newID
                        lc['P2'] = 0
                        lc['L2'] = 0
                    elif Li < L2:
                        ## P2はtarget_member, P1はnew_member
                        new_load['P2'] = 0
                        new_load['L2'] = 0
                        lc['m'] = newID
                        lc['P1'] = 0
                        lc['L1'] = 0
                    else:
                        new_load['P1'] = 0
                        new_load['L1'] = 0
                        new_load['P2'] = 0
                        new_load['L2'] = 0

                    if new_load['P1'] != 0 or new_load['P2'] != 0:
                        append_list.append(new_load)# 追加登録

        # 節点荷重に変換した荷重を削除
        for lc in delete_list:
            load_member.remove(lc)

        # 追加する荷重
        for lc in append_list:
            load_member.append(lc)

        return


    # 着目点で要素バネを分割する
    def SplitMember2of3FixMember(self, targetID, newID, fix_member):
        for fm in fix_member[:]:
            if targetID != fm['m']:
                continue
            # 新しい要素バネを追加
            newFm = copy.copy(fm)
            newFm['m'] = newID
            fix_member.append(newFm)
        return


    # 着目点で要素バネ, 荷重以外を分割する
    def SplitMember1of3Base(self, targetID: str, pos: int, newID: str):

        targetMember = self.member[targetID]

        # 新しいポイントを追加 --------------------------------------------------
        self.node[newID] =  self.GetMidPoint(targetMember, pos/1000)

        # 新しい部材を追加 -----------------------------------------------------
        newMember = {
            'ID': newID,
            'ni': newID,
            'nj': targetMember['nj'],
            'e': targetMember['e'],
            'xi': 1,
            'yi': 1,
            'zi': 1,
            'xj': targetMember['xj'],
            'yj': targetMember['yj'],
            'zj': targetMember['zj'],
            'cg': targetMember['cg']
        }
        # 部材長さを計算しておく
        newMember['L'] = round(self.GetLength(newMember) *1000)
        # 登録
        self.member[newID] = newMember

        # もともとの部材の情報を更新
        targetMember['ID'] = targetID
        # j端の接続情報を更新
        targetMember['nj'] = newID
        # 結合情報を更新
        targetMember['xj'] = 1
        targetMember['yj'] = 1
        targetMember['zj'] = 1
        # 部材長さを計算しておく
        targetMember['L'] = round(self.GetLength(targetMember)*1000)

        return targetMember, newMember


    # 部材の長さを計算する
    def GetLengthFromID(self, ID: str):
        m = self.member[ID]
        return self.GetLength(m)

    # 部材の長さを計算する
    def GetLength(self, target: dict):
        IDi   = target['ni']
        IDj   = target['nj']
        return self.GetDistanceFromID(IDi, IDj)

    # 2点間の距離を計算する
    def GetDistanceFromID(self, IDi: str, IDj: str):
        pi = self.node[IDi]
        pj = self.node[IDj]
        return self.GetDistance(pi, pj)

    # 2点間の距離を計算する
    def GetDistance(self, pi: dict, pj: dict):
        return self.base.GetDistance(pi, pj)

    # 部材途中の点の座標を計算する
    def GetMidPoint(self, target: dict, L1: float):
        pi = self.node[target['ni']]
        pj = self.node[target['nj']]

        x1 = pi['x']
        y1 = pi['y']
        z1 = pi['z']

        x2 = pj['x']
        y2 = pj['y']
        z2 = pj['z']
        L = self.GetDistance(pi, pj)

        n = L1 / L
        return {
            'x': x1 + (x2 - x1) * n,
            'y': y1 + (y2 - y1) * n,
            'z': z1 + (z2 - z1) * n
        }

    # 部材荷重を節点荷重に変換する
    def changeNodeLoad(self, load_member):

        target_member = self.member[load_member['m']]
        new_load = {
            'n': load_member['n'],
            'tx': 0, 'ty': 0, 'tz': 0,
            'rx': 0, 'ry': 0, 'rz': 0,
            'dx': 0, 'dy': 0, 'dz': 0,
            'ax': 0, 'ay': 0, 'az': 0
        }

        load_direct = load_member['direction']
        P1 = load_member['Pd']
        mk = load_member['mark']

        if load_direct == "gx":
            new_load['tx'] = P1 if mk==1 else 0
            new_load['rx'] = P1 if mk==11 else 0

        elif load_direct == "gy":
            new_load['ty'] = P1 if mk==1 else 0
            new_load['ry'] = P1 if mk==11 else 0

        elif load_direct == "gz":
            new_load['tz'] = P1 if mk==1 else 0
            new_load['rz'] = P1 if mk==11 else 0

        else:
            # 部材の角度を計算する
            pi = self.node[target_member['ni']]
            pj = self.node[target_member['nj']]
            t = tMatrix(pi['x'],pi['y'],pi['z'],pj['x'],pj['y'],pj['z'], target_member['cg'])

            # 荷重を絶対座標に換算する
            if load_direct == "x":
                re = t.get_world_vector([P1, 0, 0])
            elif load_direct == "y":
                re = t.get_world_vector([0, P1, 0])
            elif load_direct == "z":
                re = t.get_world_vector([0, 0, P1])

            # 節点荷重を登録
            if mk==1:
                new_load['tx'] = re[0]
                new_load['ty'] = re[1]
                new_load['tz'] = re[2]
            if mk==11:
                new_load['rx'] = re[0]
                new_load['ry'] = re[1]
                new_load['rz'] = re[2]

        return new_load


    #######################################################################
    # 解析結果の集計
    #######################################################################
    def getResultData(self, gmat, disg: dict, reac: dict, fsec: dict, shell_fsec: dict):
        result = dict()
        result['disg'] = self.getDisgJson(gmat, disg)
        result['reac'] = self.getReacJson(gmat, reac)
        result['fsec'] = self.getForceJson(gmat, fsec) # Bar要素の断面力
        result['shell_fsec'] = self.getShellForceJson(shell_fsec) # シェルの断面力   
        result['size'] = self.price # 計算ポイントの計算に用いる

        return result


    def getDisgJson(self, gmat, disg):
        # 変位
        dict_disg = OrderedDict()
        # オリジナルの入力データにある節点のみ返す
        for ID in self.node:
            iz = gmat.Gindex[ID]
            dict_disg[ID] = {
                "dx":disg[iz+0],
                "dy":disg[iz+1],
                "dz":disg[iz+2],
                "rx":disg[iz+3],
                "ry":disg[iz+4],
                "rz":disg[iz+5]
                }
        return dict_disg

    def getReacJson(self, gmat, reac):
        #反力
        dict_reac = OrderedDict()
        # オリジナルの入力データにある節点のみ返す
        for fn in self.fix_node:
            ID = fn['n']
            iz = gmat.Gindex[ID]
            dict_reac[ID] = {
                "tx":reac[iz+0],
                "ty":reac[iz+1],
                "tz":reac[iz+2],
                "mx":reac[iz+3],
                "my":reac[iz+4],
                "mz":reac[iz+5]
                }

        # 支点に節点荷重がある場合は荷重値を反力に加算する
        for fn in self.fix_node:
            ID = fn['n']
            # 支点と同じ点にある荷重を探す
            Lo = None
            for ln in self.load_node:
                if ID == ln['n']:
                    Lo = ln
                    break
            if Lo == None:
                continue
            # 反力に節点荷重を加える
            re = dict_reac[ID]
            if fn['tx'] != 0:
                re['tx'] += Lo['tx']
            if fn['ty'] != 0:
                re['ty'] += Lo['ty']
            if fn['tz'] != 0:
                re['tz'] += Lo['tz']
            if fn['rx'] != 0:
                re['mx'] += Lo['rx']
            if fn['ry'] != 0:
                re['my'] += Lo['ry']
            if fn['rz'] != 0:
                re['mz'] += Lo['rz']

        return dict_reac

    def getForceJson(self, gmat, fsec: dict):

        dict_fsec = OrderedDict() # 親：断面力

        for org_id in self.base.member:

            if self.SplitChar1 in org_id:
                continue
            if self.SplitChar2 in org_id:
                continue

            base_member = self.base.member[org_id]

            i = 1
            dict1 = OrderedDict()

            # 親部材の断面力を集計する
            dict2 = dict()
            Px = 'P' + str(i)
            f = fsec[org_id]
            dict2['fxi'] = f[0]
            dict2['fyi'] = f[1]
            dict2['fzi'] = f[2]
            dict2['mxi'] = f[3]
            dict2['myi'] = f[4]
            dict2['mzi'] = f[5]

            def_member = self.member[org_id]
            for child_id in def_member['child']:
                child_member = self.member[child_id]
                if base_member['nj'] == child_member['nj']:
                    f = fsec[child_id]
                    break

            dict2['fxj'] = f[ 6]
            dict2['fyj'] = f[ 7]
            dict2['fzj'] = f[ 8]
            dict2['mxj'] = f[ 9]
            dict2['myj'] = f[10]
            dict2['mzj'] = f[11]
            dict2['L'] = self.GetDistanceFromID(base_member['ni'], base_member['nj'])

            # dict1[Px] = dict2

            # 着目点 child の断面力を集計する
            # for child1_id in base_member['child']:
            #     base_member = self.base.member[child1_id]

            #     i += 1
            #     dict2 = dict()
            #     Px = 'P' + str(i)
            #     f = fsec[child1_id]
            #     dict2['fxi'] = f[0]
            #     dict2['fyi'] = f[1]
            #     dict2['fzi'] = f[2]
            #     dict2['mxi'] = f[3]
            #     dict2['myi'] = f[4]
            #     dict2['mzi'] = f[5]

            #     for child2_id in def_member['child']:
            #         child2_member = self.member[child2_id]
            #         if base_member['nj'] == child2_member['nj']:
            #             f = fsec[child2_id]
            #             break

            #     dict2['fxj'] = f[ 6]
            #     dict2['fyj'] = f[ 7]
            #     dict2['fzj'] = f[ 8]
            #     dict2['mxj'] = f[ 9]
            #     dict2['myj'] = f[10]
            #     dict2['mzj'] = f[11]
            #     dict2['L'] = self.GetDistanceFromID(base_member['ni'], base_member['nj'])

            #     dict1[Px] = dict2

            dict_fsec[org_id] = dict2

        return dict_fsec


    def getShellForceJson(self, fsec):

        dict_fsec = OrderedDict() # 親：断面力

        for org_id in fsec:

            dict2 = dict()

            # 親部材の断面力を集計する
            f = fsec[org_id]

            if org_id in dict_fsec:
                f2 = dict_fsec[org_id]
            else:
                f2 = np.zeros(12)

            dict2['fxi'] = f[ 0] + f2[ 0]
            dict2['fyi'] = f[ 1] + f2[ 1]
            dict2['fzi'] = f[ 2] + f2[ 2]
            dict2['mxi'] = f[ 3] + f2[ 3]
            dict2['myi'] = f[ 4] + f2[ 4]
            dict2['mzi'] = f[ 5] + f2[ 5]

            dict2['fxj'] = f[ 6] + f2[ 6]
            dict2['fyj'] = f[ 7] + f2[ 7]
            dict2['fzj'] = f[ 8] + f2[ 8]
            dict2['mxj'] = f[ 9] + f2[ 9]
            dict2['myj'] = f[10] + f2[10]
            dict2['mzj'] = f[11] + f2[11]
            
            dict_fsec[org_id] = dict2

        return dict_fsec
        