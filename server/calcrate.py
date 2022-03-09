
import math
import json


class FrameCalc:

    def __init__(self, js=None, filename=None):
        if js == None:
            with open(filename) as f:
                s = f.read()
            js = json.loads(s, object_pairs_hook=dict)

        # ユーザー情報
        self.production = js['production'] if 'production' in js else True
        self.uid = js['uid'] if 'uid' in js else None
        self.username = js['username'] if 'username' in js else ''
        self.password = js['password'] if 'password' in js else ''
        # 点データ
        self.node = js["node"]
        # 地盤情報
        self.strana = js["strana"]
        # 水位情報
        self.waterlevel = js["waterlevel"]
        # 荷重情報
        self.load = js["load"]
        # 設計条件
        self.init = js["init"]


    # 全ての荷重ケースの解析を開始する
    def calcrate(self):
        # 入力を計算用に変換
        # 地盤条件
        plygons = []
        for k in self.strana:
            st = self.strana[k]
            tmp = []
            for i in st['organization']:
                pos = self.node[ i - 1 ]
                tmp.append(pos)
            poly = Polygon(tmp)
            plygons.append({
                'organization': poly,
                'strana': k,
            })
        # 水位線
        segments = LineString(self.waterlevel)

        # 斜面の安定計算
        """将来は、以下をトライアル計算することになる"""
        # 円弧の中心位置と半径を仮定
        center = Point(17, 12)   # 円弧中心
        radius = 8 # 半径
        # 円を生成
        circle = center.buffer(radius)
        # 円の内側の地層、水位線を求める
        plygons2 = []
        # 地盤条件をトリム
        for p in plygons:
            poly = p['organization']
            if circle.intersects(poly):
                plygons2.append({
                    'organization': circle.intersection(poly),
                    'strana': p['strana'],
                })
        # 水位線をトリム
        if circle.intersects(segments):
            segments2 = circle.intersection(segments)
        else:
            segments2 = None

        # 水位で地層を分割する
        plygons3 = []
        for p3 in plygons2:
            poly = p3['organization']
            for p in split(poly, segments).geoms: # 水位線で地層を分割する
                # 水位線と地層の最も近い点
                p1, _ = nearest_points(segments, p)
                _, y1 = p1.xy
                # 地層の重心位置
                p2 = p.centroid
                _, y2 = p2.xy

                plygons3.append({
                    'organization': p,
                    'strana': p3['strana'],
                    'wet': y2[0] < y1[0]  # 水位より地層のy座標が小さかったら True
                })


        # 円と地表、地層、水位折れ線との交点を求める
        # 交点
        gx = []
        gy = []

        # ポリゴンを１つずつ
        for p3 in plygons3:
            poly = p3['organization']
            x, y = poly.exterior.xy
            gx.extend(x)
            gy.extend(y)

        # 図化する
        cross = []
        for i in range(len(gx)):
            cross.append((gx[i],gy[i]))

        # 荷重の載荷点
        for lo in self.load:
            gx.append(lo['x_s'])
            gx.append(lo['x_d'])

        # 重複する点を削除
        splitX = list(set(gx))

        splitX.sort()

        # 地層をsplitX(交点のx座標)で分割する
        cx, cy = center.xy
        top = cy[0] + radius * 2
        bottom = cy[0] - radius * 2

        plygons4 = []
        delX = [] # splitX の中から不要な点を削除するリスト
        for i in range(1, len(splitX)):
            x1 = splitX[i-1]
            x2 = splitX[i]
            splitter = Polygon([(x1, top),(x2, top),(x2, bottom),(x1, bottom), (x1, top)])
            # ポリゴンを１つずつ
            p4 = []
            for p3 in plygons3:
                poly = p3['organization']
                if splitter.intersects(poly):
                    p = splitter.intersection(poly)
                    if type(p) == Polygon:
                        x, y = p.exterior.xy
                        p4.append({
                            'organization': p,
                            'strana': p3['strana'],
                            'wet': p3['wet'],
                            'top': max(y), # 最も高い点
                            'bottom': min(y) # 最も低い点
                        })
            if len(p4) > 0:
                plygons4.append(p4)
            else:
                delX.append(i)

        # splitX の中から不要な点を削除する
        for i in reversed(delX):
            splitX.pop(i)


        # 円を分割する
        # 分割位置と円の交点
        Edge = []
        x, y = center.xy
        for i in range(1, len(splitX)):
            x1 = splitX[i-1]
            y1 = cli(x1, x[0], y[0], radius)
            x2 = splitX[i]
            y2 = cli(x2, x[0], y[0], radius)

            Edge.append([
                LineString([(x1, y1[0]),(x2, y2[0])]),
                LineString([(x1, y1[1]),(x2, y2[1])])
            ])


        # 安定計算
        # 区間ごとの地層の情報を集計
        block_info = []
        for i in range(len(plygons4)):
            info = { "weight" : 0, "edge": [] } # 区間ごとの地層情報
            for p in plygons4[i]:
                # 土質条件
                s = self.strana[p['strana']]
                # 単位体積重量
                g = s['gamma']
                if p['wet']:
                    g = g - 10 # 水中重量
                # ブロック面積
                poly = p['organization']
                info["weight"] += poly.area * g
                # 層境界
                for edge in Edge[i]:
                    mid = midPoint(edge) # 境界の中点Y座標
                    top = p['top']
                    bottom = p['bottom']

                    if bottom < mid and mid < top:
                        # 境界層だった場合
                        a = bowArea(radius, edge.length) # 弦の面積
                        info["weight"] += a * g
                        info["edge"].append({
                            'c' :s['c'],     # 底面の地層の粘着力
                            'fai': s['fai'], # 底面の地層の内部摩擦角
                            'rad': getAngle(edge), # 弦の角度
                            'len': edge.length  # 底面の長さ
                        })

        block_info.append(info)

        # 安定計算
        P = 0 # 発生力
        R = 0 # 抵抗力
        for info in block_info:
            w = info['weight']
            for edge in info['edge']:
                c = edge['c']
                fai = math.radians(edge['fai'])
                rad = edge['rad']
                l = edge['len']

                e = math.tan(fai)
                P += w * math.sin(rad)
                R += e * w * math.cos(rad) + c * l

        P = abs(P)
        A = abs(R) / P # 安全率

        return {
            '安全率': A,
            '判定': 'OK' if A>1 else 'NG',
            '発生力': P,
            '抵抗力': R
        }

