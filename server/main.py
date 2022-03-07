import json
from calcrate import FrameCalc

def main(request):

    # Set CORS headers for the preflight request
    if request.method == 'OPTIONS':
        headers = {
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Methods': 'GET, POST',
            'Access-Control-Allow-Headers': 'Origin, X-Requested-With, Content-Type, Accept, Content-Encoding',
            'Access-Control-Max-Age': '3600'
        }

        return ('', 204, headers)

    # Set CORS headers for the main request
    headers = {
        'Access-Control-Allow-Origin': '*',
        'Access-Control-Allow-Headers': 'Origin, X-Requested-With, Content-Type, Accept, Content-Encoding'
    }

    if request.method == 'GET':
        return (json.dumps({ 'body': 'Hello World!'}), 200, headers) # テスト用コード

    # 計算開始
    # 入力データを読み込む
    calc = FrameCalc(js=request.get_json())

    # 計算開始
    error, result = calc.calcrate()
    if error != None:
        return (json.dump({'error': error}), 500, headers)


    # 結果を返送する
    response = json.dumps(result)
    return (response, 200, headers)





