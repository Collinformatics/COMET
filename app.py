from flask import (Flask, jsonify, render_template, request,
                   send_file, send_from_directory, abort)
from flask_wtf.csrf import CSRFProtect, generate_csrf
from functions import WebApp
import io
import os
import sys
import threading
import time
import zipfile


# Set up the app
app = Flask(__name__)
csrf = CSRFProtect(app)

# Set app limits
app.config['MAX_FORM_MEMORY_SIZE'] = 16 * 1024 * 1024  # 16 MB (or higher)
app.config['MAX_FORM_PARTS'] = 10000  # Max number of fields in a multipart/form-data body
# print(f'Max content length: {app.config['MAX_CONTENT_LENGTH']}\n'
#       f'Max number of fields: {app.config['MAX_FORM_PARTS']:,}\n'
#       f'Max form memory: {app.config["MAX_FORM_MEMORY_SIZE"]:,}\n')

# Initialize: Application
webapp: WebApp = WebApp()
webapp.getKey(app)

# Figure storage
figures = {}


def parseForm():
    data = {}

    # Parse form
    keys = set(request.form.keys())
    for key in keys:
        values = request.form.getlist(key)
        data[key] = values
        if len(values) == 1:
            data[key] = values[0]
        else:
            data[key] = values

    # Parse form
    for key in request.files.keys():
        files = request.files.getlist(key)
        buffer = []
        for value in files:
            if value:
                buf = io.BytesIO(value.read())
                buf.filename = value.filename # type: ignore[attr-defined]
                print(f'Filename: {buf.filename}')
                buffer.append(buf)
        data[key] = buffer if len(buffer) > 1 else buffer[0]

    return data


def getValue(json, key, default):
    val = json.get(key)
    if val is None or val == "":
        return default
    return float(val)


@app.route('/')
def home():
    webapp.jobDone = True
    # return render_template('home.html')
    w = 'home.html'
    x = 'processDNA.html'
    y = 'combineProfiles.html'
    z = 'prediction.html'
    return render_template(
        w, csrf_token=generate_csrf()
    )


@app.route('/processDNA')
def pProcessDNA():
    webapp.jobDone = True
    return render_template(
        'processDNA.html', csrf_token=generate_csrf()
    )


@app.route('/filterAA')
def pFilterAA():
    webapp.jobDone = True
    return render_template(
        'filterAA.html', csrf_token=generate_csrf()
    )


@app.route('/filterMotif')
def pFilterMotif():
    webapp.jobDone = True
    return render_template(
        'filterMotif.html', csrf_token=generate_csrf()
    )


@app.route('/combineProfiles')
def pCombineProfiles():
    webapp.jobDone = True
    return render_template(
        'combineProfiles.html', csrf_token=generate_csrf()
    )


@app.route('/predictions')
def pPredictions():
    webapp.jobDone = True
    return render_template(
        'prediction.html', csrf_token=generate_csrf()
    )


@app.route('/evalFormDNA', methods=['POST'])
def evalDNA():
    keys = ['fileExp', 'fileExpRev', 'fileBg', 'fileBgRev']
    for key in keys:
        files = request.files.getlist(key)
        files = [f for f in files if f and f.filename]
    thread = threading.Thread(target=webapp.evalDNA,
                              args=(parseForm(),))
    thread.start()
    time.sleep(1)
    return render_template(
        'results.html', parameters=webapp.jobParams
    )


@app.route('/evalFormFilterAA', methods=['POST'])
def filterSubs():
    thread = threading.Thread(target=webapp.evalData,
                              args=(parseForm(),))
    thread.start()
    time.sleep(1)
    return render_template(
        'results.html', parameters=webapp.jobParams
    )


@app.route('/evalFormFilterMotif', methods=['POST'])
def filterMotif():
    thread = threading.Thread(target=webapp.evalData,
                              args=(parseForm(),True,))
    thread.start()
    time.sleep(2)
    return render_template(
        'setEntropy.html', parameters=webapp.jobParams,
        minS=webapp.minS, motifPos=list(webapp.motifPos.items())
    )


@app.route('/evalFormCombineProfiles', methods=['POST'])
def combineProfiles():
    thread = threading.Thread(target=webapp.evalData,
                              args=(parseForm(),False,True,))
    thread.start()
    time.sleep(1)
    return render_template(
        'combineProfiles.html', parameters=webapp.jobParams,
        motifPos=list(webapp.motifPos.items())
    )


@app.route('/evalFormPredActivity', methods=['POST'])
def predActivity():
    thread = threading.Thread(target=webapp.evalData,
                              args=(parseForm(), False , False, True,))
    thread.start()
    time.sleep(1)
    return render_template(
        'results.html', parameters=webapp.jobParams,
    )


@app.route(f'/<filename>')
def getFigure(filename):
    response = send_from_directory(webapp.pathFigs, filename)
    response.headers['Cache-Control'] = 'no-store'
    return response


@app.route('/checkFigures')
def checkFigures():
    return jsonify(webapp.figures)


@app.route('/refreshCSRF', methods=['GET'])
def refreshCSRF():
    return jsonify({'csrf_token': generate_csrf()})


@app.route('/jobSummary')
def jobSummary():
    print('Job Summary')
    return render_template(
        'results.html', parameters=webapp.jobParams()
    )


@app.route('/results')
def results():
    return render_template(
        'results.html', figures=webapp.figures,
        parameters=webapp.jobParams, csrf_token=generate_csrf()
    )


@app.route('/download', methods=['POST'])
@csrf.exempt
def download():
    try:
        # Define file tag
        if webapp.datasetTagMotif:
            tag = webapp.datasetTagMotif.replace(' ', '_')
        else:
            tag = webapp.datasetTag.replace(' ', '_')
        dir = f'{webapp.enzymeName}-{tag}'

        print(f'Downloading Directory: {webapp.pathDir}')
        memoryFile = io.BytesIO()
        with zipfile.ZipFile(memoryFile, 'w', zipfile.ZIP_DEFLATED) as zf:
            for root, dirs, files in os.walk(webapp.pathDir):
                for file in files:
                    filePath = os.path.join(root, file)
                    arcName = os.path.relpath(filePath, webapp.pathDir)
                    # print(f'File: {filePath}\n* arcName: {arcName}')
                    arcName = os.path.join(dir, arcName)
                    # print(f'* arcName: {arcName}\n')
                    zf.write(filePath, arcName)
        memoryFile.seek(0)
        # zip_bytes = memoryFile.getvalue()  # Extract raw bytes

        return send_file(
            memoryFile,
            mimetype='application/zip',
            as_attachment=True,
            download_name=f'comet-{dir}.zip'
        )
    except Exception as e:
        webapp.logError(f'ERROR: download()\n\n{e}')


@app.route('/updateFig', methods=['POST'])
def updateFig():
    json = request.get_json()
    webapp.minS = float(json.get('minS'))
    webapp.jobParams['Minimum ∆S'] = webapp.minS
    webapp.minES = float(json.get('minES'))
    webapp.minESRel = float(json.get('minESRel'))
    webapp.selectMotifPos()
    webapp.plotEntropy()
    data = {
        'status': 'success',
        'minS': webapp.figures.get('minS'),
        'motifPos': list(webapp.motifPos.items())
    }
    return jsonify(data)


@app.route('/setEntropy', methods=['GET'])
def setEntropy():
    return render_template(
        'setEntropy.html', minS=webapp.minS,
        minES=webapp.minES, minESRel=webapp.minESRel,
        parameters=webapp.jobParams, motifPos=list(webapp.motifPos.items())
    )


@app.route('/updateMinS', methods=['POST'])
def updateMinS():
    json = request.get_json()
    if not json:
        return jsonify({'error': 'No data'}), 400

    webapp.minS = getValue(json, 'minS', webapp.minS)
    webapp.minS = float(json.get('minS'))
    webapp.jobParams['Minimum ∆S'] = webapp.minS
    webapp.jobParams['Minimum ES Filter'] = getValue(json, 'minES',
                                                     webapp.minES)
    webapp.jobParams['Minimum ES Release'] = getValue(json, 'minESRel',
                                                      webapp.minESRel)
    webapp.selectMotifPos()
    return jsonify({
        'motifPos': list(webapp.motifPos.items()),
        'entropy': webapp.figures['entropy']
    })


@app.route('/motifPos')
def motifPos():
    return jsonify(list(webapp.motifPos.items()))


@app.route('/comet', methods=['POST'])
def comet():
    thread = threading.Thread(target=webapp.comet,
                              args=(parseForm(),))
    thread.start()
    time.sleep(1)
    return render_template(
        'results.html', parameters=webapp.jobParams,
        motifPos=list(webapp.motifPos.items())
    )


@app.route('/jobStatus')
def jobStatus():
    # print(f'Job Done: {webapp.jobDone}')
    return {'jobStatus': webapp.jobDone}


@app.route('/error')
def error():
    return render_template('error.html')


@app.route('/downloadFiles', methods=['POST'])
def downloadFiles():
    try:
        dir = 'TemplateData'
        print(f'Downloading Directory: {dir}')
        memoryFile = io.BytesIO()
        with zipfile.ZipFile(memoryFile, 'w', zipfile.ZIP_DEFLATED) as zf:
            for root, dirs, files in os.walk(dir):
                for file in files:
                    filePath = os.path.join(root, file)
                    arcName = os.path.relpath(filePath, dir)
                    # print(f'File: {filePath}\n* arcName: {arcName}')
                    arcName = os.path.join(dir, arcName)
                    # print(f'* arcName: {arcName}\n')
                    zf.write(filePath, arcName)
        memoryFile.seek(0)
        # zip_bytes = memoryFile.getvalue()  # Extract raw bytes

        return send_file(
            memoryFile,
            mimetype='application/zip',
            as_attachment=True,
            download_name=f'comet-{dir}.zip'
        )
    except Exception as e:
        webapp.logError(f'ERROR: download()\n\n{e}')


# Run the app
if __name__ == '__main__':
    app.run(threaded=True, debug=False, use_reloader=False, port=9090)
