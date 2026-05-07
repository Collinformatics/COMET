from flask import (Flask, jsonify, render_template, request,
                   send_file, send_from_directory)
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

# Initialize: Application
webapp = WebApp()
webapp.getKey(app)

# Figure storage
figures = {}


def parseForm():
    data = {}

    # Parse form
    keys = set(request.form.keys())
    # print('Data:')
    for key in keys:
        values = request.form.getlist(key)
        data[key] = values
        if len(values) == 1:  # Check if there are multiple values
            data[key] = values[0]
        else:
            data[key] = values  # Store as a list if multiple values exist
        # print('*', key, data[key])

    # Parse form
    for key in request.files.keys():
        files = request.files.getlist(key)
        bufs = []
        for value in files:
            if value:
                buf = io.BytesIO(value.read())
                buf.filename = value.filename
                bufs.append(buf)
        if bufs:
            data[key] = bufs if len(bufs) > 1 else bufs[0]

    return data


@app.route('/run', methods=['POST'])
def run():
    data = request.get_json() # Get the JSON body
    message = data.get('message') # Extract the message

    # Call method
    data = webapp.pressButton(message)

    return jsonify(data)


@app.route('/')
def home():
    # return render_template('home.html')
    return render_template(
        'filterMotif.html', csrf_token=generate_csrf()
    )


@app.route('/processDNA')
def pProcessDNA():
    return render_template(
        'processDNA.html', csrf_token=generate_csrf()
    )


@app.route('/filterAA')
def pFilterAA():
    return render_template(
        'filterAA.html', csrf_token=generate_csrf()
    )


@app.route('/filterMotif')
def pFilterMotif():
    return render_template(
        'filterMotif.html', csrf_token=generate_csrf()
    )


@app.route('/combineProfiles')
def pCombineProfiles():
    return render_template(
        'combineProfiles.html', csrf_token=generate_csrf()
    )


@app.route('/resources')
def resources():
    return render_template('resources.html',)


@app.route('/evalFormDNA', methods=['POST'])
def evalDNA():
    # Process the dset
    # webapp.evalDNA(parseForm())
    # print('Job Done: Eval DNA')
    thread = threading.Thread(target=webapp.evalDNA, args=(parseForm(),))
    thread.start()
    time.sleep(2)
    return render_template(
        'results.html', parameters=webapp.jobParams
    )


@app.route('/evalFormFilterAA', methods=['POST'])
def filterSubs():
    # Parse the form
    thread = threading.Thread(target=webapp.evalSubs, args=(parseForm(),))
    thread.start()
    time.sleep(2)
    return render_template(
        'results.html', parameters=webapp.jobParams
    )


@app.route('/evalFormFilterMotif', methods=['POST'])
def filterMotif():
    # Parse the form
    thread = threading.Thread(target=webapp.evalSubs, args=(parseForm(),True,))
    thread.start()
    time.sleep(2)
    return render_template(
        'setEntropy.html', parameters=webapp.jobParams,
        minS=webapp.minS, motifPos=list(webapp.motifPos.items())
    )


@app.route('/evalProfiles', methods=['POST'])
def combineProfiles():
    print('Combine Profiles')
    thread = threading.Thread(target=webapp.combineProfiles, args=(parseForm(),))
    thread.start()
    time.sleep(2)
    return render_template(
        'combineProfiles.html', parameters=webapp.jobParams,
        motifPos=list(webapp.motifPos.items())
    )


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
        parameters=webapp.jobParams
    )


@app.route(f'/<filename>')
def getFigure(filename):
    response = send_from_directory(webapp.pathFigs, filename)
    response.headers['Cache-Control'] = 'no-store'
    return response


@app.route('/checkFigures')
def checkFigures():
    return jsonify(webapp.figures)


"""
Archive:  comet.zip
  End-of-central-directory signature not found.  Either this file is not
  a zipfile, or it constitutes one disk of a multi-part archive.  In the
  latter case the central directory and zipfile comment will be found on
  the last disk(s) of this archive.
unzip:  cannot find zipfile directory in one of comet.zip or
        comet.zip.zip, and cannot find comet.zip.ZIP, period.
"""

@app.route('/download', methods=['POST'])
def download():
    print(f'Downloading Directory: {webapp.pathDir}')
    memory_file = io.BytesIO()
    with zipfile.ZipFile(memory_file, 'w', zipfile.ZIP_DEFLATED) as zf:
        for root, dirs, files in os.walk(webapp.pathDir):
            for file in files:
                filePath = os.path.join(root, file)
                arcname = os.path.relpath(filePath, webapp.pathDir)
                print(f'File: {filePath}\narcname: {arcname}')
                zf.write(filePath, arcname)
    memory_file.seek(0)

    # Define file tag
    if webapp.datasetTagMotif:
        tag = webapp.datasetTagMotif.replace(' ', '_')
    else:
        tag = webapp.datasetTag.replace(' ', '_')

    return send_file(
        memory_file,
        mimetype='application/zip',
        as_attachment=True,
        download_name=f'{webapp.enzymeName}-{tag}.zip'
    )


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
    webapp.minS = float(json.get('minS'))
    webapp.selectMotifPos()
    webapp.figures['entropy'] = webapp.plotEntropy()
    return jsonify({
        'motifPos': list(webapp.motifPos.items()),
        'entropy': webapp.figures['entropy']
    })


@app.route('/motifPos')
def motifPos():
    return jsonify(list(webapp.motifPos.items()))


@app.route('/comet', methods=['POST'])
def comet():
    webapp.jobDone = False
    print('COMET')
    thread = threading.Thread(target=webapp.filterMotifs, args=(parseForm(),))
    thread.start()
    # time.sleep(2)
    return render_template(
        'results.html', parameters=webapp.jobParams,
        motifPos=list(webapp.motifPos.items())
    )


@app.route('/jobStatus')
def jobStatus():
    print(f'Job Done: {webapp.jobDone}')
    return {'jobStatus': webapp.jobDone}


# Run the app
if __name__ == '__main__':
    app.run(threaded=True, debug=True, use_reloader=False, port=9090)
    # sys.exit(0)
