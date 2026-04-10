from flask import (Flask, jsonify, render_template, request,
                   send_file, send_from_directory)
from flask_wtf.csrf import CSRFProtect, generate_csrf
from functions import WebApp
import io
import os
import sys
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
        if len(values) == 1:  # Check if there are multiple values
            data[key] = values[0]
        else:
            data[key] = values  # Store as a list if multiple values exist
        # print('*', key, data[key])
    # Parse files
    for key, value in request.files.items():
        if value:
            data[key] = value

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
    return render_template('filterMotif.html',
                           csrf_token=generate_csrf())


@app.route('/processDNA')
def pProcessDNA():
    return render_template('processDNA.html',
                           csrf_token=generate_csrf())


@app.route('/filterAA')
def pFilterAA():
    return render_template('filterAA.html',
                           csrf_token=generate_csrf())


@app.route('/filterMotif')
def pFilterMotif():
    return render_template('filterMotif.html',
                           csrf_token=generate_csrf())


@app.route('/resources')
def resources():
    return render_template('resources.html',)


@app.route('/results')
def results():
    return render_template('results.html',
                           figures=webapp.figures,
                           parameters=webapp.jobParams)


@app.route(f'/<filename>')
def getFigure(filename):
    return send_from_directory(webapp.pathFigs, filename)


@app.route('/checkFigures')
def checkFigures():
    if webapp.done:
        return jsonify(webapp.figures)
    else:
        return jsonify({}) # No figures yet


@app.route('/download', methods=['POST'])
def download():
    print(f'Downloading Directory: {webapp.pathDir}')
    memory_file = io.BytesIO()
    with zipfile.ZipFile(memory_file, 'w', zipfile.ZIP_DEFLATED) as zf:
        for root, dirs, files in os.walk(webapp.pathDir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, webapp.pathDir)
                zf.write(file_path, arcname)
    memory_file.seek(0)

    # Define file tag
    if webapp.datasetTagMotif:
        tag = webapp.datasetTagMotif
    else:
        tag = webapp.datasetTag

    return send_file(
        memory_file,
        mimetype='application/zip',
        as_attachment=True,
        download_name=f'{webapp.enzymeName}-{tag}'
    )


@app.route('/jobSummary')
def jobSummary():
    print('Job Summary')
    return render_template('results.html',
                           parameters=webapp.jobParams())


@app.route('/evalFormDNA', methods=['POST'])
def evalDNA():
    # Process the dset
    webapp.evalDNA(parseForm())
    print('Job Done: Eval DNA')
    webapp.done = True
    return render_template('results.html', parameters=webapp.jobParams)


@app.route('/evalFormFilterAA', methods=['POST'])
def filterSubs():
    # Parse the form
    error = webapp.evalSubs(parseForm())
    print('Job Done: Fix AA')
    webapp.done = True
    return render_template('results.html', parameters=webapp.jobParams)


@app.route('/evalFormFilterMotif', methods=['POST'])
def filterMotif():
    webapp.evalSubs(parseForm(), filterMotifs=True)
    webapp.done = True
    print(f'Motif:\n{list(webapp.motifPos.items())}')
    return render_template(
        'setEntropy.html', parameters=webapp.jobParams,
        minS=webapp.minS, motifPos=list(webapp.motifPos.items()))


@app.route('/updateFig', methods=['POST'])
def updateFig():
    print('Update Figures:')
    json = request.get_json()
    webapp.minS = float(json.get('entropy'))
    webapp.minES = float(json.get('minES'))
    webapp.minESRel = float(json.get('minESRel'))
    webapp.selectMotifPos()
    webapp.plotEntropy()
    data = {
        'status': 'success',
        'entropy': webapp.figures.get('entropy'),
        'motifPos': list(webapp.motifPos.items())  # ← list of pairs
    }
    print(f'Motif: {data["motifPos"]}')
    return jsonify(data)


@app.route('/setEntropy', methods=['GET'])
def setEntropy():
    print('Set Entropy')
    print(f'Motif (SE): {list(webapp.motifPos.items())}')
    return render_template('setEntropy.html', minS=webapp.minS,
                           minES=webapp.minES, minESRel=webapp.minESRel,
                           parameters=webapp.jobParams,
                           motifPos=list(webapp.motifPos.items()))


@app.route('/comet', methods=['POST'])
def comet():
    print('Start Job: COMET')
    json = request.get_json()
    webapp.filterMotifs() ##
    return render_template('results.html', parameters=webapp.jobParams,
                           motifPos=list(webapp.motifPos.items()))


# Run the app
if __name__ == '__main__':
    app.run(threaded=True, debug=True, use_reloader=False, port=9090)
    # sys.exit(0)
