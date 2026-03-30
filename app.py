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
    print('Data:')
    for key in keys:
        values = request.form.getlist(key)
        if len(values) == 1:  # Check if there are multiple values
            data[key] = values[0]
        else:
            data[key] = values  # Store as a list if multiple values exist
        print('*', key, data[key])
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


@app.route('/jobSummary')
def jobSummary():
    print('Job Summary')
    return render_template('results.html',
                           parameters=webapp.jobParams())


@app.route('/evalFormDNA', methods=['POST'])
def evalDNA():
    # Process the data
    webapp.evalDNA(parseForm())
    print('Job Done: Eval DNA')
    webapp.done = True
    return render_template('results.html',
                           parameters=webapp.jobParams)


@app.route('/fixAA', methods=['POST'])
def fixAA():
    # Parse the form
    webapp.fixAA(parseForm())
    print('Job Done: Fix AA')
    return render_template('results.html',
                           parameters=webapp.jobParams)


@app.route('/')
def home():
    # return render_template('home.html')
    return render_template('processDNA.html',
                           csrf_token=generate_csrf())


@app.route('/processDNA')
def processDNA():
    return render_template('processDNA.html',
                           csrf_token=generate_csrf())


@app.route('/filterAminoAcids')
def filterAA():
    return render_template('filterAA.html',
                           csrf_token=generate_csrf())


@app.route('/filterMotif')
def filterMotif():
    return render_template('filterMotif.html',
                           csrf_token=generate_csrf())


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
    print(f'Downloading file: {webapp.pathData}')
    dir = webapp.pathData
    memory_file = io.BytesIO()

    with zipfile.ZipFile(memory_file, 'w', zipfile.ZIP_DEFLATED) as zf:
        for root, dirs, files in os.walk(dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, dir)
                zf.write(file_path, arcname)

    memory_file.seek(0)
    return send_file(
        memory_file,
        mimetype='application/zip',
        as_attachment=True,
        download_name=f'{webapp.enzymeName}.zip'
    )


# Run the app
if __name__ == '__main__':
    app.run(threaded=True, debug=True, use_reloader=False, port=9090)
    # sys.exit(0)
