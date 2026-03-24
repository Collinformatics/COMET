from flask import (Flask, jsonify, redirect, render_template, request,
                   send_from_directory, url_for)
from flask_wtf.csrf import CSRFProtect, generate_csrf
from functions import WebApp
import sys
import threading


# Set up the app
app = Flask(__name__)
csrf = CSRFProtect(app)


# Initialize: Application
webapp = WebApp()
webapp.genKey(app)

# Figure storage
figures = {}


def parseForm():
    data = {}

    # Parse form
    keys = set(request.form.keys())
    for key in keys:
        values = request.form.getlist(key)
        if len(values) == 1:  # Check if there are multiple values
            data[key] = values[0]
        else:
            data[key] = values  # Store as a list if multiple values exist

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



def jobID(form, analysis):
    webapp.jobID(form, analysis)
    print('Redirect')
    return redirect(url_for('jobSummary'))



@app.route('/jobSummary')
def jobSummary():
    print('Job Summary')
    return render_template('results.html',
                           parameters=webapp.jobParams())


@app.route('/evalFormDNA', methods=['POST'])
def evalDNA():
    form = parseForm() # Parse the form
    jobID(form, 'evalDNA') # Evaluate job request

    # Process the data
    webapp.evalDNA()
    print('Done')
    return render_template('results.html',
                           parameters=webapp.jobParams)



# @app.route('/evalFormDNA', methods=['POST'])
# def evalDNA():
#     # Parse the form
#     form = parseForm()
#
#     # Kick off background job
#     threading.Thread(target=webapp.evalDNA, args=(form,)).start()
#
#     # Immediately return "job started" page
#     return render_template('results.html', parameters=webapp.jobParams)



@app.route('/data/figures/<filename>')
def getFigure(filename):
    return send_from_directory('data/figures', filename)



@app.route('/checkFigures')
def checkFigures():
    # figures = session.get('figures', {})
    figures = webapp.figures
    return jsonify(figures)



@app.route('/results')
def results():

    return render_template('results.html',
                           figures=webapp.figures,
                           parameters=webapp.jobParams)



@app.route('/')
def home():
    # return render_template('home.html')
    token = generate_csrf()
    return render_template('processDNA.html', csrf_token=token)



@app.route('/processDNA')
def processDNA():
    return render_template('processDNA.html')

@app.route('/filterAminoAcids')
def filterAA():
    return render_template('filterAA.html')

@app.route('/filterMotif')
def filterMotif():
    return render_template('filterMotif.html')

# Run the app
if __name__ == '__main__':
    app.run(debug=True, use_reloader=False, port=9090)
    # sys.exit(0)
