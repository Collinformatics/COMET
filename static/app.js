// Fix AA
function updateFixedAA() {
    const container = document.getElementById('fixedAAContainer');
    const aminoAcids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"];
    const seqLength = parseInt(document.getElementById('seqLength').value);

    container.innerHTML = '';
    container.style.display = 'flex';
    container.style.flexWrap = 'wrap'; // Wrap to next line if needed
    container.style.gap = '12px'; // spacing between each checkbox
    container.paddingLeft = '5px';

    for (let i = 1; i <= seqLength.value; i++) {
        const wrapper = document.createElement('div');
        wrapper.style.display = 'flex';
        wrapper.style.flexDirection = 'column';
        wrapper.style.flex = '0 0 60px';

        const label = document.createElement('label');
        label.style.display = 'flex';
        label.style.alignItems = 'center';
        label.style.color = '#FFF'; // '#FA8128'
        label.style.flex = '0 0 30px';

        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.name = `filterPos`;
        checkbox.value = `R${i}`;

        const text = document.createTextNode(`R${i}`);
        label.appendChild(checkbox);
        label.appendChild(text);
        wrapper.appendChild(label); // Append label to wrapper

        const aaGroup = document.createElement('div');
        aaGroup.style.display = 'none';
        aaGroup.style.flexWrap = 'wrap';
        aaGroup.style.gap = '0px';
        aaGroup.style.marginLeft = '10px';
        aaGroup.style.marginTop = '0px';

        aminoAcids.forEach(aa => {
            const aaLabel = document.createElement('label');
            aaLabel.style.display = 'flex';
            aaLabel.style.alignItems = 'center';
            aaLabel.style.color = 'white';
            aaLabel.style.fontSize = '12px';

            const aaCheckbox = document.createElement('input');
            aaCheckbox.type = 'checkbox';
            aaCheckbox.name = `fixR${i}`;
            aaCheckbox.value = aa;

            const aaText = document.createTextNode(aa);
            aaLabel.appendChild(aaCheckbox);
            aaLabel.appendChild(aaText);
            aaGroup.appendChild(aaLabel);
        });

        checkbox.addEventListener('change', function () {
            aaGroup.style.display = checkbox.checked ? 'flex' : 'none';
        });

        wrapper.appendChild(aaGroup);      // Append AA group to wrapper
        container.appendChild(wrapper);    // Append full wrapper to container
    }
}

// Create listener to inspect html input after the form has been loaded and parsed
document.querySelectorAll('input').forEach(input => {
  input.addEventListener('change', yourFunction);
});


async function processForm(formData) {
    const json = {}; // Dont send files as a JSON

    // Process the input form
    for (const [key, value] of formData.entries()) {
        if (key === 'filterPos') {
            selectedFixPositions.push(value);  // e.g., ['R2']
        }
        if (json[key]) {
            // When you have more that one value or a key, put the values in a list
            if (!Array.isArray(json[key])) {
                json[key] = [json[key]]; // Convert to array
            }
            json[key].push(value); // Push another value into the list
        } else {
            json[key] = value;
        }
    }

    // Clean out fixR* keys not selected
    Object.keys(json).forEach(key => {
        if (key.startsWith('fix') && !selectedFixPositions.includes(key.replace('fix', ''))) {
            delete json[key];
        }
    });

    // Evaluate job ID
    let jobID = '';
    // console.log('Form:');
    for (let [key, value] of formData.entries()) {
        if (value instanceof File) {
            if (value.name) {
                // console.log('* File name:', value.name);
                jobID += key + '_' + value.name + ' ';
            }
        } else if (value) {
            // console.log(key, value);
            jobID += key + '_' + value + ' ';
        }
    }
    const date = new Date().toUTCString();
    jobID += date; // jobID.slice(0, -1);
    const encoder = new TextEncoder();
    const data = encoder.encode(jobID);
    const hashBuffer = await crypto.subtle.digest('SHA-256', data);
    const hashArray = Array.from(new Uint8Array(hashBuffer));
    jobID = hashArray.map(b => b.toString(16).padStart(2, '0')).join('');

    // Log the form
    console.log('Input Form:', json);
    return jobID;
}


// Define button function
async function buttonProcessDNA() {
    // Disable to prevent double click
    const button = document.querySelector('button[onclick="buttonProcessDNA()"]');
    button.disabled = true;
    button.textContent = 'Processing';

    // const csrfToken = document.querySelector('input[name="csrf_token"]').value;
    const form = document.getElementById("formDNA");
    const csrfToken = form.querySelector('input[name="csrf_token"]').value;
    const formData = new FormData(form);
    formData.delete('csrf_token');
    const selectedFixPositions = [];

    // Evaluate the form
    jobID = await processForm(formData);
    formData.append('jobID', jobID);

    // POST the raw formData to Flask
    fetch('/evalFormDNA', {
        body: formData,  // Send the actual FormData object, not a JSON
        headers: { 'X-CSRFToken': csrfToken },
        method: 'POST'
    })
    .then(response => {
        if (response.ok) {
            // Redirect
            console.log('Redirect:')
            window.location.href = '/results';
        } else {
            console.log('Error processing DNA.');
            alert("Error processing DNA.");
        }
    })
    .catch(error => {
        console.error('Error:', error);
        alert("An error occurred.");
    });
}


//
function pageHome() {
    window.location.href = "/"
}


function pageProcessDNA() {
    window.location.href = "/processDNA";
}


function pageFilterAA() {
    window.location.href = "/filterAminoAcids";
}


function pageFilterMotif() {
    window.location.href = "/filterMotif";
}


function clickButton() {
    const message = 'your data';  // Pass data to app.py

    // POST request to app.py
    fetch('/run', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json' // Send JSON (Optional)
        },
        body: JSON.stringify({ message: message }) // Send data to app.py
    })
    .then(response => response.text())
    .then(data => {
        console.log('Server Response:\n', data);
    })
    .catch(error => {
        console.error('ERROR: ', error);
    });
}


function addFigure(container, label, filename) {
    const p = document.createElement('p');
    p.className = 'p2';
    p.textContent = label;
    container.appendChild(p);

    const img = document.createElement('img');
    img.src = filename // dir; /* ## */
    img.style.maxWidth = '80vw';
    img.style.height = 'auto';
    container.appendChild(img);
    container.appendChild(document.createElement('br'));
    container.appendChild(document.createElement('br'));
}

// Get figures
function getFigures() {
    const container = document.getElementById("figures-container");
    const csrfToken = document.querySelector('meta[name="csrf-token"]').content;

    const interval = setInterval(() => {
        // new Flask route returning JSON with filenames
        fetch('/checkFigures', {
            method: 'GET',
            headers: { 'X-CSRFToken': csrfToken },
            credentials: 'same-origin'
        })
        .then(res => res.json())
        // Verify if one figure is ready
        .then(data => {
            if (data.eMap || data.eMapSc || data.exp_counts || data.bg_counts || data.eLogo || data.wLogo || data.words) {
                clearInterval(interval); // Stop polling
                container.innerHTML = ''; // Clear loading message

                /* Download  */
                const buttonWrapper = document.createElement('div');
                buttonWrapper.className = 'button-wrapper';
                const button = document.createElement('button');
                button.className = 'btn';
                button.textContent = 'Download';
                button.onclick = download;
                buttonWrapper.appendChild(button);
                container.appendChild(buttonWrapper);

                if (data.eMap) {
                    addFigure(container, "Enrichment Map", data.eMap);
                }
                if (data.eMapSc) {
                    addFigure(container, "Scaled Enrichment Map", data.eMapSc);
                }
                if (data.eLogo) {
                    addFigure(container, "Enrichment Logo", data.eLogo);
                }
                if (data.wLogo) {
                    addFigure(container, "WebLogo", data.wLogo);
                }
                if (data.words) {
                    addFigure(container, "Word Cloud", data.words);
                }
                if (data.barCounts) {
                    addFigure(container, "Substrate Counts", data.barCounts);
                }
                if (data.barRF) {
                    addFigure(container, "Substrate Frequency", data.barRF);
                }
                if (data.exp_counts) {
                    addFigure(container, "Experimental Counts", data.exp_counts);
                }
                if (data.bg_counts) {
                    addFigure(container, "Background Counts", data.bg_counts);
                }
            }
        });
    }, 1000); // poll every 1 second
}

async function download() {
    const enzymeInput = document.getElementById("enzymeName");
    const dirName = enzymeInput ? enzymeInput.value + '.zip' : 'comet.zip';
    console.log('Dir Name:', dirName)
    const csrfToken = document.querySelector('meta[name="csrf-token"]').content;
    const response = await fetch('/download', {
        body: JSON.stringify({}),
        headers: { 'X-CSRFToken': csrfToken },
        method: 'POST'
    });

    const blob = await response.blob();


    // Extract filename from Content-Disposition header
    const disposition = response.headers.get('Content-Disposition');
    let filename = dirName;
    if (disposition) {
        console.log('Dis: True');
        const match = disposition.match(/filename[^;=\n]*=([^;\n]*)/);
        if (match && match[1]) {
            filename = match[1].trim().replace(/['"]/g, '');
        }
    } else {
        console.log('Dis: False');
    }
    console.log('File Name:', filename);

    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}
