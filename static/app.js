// Filter options
function createAAContainer(containerId) {
    const container = document.getElementById(containerId);
    if (!container) {
        return;
    }
    const seqLength = parseInt(document.getElementById('seqLength').value);
    const aminoAcids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"];

    container.innerHTML = '';


    for (let i = 1; i <= seqLength; i++) {
        const wrapper = document.createElement('div');
        wrapper.style.flex = '0 0 60px';

        const label = document.createElement('label');
        label.style.color = '#FFF';
        label.appendChild(document.createTextNode(`R${i}`));

        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.name = `filterPos`;
        checkbox.value = `R${i}`;

        label.prepend(checkbox);
        wrapper.appendChild(label);

        const aaGroup = document.createElement('div');
        aaGroup.style.display = 'none';
        aaGroup.style.flexWrap = 'wrap';
        aaGroup.style.marginLeft = '10px';

        aminoAcids.forEach(aa => {
            const aaLabel = document.createElement('label');
            aaLabel.style.color = 'white';
            aaLabel.style.fontSize = '12px';

            const aaCheckbox = document.createElement('input');
            aaCheckbox.type = 'checkbox';
            aaCheckbox.name = `${containerId === 'exclAAContainer' ? 'exclR' : 'fixR'}${i}`;
            aaCheckbox.value = aa;

            aaLabel.appendChild(aaCheckbox);
            aaLabel.appendChild(document.createTextNode(aa));
            aaGroup.appendChild(aaLabel);
        });

        checkbox.addEventListener('change', () => {
            aaGroup.style.display = checkbox.checked ? 'flex' : 'none';
        });

        wrapper.appendChild(aaGroup);
        container.appendChild(wrapper);
    }
}
document.addEventListener('DOMContentLoaded', () => {
    // Monitor AA selection boxes for changes in seqLength
    createAAContainer('fixAAContainer');
    createAAContainer('exclAAContainer');
});
function updateFixedAA() {
    createAAContainer('fixAAContainer');
    createAAContainer('exclAAContainer');
}



async function processForm(formData) {
    const json = {}; // Dont send files as a JSON
    const selectedFixPositions = [];

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
    jobID += date;
    const encoder = new TextEncoder();
    const data = encoder.encode(jobID);
    const hashBuffer = await crypto.subtle.digest('SHA-256', data);
    const hashArray = Array.from(new Uint8Array(hashBuffer));
    jobID = hashArray.map(b => b.toString(16).padStart(2, '0')).join('');

    // Log the form
    console.log('Input Form:', json);
    return jobID;
}


//
function pageHome() {
    window.location.href = "/"
}


function pageProcessDNA() {
    window.location.href = "/processDNA";
}


function pageFilterAA() {
    window.location.href = "/filterAA";
}


function pageFilterMotif() {
    window.location.href = "/filterMotif";
}

function pageResources() {
    window.location.href = "/resources";
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
            console.log('ERROR: Processing DNA.');
            alert("ERROR: Processing DNA.");
        }
    })
    .catch(error => {
        console.error('Error:', error);
        alert("An error occurred.");
    });
}

async function buttonFilterSubs(filter) {
    // Disable to prevent double click
    const button = document.querySelector('.button');
    button.disabled = true;
    button.textContent = 'Processing';

    //
    const form = document.getElementById("filterSubs");
    const csrfToken = form.querySelector('input[name="csrf_token"]').value;
    const formData = new FormData(form);
    formData.delete('csrf_token');

    // Evaluate the form
    jobID = await processForm(formData);
    formData.append('jobID', jobID);


    // POST the raw formData to Flask
    if (filter == 'aa') {
        fetch('/evalFormFilterAA', {
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
                console.log("ERROR: Filtering substrates.");
                alert("ERROR: Filtering substrates.");
            }
        })
        .catch(error => {
            console.error('Error:', error);
            alert("An error occurred.");
        });
    } else {
       fetch('/evalFormFilterMotif', {
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
                console.log("ERROR: Filtering substrates.");
                alert("ERROR: Filtering substrates.");
            }
        })
        .catch(error => {
            console.error('Error:', error);
            alert("An error occurred.");
        });
    }
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


function addFigure(container, label, fig, fig2 = null) {
    const p = document.createElement('p');
    p.className = 'p2';
    p.textContent = label;
    container.appendChild(p);

    // Add figure
    const img1 = document.createElement('img');
    img1.src = fig;
    img1.style.maxWidth = '80vw';
    img1.style.height = 'auto';
    container.appendChild(img1);

    // Add the second figure
    if (fig2) {
        container.appendChild(document.createElement('br'));
        const img2 = document.createElement('img');
        img2.src = fig2;
        img2.style.maxWidth = '80vw';
        img2.style.height = 'auto';
        container.appendChild(img2);
    }
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

                // Append download button to the box
                const box = document.querySelector('.box');
                const buttonWrapper = document.createElement('div');
                buttonWrapper.className = 'button-wrapper';
                const button = document.createElement('button');
                button.className = 'button';
                button.textContent = 'Download';
                button.onclick = download;
                buttonWrapper.appendChild(button);
                document.getElementById('button-container').appendChild(buttonWrapper);
                //container.appendChild(buttonWrapper);

                if (data.entropy) {
                    addFigure(container, "Entropy", data.entropy);
                }
                if (data.eMap) {
                    addFigure(container, "Enrichment Map", data.eMap);
                }
                if (data.eMapSc) {
                    addFigure(container, "Scaled Enrichment Map", data.eMapSc);
                }
                if (data.eLogo && data.eLogoMin) {
                    addFigure(container, "Enrichment Logo", data.eLogo, data.eLogoMin);
                } else {
                    if (data.eLogo) {
                        addFigure(container, "Enrichment Logo", data.eLogo);
                    }
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
    const dirName = 'comet.zip';
    const csrfToken = document.querySelector('meta[name="csrf-token"]').content;
    const response = await fetch('/download', {
        body: JSON.stringify({}),
        headers: { 'X-CSRFToken': csrfToken },
        method: 'POST'
    });

    const blob = await response.blob();

    // Extract filename from Content-Disposition header
    let filename = 'comet.zip';
    const disposition = response.headers.get('Content-Disposition');
    if (disposition) {
        const match = disposition.match(/filename[^;=\n]*=([^;\n]*)/);
        if (match && match[1]) {
            filename = match[1].trim().replace(/['"]/g, '');
        }
    }
    console.log('File Name:', filename);

    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}
