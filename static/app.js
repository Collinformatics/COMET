// Filter options
function createAAContainer(containerId) {
    const container = document.getElementById(containerId);
    if (!container) {
        return;
    }
    container.className = 'label-sub';
    const seqLength = parseInt(document.getElementById('seqLength').value);
    const aminoAcids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"];

    container.innerHTML = '';


    for (let i = 1; i <= seqLength; i++) {
        const wrapper = document.createElement('div');
        wrapper.style.flex = '0 0 60px';

        const label = document.createElement('label');
        label.style.color = '#FFF';
        label.style.marginLeft = '1px';
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
        aaGroup.style.marginLeft = '18px';


        aminoAcids.forEach(aa => {
            const aaLabel = document.createElement('label');
            aaLabel.style.color = 'white';
            aaLabel.style.fontSize = '14px';
            aaLabel.style.display = 'flex';
            aaLabel.style.flexDirection = 'column';
            aaLabel.style.alignItems = 'center';
            aaLabel.style.textAlign = 'center';

            const aaCheckbox = document.createElement('input');
            aaCheckbox.type = 'checkbox';
            aaCheckbox.name = `${containerId === 'exclAAContainer' ? 'exclR' : 'fixR'}${i}`;
            aaCheckbox.value = aa;
//            aaCheckbox.style.display = 'flex';

            aaLabel.appendChild(document.createTextNode(aa));
//            aaLabel.appendChild(document.createElement('br')); // Optional: add line break
            aaLabel.appendChild(aaCheckbox);
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


// Combine substrate profiles
function createProfileContainer(containerId = 'profileContainer') {
    const container = document.getElementById(containerId);
    if (!container) return;

    const nProfiles = parseInt(document.getElementById('nProfiles').value);
    const seqLength = parseInt(document.getElementById('seqLength').value);
    const motifLength = parseInt(document.getElementById('motifLength').value);
    container.innerHTML = '';

    // Header label
    container.innerHTML = `
            <div class="help">
                <label>Experimental Data:</label>
                <div class="help-icon">?
                    <span class="help-tooltip">
                        Substrate profiles obtained from "Filter Motif".<br>
                        Acceptable file extension: .pkl<br><br>
                        The "Motif Index" defines the index of the first AA in the motif within the full substrate sequence
                    </span>
                </div>
            </div>
    `;

    for (let i = 1; i <= nProfiles; i++) {
        if (i > seqLength-motifLength) {
            break; // Enforce data boundaries
        }

        // Upload file
        const wrapper1 = document.createElement('div');
        wrapper1.className = 'form-wrapper';
        wrapper1.innerHTML = `
            <label class="label-sub" for="fileExp${i}">Profile ${i}:</label>
            <input type="file" id="fileExp${i}" name="fileExp${i}" accept=".pkl">
        `;
        container.appendChild(wrapper1);

        // Motif start index
        const wrapper2 = document.createElement('div');
        wrapper2.className = 'form-wrapper';
        wrapper2.style.marginBottom = '1px';
        wrapper2.innerHTML = `
            <label class="label-sub" for="fileExp${i}">Motif Index:</label>
            <input type="number" id="idxStart${i}" value="${i}"
                   min="1" max="${seqLength-motifLength}" style="width: 60px;"
                   required>
        `;
        container.appendChild(wrapper2);
    }
}

document.addEventListener('DOMContentLoaded', () => createProfileContainer());
function updateNumProfiles() {
    createProfileContainer();
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
    //console.log('Input Form:', json);
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

function pageCombineProfiles() {
    window.location.href = "/combineProfiles";
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
    console.log('Data:\n', formData);

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

    // Process form
    const form = document.getElementById("filterSubs");
    const csrfToken = form.querySelector('input[name="csrf_token"]').value;
    const formData = new FormData(form);
    formData.delete('csrf_token');
    console.log('Data:\n', formData);

    // Evaluate the form
    if (filter != 'comet') {
        jobID = await processForm(formData);
        console.log('jobID', jobID);
        formData.append('jobID', jobID);
    }

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
                console.log('Redirect:');
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
    } else if (filter == 'motif') {
        console.log('Run: Eval Motif');

        fetch('/evalFormFilterMotif', {
            body: formData,  // Send the actual FormData object, not a JSON
            headers: { 'X-CSRFToken': csrfToken },
            method: 'POST'
        })
        .then(response => {
            if (response.ok) {
                // Redirect
                console.log('Redirect:');
                window.location.href = '/setEntropy';
            } else {
                console.log("ERROR: Filtering motif.");
                alert("ERROR: Filtering motif.");
                }
        })
        .catch(error => {
            console.error('Error:', error);
            alert("An error occurred.");
        });
    } else if (filter == 'comet') {
        console.log('Run: COMET');
        fetch('/comet', {
            body: formData,  // Send the actual FormData object, not a JSON
            headers: { 'X-CSRFToken': csrfToken },
            method: 'POST'
        })
        .then(response => {
            if (response.ok) {
                // Redirect
                console.log('Redirect:');
                window.location.href = '/results';
            } else {
                console.log("ERROR: Running COMET.");
                alert("ERROR: Running COMET.");
                }
        })
        .catch(error => {
            console.error('Error:', error);
            alert("An error occurred.");
        });
    } else {
        alert('No valid filter type was given to: buttonFilterSubs(filter)')
    }
}

// Get figures
function getFigures() {
    const csrfToken = document.querySelector('meta[name="csrf-token"]').content;
    const container = document.getElementById("figures-container");
    if (!container) return;
    const displayed = new Map();  // track what's already shown

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
            console.log('checkFigures comet response:', data);
            if (data.entropy || data.subProfile || data.eMap || data.eMapSc ||
                data.eLogo || data.eLogoMin || data.wLogo || data.words ||
                data.barCounts || data.barRF || data.exp_counts || data.bg_counts) {
                container.innerHTML = ''; // Clear loading message

                if (data.entropy) {
                    addFigure(container, "Entropy", data.entropy);
                }
                if (data.subProfile) {
                    addFigure(container, "Substrate Profile", data.subProfile);
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
                    if (data.eLogoMin) {
                        addFigure(container, "Enrichment Logo", data.eLogoMin);
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
        // Only stop polling when job is done
        fetch('/jobStatus')
            .then(r => r.json())
            .then(status => {
                console.log('Status:', status.jobStatus);
                if (status.jobStatus) {
                    clearInterval(interval);
                    // Append download button to the box
                    const containerBtn = document.getElementById("button-container");
                    if (containerBtn) {
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
                    };
                }
            });
    }, 1000); // poll: 1000 = 1 second
}

function getFiguresCOMET() {
    const csrfToken = document.querySelector('meta[name="csrf-token"]').content;
    const container = document.getElementById("figures-container");
    if (!container) return;

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
            console.log('checkFigures comet response:', data);
            if (data.entropy || data.eMap || data.eMapSc || data.eLogo || data.eLogoMin ||
            data.wLogo || data.words || data.barCounts || data.barRF || data.exp_counts ||
            data.bg_counts) {
                container.innerHTML = ''; // Clear loading message

                if (data.entropy) {
                    addFigure(container, "Entropy", data.entropy);
                }
                if (data.subProfile) {
                    addFigure(container, "Substrate Profile", data.subProfile);
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
                    if (data.eLogoMin) {
                        addFigure(container, "Enrichment Logo", data.eLogoMin);
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

        // Only stop polling when job is done
        fetch('/jobStatus')
            .then(r => r.json())
            .then(status => {
                if (status.done) {
                    clearInterval(interval);
                    // Append download button to the box
                    const containerBtn = document.getElementById("button-container");
                    if (containerBtn) {
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
                    };
                }
                });
    }, 1000); // poll: 1000 = 1 second
}

function addFigure(container, label, fig, fig2 = null) {
    const p = document.createElement('p');
    p.className = 'p2';
    p.textContent = label;
    container.appendChild(p);

    // Add figure
    const img1 = document.createElement('img');
    img1.src = fig;
    //console.log('fig:', img1.src);
    img1.style.maxWidth = '80vw';
    img1.style.height = 'auto';
    img1.style.marginBottom = '20px';
    container.appendChild(img1);

    // Add a second figure
    if (fig2) {
        const img2 = document.createElement('img');
        img2.src = fig2;
        img2.style.maxWidth = '80vw';
        img2.style.height = 'auto';
        img2.style.marginBottom = '20px';
        container.appendChild(img2);
    }
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

function updateMinS() {
    const csrfToken = document.querySelector('meta[name="csrf-token"]').content;
    const minS = document.getElementById('minS').value;
    const minES = document.getElementById('minES').value;
    const minESRel = document.getElementById('minESRel').value;
    console.log('Sending entropy value:', minS); // Debugging statement
    const data = JSON.stringify({
        minS: minS,
        minES: minES,
        minESRel: minESRel
    })

    fetch('/updateFig', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-CSRFToken': csrfToken
        },
        body: data
    })
    .then(response => {
        if (!response.ok) {
            throw new Error('Network response was not ok ' + response.statusText);
        }
        return response.json();
    })
    .then(data => {
        if (data.status === 'success') {
            console.log('Success:', data);

            // update the entropy image
            const imgs = document.querySelectorAll('img');
            imgs.forEach(img => {
                if (img.src.includes('entropy')) {
                    img.src = img.src.split('?')[0] + '?t=' + Date.now();
                }
            });

            // Update motifPos list
            const motifContainer = document.getElementById('motifPos');
            if (motifContainer) {
                const label = motifContainer.querySelector('label');
                motifContainer.innerHTML = '';
                motifContainer.appendChild(label);

                if (data.motifPos && data.motifPos.length > 0) {
                    data.motifPos.forEach(([pos, val]) => {
                        const p = document.createElement('p');
                        p.className = 'p3';
                        p.innerHTML = `${pos}: ∆S=<span class="param-value">${val.toFixed(2)}</span>`;
                        motifContainer.appendChild(p);
                    });
                } else {
                    const p = document.createElement('p');
                    p.className = 'p3';
                    p.textContent = 'No positions selected.';
                    motifContainer.appendChild(p);
                }
            }
        }
    })
    .catch((error) => {
        console.error('Error:', error);
    });
}
