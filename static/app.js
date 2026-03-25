// Fix AA
function updateFixedAA() {
    const seqLength = parseInt(document.getElementById('seqLength').value);
    const container = document.getElementById('fixedAAContainer');
    const aminoAcids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"];

    container.innerHTML = '';
    container.style.display = 'flex';
    container.style.flexWrap = 'wrap'; // Wrap to next line if needed
    container.style.gap = '12px'; // spacing between each checkbox
    container.paddingLeft = '5px';

    for (let i = 1; i <= seqLength; i++) {
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
document.addEventListener('DOMContentLoaded', function() {
    updateFixedAA();

    // Inspect File Upload: Experimental data
    document.getElementById('fileExp').addEventListener('change', function(event) {
    const file = event.target.files[0];
    if (!file) return;

    const allowed = ['.fasta', '.fa', '.fasta.gz', '.fa.gz', '.fastq', '.fq', '.fastq.gz', '.fq.gz'];
    const filename = file.name.toLowerCase();

    const valid = allowed.some(ext => filename.endsWith(ext));
    if (!valid) {
        alert(`ERROR:\nYou can only upload files with these extensions:
        ${allowed.join(' ')}`);
        event.target.value = ''; // Clear file input
    }
    });

    // Inspect File Upload: Background data
    document.getElementById('fileBg').addEventListener('change', function(event) {
    const file = event.target.files[0];
    if (!file) return;

    const allowed = ['.fasta', '.fa', '.fasta.gz', '.fa.gz', '.fastq', '.fq', '.fastq.gz', '.fq.gz'];
    const filename = file.name.toLowerCase();

    const valid = allowed.some(ext => filename.endsWith(ext));
    if (!valid) {
        alert(`ERROR:\nYou can only upload files with these extensions:
        ${allowed.join(' ')}`);
        event.target.value = ''; // Clear file input
    }
    });
});


// Define button function
function buttonProcessDNA() {
    // const csrfToken = document.querySelector('input[name="csrf_token"]').value;
    const form = document.getElementById("formDNA");
    const csrfToken = form.querySelector('input[name="csrf_token"]').value;
    const formData = new FormData(form);
    formData.delete('csrf_token');
    const json = {}; // Dont send files as a JSON
    const selectedFixPositions = [];

    // Redirect
//    window.location.href = '/results';

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

            // Push another value into the list
            json[key].push(value);
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

    // Log the form
    console.log('Form:');
    for (let [key, value] of formData.entries()) {
        console.log(key, value);
    }

    // Log the form
    console.log('Input Form:', json);


    // POST the raw FormData to Flask
    fetch('/evalFormDNA', {
        method: 'POST',
        headers: { 'X-CSRFToken': csrfToken },
        body: formData  // Send the actual FormData object, not a JSON
    })
    .then(response => {
        if (response.ok) {
            // Go to the results page after successful processing
            window.location.assign('/results')
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
        .then(data => {
            if (data.exp_counts || data.bg_counts) {
                clearInterval(interval); // Stop polling
                container.innerHTML = ''; // Clear loading message


                if (data.exp_counts) {
                    // Figure label
                    const label = document.createElement('p');

                    label.className = 'p2';
                    label.textContent = "Experimental Counts";
                    container.appendChild(label);

                    // Add figure
                    const img = document.createElement('img');
                    img.src = `/data/figures/${data.exp_counts}`;
                    img.style.maxWidth = '80vw';
                    img.style.height = 'auto';
                    container.appendChild(img);
                    container.appendChild(document.createElement('br'));
                    container.appendChild(document.createElement('br'));
                }

                if (data.bg_counts) {
                    // Figure label
                    const label = document.createElement('p');
                    label.className = 'p2';
                    label.textContent = "Background Counts";
                    container.appendChild(label);

                    // Add figure
                    const img = document.createElement('img');
                    img.src = `/data/figures/${data.bg_counts}`;
                    img.style.maxWidth = '80vw';
                    img.style.height = 'auto';
                    container.appendChild(img);
                    container.appendChild(document.createElement('br'));
                    container.appendChild(document.createElement('br'));
                }
           }

            if (data.eMap) {
                clearInterval(interval); // Stop polling
                container.innerHTML = ''; // Clear loading message

                // Figure label
                const label = document.createElement('p');

                label.className = 'p2';
                label.textContent = "Enrichment Map";
                container.appendChild(label);

                // Add figure
                const img = document.createElement('img');
                img.src = `/data/figures/${data.eMap}`;
                img.style.maxWidth = '80vw';
                img.style.height = 'auto';
                container.appendChild(img);
                container.appendChild(document.createElement('br'));
                container.appendChild(document.createElement('br'));
            }

            if (data.eMapSc) {
                clearInterval(interval); // Stop polling
                container.innerHTML = ''; // Clear loading message

                // Figure label
                const label = document.createElement('p');

                label.className = 'p2';
                label.textContent = "Enrichment Map";
                container.appendChild(label);

                // Add figure
                const img = document.createElement('img');
                img.src = `/data/figures/${data.eMapSc}`;
                img.style.maxWidth = '80vw';
                img.style.height = 'auto';
                container.appendChild(img);
                container.appendChild(document.createElement('br'));
                container.appendChild(document.createElement('br'));
            }

            if (data.eLogo) {
                clearInterval(interval); // Stop polling
                container.innerHTML = ''; // Clear loading message

                // Figure label
                const label = document.createElement('p');

                label.className = 'p2';
                label.textContent = "Enrichment Logo";
                container.appendChild(label);

                // Add figure
                const img = document.createElement('img');
                img.src = `/data/figures/${data.exp_counts}`;
                img.style.maxWidth = '80vw';
                img.style.height = 'auto';
                container.appendChild(img);
                container.appendChild(document.createElement('br'));
                container.appendChild(document.createElement('br'));
            }

            
            /* self.figures = {
                'eMap': False, 'eMapSc': False, 'eLogo': False, 'wLogo': False,
                'words': False, 'barCounts': False, 'barRF': False,
                'exp_counts': False, 'bg_counts': False
            } */
        });
    }, 1000); // poll every 1 second
}
