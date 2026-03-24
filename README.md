# Purpose:

COMET (COmprehensive Motif Evaluation Toolkit) was developed to process high-throughput datasets for the evaluation of enzymatic specificity.

# Installation:

Clone the GitHub with the terminal command:

    git clone https://github.com/Collinformatics/COMET

Then move to the COMET directory:

    cd COMET

Create conda environment:

    conda env create -f environment.yml


# Host Website

Activate conda environment:

    conda activate comet

You can start up the website with:

    python app.py

- This will generate a URL with an ip address and port, likely http://127.0.0.1:9090
  - Click on the link to go to the website, or enter the address in a browser.


# Terminating Processes:

Unterminated processes can result in the website not starting correctly.

To fix this run this script to terminate the python processes at the port used to host the server:

    ./killServer.sh
