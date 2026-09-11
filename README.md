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

- This will generate a URL with an ip address and port, it should be http://127.0.0.1:9090

  - Click on the link to go to the website, or enter the address in a browser.


# Testing The Program:

To test the program a trial dataset is available in the "TemplateData" folder, this includes:

- Fastq files

- Translated protein substrates, and AA counts

- Substrates filtered for Q@R5

- Figures from Process DNA, and Filter AA 

Additional instructions can be found on the website's home page.

# Troubleshooting:

Unterminated processes can result in the website not starting correctly.

Two possible solutions are:

1) Automated Fix:

    To fix this problem, execute this command to terminate the python processes:
    
        ./killServer.sh

    If you host the website at a port other than 9090, add the port to the command:
    
        ./killServer.sh <port>

2) Manual Fix:

    List Open Files at port 9090 to find relevant process IDs:

        lsof -i :9090

    Kill these processes:

        kill <process ID>

