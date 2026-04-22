# Purpose:

COMET (COmprehensive Motif Evaluation Toolkit) was developed to process high-throughput datasets for the evaluation of enzymatic specificity.

<img width="9900" height="3300" alt="image" src="https://github.com/user-attachments/assets/1617cf0c-8ecb-4af4-8b2f-5ba292d1f3fa" />

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

Processing DNA:

- To test the program a trial dataset is available in the folder:

      data/validation/



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

