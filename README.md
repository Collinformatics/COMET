# Purpose:

No

COMET (COmprehensive Motif Evaluation Toolkit) was developed to process high-throughput datasets for the evaluation of enzymatic specificity.

# Installation:

Create Environment:

    conda env create -f environment.yml

Activate Environment

    conda activate comet

# Host Website

You can use your computer to host the website.

    python app.py

- This will generate a URL with ip address and port, likely http://127.0.0.1:9090
  - Click on it to go to the site, or enter the address in a browser.

# Terminating The Server:

Leaving the process running can be a significant drain on your battery, so its best to turn things off after you are done.

The commands you enter into the terminal will depend upon your OS:

- MacOS:

    List open files on port 5000:

        lsof -i :5000

    - Inspect output, and look for "Python" in the "COMMAND" column:

            COMMAND     PID     NODE NAME
            Python    49231 ... 0t0  TCP localhost:commplex-main (LISTEN)
            Python    49288 ... 0t0  TCP localhost:commplex-main (LISTEN)
            Python    49288 ... 0t0  TCP localhost:commplex-main (LISTEN)

    Use the Process ID (49231, 49288) to terminate these processes:

         kill 49231 49288 

- Windows

    Inspect network connections on port 5000

        netstat -ano | findstr :5000
    
     - Inspect output, and look for the process ID (20352):

            TCP    127.0.0.1:5000         0.0.0.0:0              LISTENING       20352

    Lists all running tasks with matching process ID:

        tasklist /FI "PID eq 20352"

    - Inspect output, and look for "python.exe" in the "Image Name" column:
  
            Image Name                     PID Session Name        Session#    Mem Usage
            ========================= ======== ================ =========== ============
            python.exe                   20352 Console                    3     39,732 K

    Kill the process:

        taskkill /PID 20352 /F
