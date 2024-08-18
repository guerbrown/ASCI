import subprocess
import sys
import os

if __name__ == "__main__":
    # Get the directory of the executable
    if getattr(sys, 'frozen', False):
        application_path = sys._MEIPASS
    else:
        application_path = os.path.dirname(os.path.abspath(__file__))

    # Path to the ASCI executable
    asci_path = os.path.join(application_path, 'ASCI')

    # Run the ASCI executable
    subprocess.run([asci_path])
