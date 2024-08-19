import os
import sys
import subprocess

if getattr(sys, 'frozen', False):
    # We are running in a PyInstaller bundle
    bundle_dir = sys._MEIPASS
else:
    # We are running in a normal Python environment
    bundle_dir = os.path.dirname(os.path.abspath(__file__))

# Set the PATH to include the bundled executables
os.environ['PATH'] = os.path.join(bundle_dir, 'mamba_env', 'bin') + os.pathsep + os.environ['PATH']

# Run the ASCI executable
asci_path = os.path.join(bundle_dir, 'ASCI')
subprocess.run([asci_path])
