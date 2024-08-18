# -*- mode: python ; coding: utf-8 -*-
import sys
import os
from PyQt5 import QtCore

block_cipher = None

# Path to your Mamba environment
mamba_env = '/home/guerbrown/software/miniforge3/envs/asci'

a = Analysis(['ASCI_launcher.py'],
             pathex=['.', mamba_env],
             binaries=[],
             datas=[
                 ('sequence_aligner.py', '.'), 
                 ('main.py', '.'),
                 ('ab1_to_fasta.py', '.'),
                 ('sequence_matcher.py', '.'),
                 ('ASCI', '.'),
                 # Include the entire Mamba environment
                 (mamba_env, 'mamba_env')
             ],
             hiddenimports=['Bio', 'Bio.Seq', 'Bio.SeqIO', 'Bio.Align', 'Bio.Align.AlignInfo'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)

# Include Qt plugins
qt_plugins = [
    ('platforms', QtCore.QLibraryInfo.location(QtCore.QLibraryInfo.PluginsPath) + '/platforms'),
    ('platforms/libqxcb.so', QtCore.QLibraryInfo.location(QtCore.QLibraryInfo.PluginsPath) + '/platforms/libqxcb.so')
]
a.datas += qt_plugins

pyz = PYZ(a.pure, a.zipped_data,
          cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='ASCI',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )  # Changed to True for debugging

# Modify the ASCI_launcher.py to use the bundled Mamba environment
with open('ASCI_launcher.py', 'a') as f:
    f.write('''
import os
import sys

# Add the bundled Mamba environment to sys.path
mamba_env_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mamba_env')
sys.path.insert(0, mamba_env_path)
os.environ['PATH'] = mamba_env_path + os.pathsep + os.environ['PATH']
''')
