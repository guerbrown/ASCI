# -*- mode: python ; coding: utf-8 -*-
import os
import sys
from PyQt5 import QtCore
from PyInstaller.utils.hooks import collect_data_files
from PyInstaller.building.build_main import Tree

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
             ],
             hiddenimports=['PyQt5.QtCore', 'PyQt5.QtGui', 'PyQt5.QtWidgets'],
             hookspath=[],
             hooksconfig={},
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)

# Add Qt5 plugins
qtbase_path = os.path.dirname(QtCore.__file__)
qt_plugins = [
    ('platforms', os.path.join(qtbase_path, 'plugins', 'platforms')),
    ('xcbglintegrations', os.path.join(qtbase_path, 'plugins', 'xcbglintegrations')),
    ('platformthemes', os.path.join(qtbase_path, 'plugins', 'platformthemes')),
]
for plugin_name, plugin_path in qt_plugins:
    if os.path.exists(plugin_path):
        a.datas += Tree(plugin_path, prefix=os.path.join('PyQt5', 'Qt', 'plugins', plugin_name))

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

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
          console=True,
          disable_windowed_traceback=False,
          target_arch=None,
          codesign_identity=None,
          entitlements_file=None )
