import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QProgressBar, QFileDialog
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from main import process_files
import os

os.environ['QT_QPA_PLATFORM'] = 'minimal'

class Worker(QThread):
    progress = pyqtSignal(int)
    finished = pyqtSignal(str)

    def __init__(self, input_dir, output_dir):
        super().__init__()
        self.input_dir = input_dir
        self.output_dir = output_dir

    def run(self):
        # Simulate progress
        for i in range(101):
            self.progress.emit(i)
            self.msleep(50)
        
        consensus_file = process_files(self.input_dir, self.output_dir)
        self.finished.emit(consensus_file)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ASCI")
        self.setGeometry(100, 100, 300, 150)

        self.button = QPushButton("Drop AB1 Files", self)
        self.button.setGeometry(75, 20, 150, 40)
        self.button.clicked.connect(self.select_files)

        self.progress = QProgressBar(self)
        self.progress.setGeometry(50, 70, 200, 25)
        self.progress.hide()

        self.download_button = QPushButton("Download Consensus", self)
        self.download_button.setGeometry(75, 100, 150, 40)
        self.download_button.clicked.connect(self.download_consensus)
        self.download_button.hide()

        self.consensus_file = None

    def select_files(self):
        input_dir = QFileDialog.getExistingDirectory(self, "Select Input Directory")
        if input_dir:
            output_dir = QFileDialog.getExistingDirectory(self, "Select Output Directory")
            if output_dir:
                self.progress.show()
                self.worker = Worker(input_dir, output_dir)
                self.worker.progress.connect(self.update_progress)
                self.worker.finished.connect(self.process_finished)
                self.worker.start()

    def update_progress(self, value):
        self.progress.setValue(value)

    def process_finished(self, consensus_file):
        self.consensus_file = consensus_file
        self.download_button.show()

    def download_consensus(self):
        if self.consensus_file:
            save_path, _ = QFileDialog.getSaveFileName(self, "Save Consensus File", "", "FASTA Files (*.fasta)")
            if save_path:
                import shutil
                shutil.copy(self.consensus_file, save_path)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
