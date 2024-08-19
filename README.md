```
                   __    ___   ___   __
                  /__\  / __) /  _/ (  )
                 /    \ \__ \ | |_  (  )
                 \_/\_/ (___/ \___\ (__)
                 
         Automated Sanger Consensus Interface
```
# ASCI: Automated Sanger Consensus Interface

ASCI is designed to automate the process of creating consensus sequences from Sanger sequencing data. It provides a user-friendly graphical interface for selecting .ab1 files, pairing forward and reverse sequences, and generating consensus sequences.

Alignments are created using the Muslce alignment algorithm (Edgar, 2004). Please note that this was not designed with custom alignment algorithms in mind but could be adapted if desired. There are no plans to add parameter customization in the future but if there is interest, please open an issue.

## Features

- Intuitive GUI for selecting and pairing .ab1 files
- Automatic forward and reverse sequence alignment
- Consensus sequence generation
- Easy-to-use interface for viewing and downloading results

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/ASCI.git
   cd ASCI
   ./lfs.sh
   ```

2. Set up the Mamba environment and install dependencies:
   ```
   # Create and activate the Mamba environment
   mamba env create -f documentation/environment.yml
   mamba activate asci

   # Install additional dependencies
   pip install -r documentation/requirements.txt
   ```

3. Install system dependencies (for Ubuntu/Debian-based systems):
   ```
   sudo apt update
   sudo apt install qt5-default qtcreator build-essential qtbase5-dev qt5-qmake qttools5-dev-tools python3-dev libx11-xcb1 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 libxcb-xinerama0 libxkbcommon-x11-0
   ```

## Usage

1. Activate the Mamba environment:
   ```
   mamba activate asci
   ```

2. Run the ASCI application:
   ```
   ./ASCI
   ```

3. Use the GUI to select your .ab1 files, create pairs, and generate consensus sequences.

4. Download the resulting consensus sequences as a FASTA file.

## Input Data Structure

- .ab1 files: Sanger sequencing output files
- File naming convention: Files should be named with a common prefix followed by a number and 'F' or 'R' for forward and reverse sequences, respectively (e.g., 'Sample-1F.ab1', 'Sample-1R.ab1')

## Troubleshooting

If you encounter any issues, please check the following:

- Ensure all dependencies are correctly installed
- Verify that your .ab1 files follow the expected naming convention
- Check that you have the necessary permissions to read/write in the application directory

For further assistance, please open an issue on the GitHub repository.

## Contributing

Contributions to ASCI are welcome! Please feel free to submit a Pull Request.

## Acknowledgments

- BioPython for sequence handling
- MUSCLE for sequence alignment
- PyQt5 for the graphical user interface
- Claude 3.5 Sonnet for assistance building the C++ application and compiling it
