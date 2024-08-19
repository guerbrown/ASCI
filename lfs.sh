#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Function to print messages
print_message() {
    echo "==== $1 ===="
}

# Check if we're in a git repository
if ! git rev-parse --is-inside-work-tree > /dev/null 2>&1; then
    echo "Error: Not in a Git repository. Please run this script from within a Git repository."
    exit 1
fi

# List LFS files
print_message "Listing LFS files"
git lfs ls-files

# Pull LFS objects
print_message "Pulling LFS objects"
git lfs pull

# Check LFS files integrity
print_message "Checking LFS files integrity"
git lfs fsck

# List LFS locks
print_message "Listing LFS locks"
git lfs locks

print_message "Git LFS operations completed"
