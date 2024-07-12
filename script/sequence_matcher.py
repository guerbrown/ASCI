import re
import os

def find_matching_sequences(directory):
    files = os.listdir(directory)
    pattern = r'.*?(\d+)(F|R).*\.ab1'
    matches = {}
    
    for file in files:
        match = re.match(pattern, file)
        if match:
            number, direction = match.groups()
            if number not in matches:
                matches[number] = {}
            matches[number][direction] = file
    
    return [pair for pair in matches.values() if len(pair) == 2]

# Usage
# matching_pairs = find_matching_sequences("path/to/ab1/files")