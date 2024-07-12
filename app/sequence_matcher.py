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

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        matching_pairs = find_matching_sequences(sys.argv[1])
        print(matching_pairs)
    else:
        print("Usage: python sequence_matcher.py /path/to/ab1/files")
