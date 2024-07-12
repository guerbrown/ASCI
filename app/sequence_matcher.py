import os
import re

def find_matching_sequences(files):
    print(f"Finding matches among files: {files}")
    pattern = r'.*?(\d+)(F|R).*\.ab1'
    matches = {}
    for file in files:
        # Use os.path.basename to get just the filename, not the full path
        filename = os.path.basename(file)
        match = re.match(pattern, filename)
        if match:
            number, direction = match.groups()
            if number not in matches:
                matches[number] = {}
            matches[number][direction] = file
    print(f"Matches found: {matches}")
    return {number: pair for number, pair in matches.items() if len(pair) == 2}

# Example usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        files = sys.argv[1:]
        matching_pairs = find_matching_sequences(files)
        print("Matching pairs:", matching_pairs)
    else:
        print("Usage: python sequence_matcher.py file1.ab1 file2.ab1 ...")
