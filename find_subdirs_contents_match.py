#!/usr/bin/env python3

import sys, os, fnmatch
if len(sys.argv) < 3:
    print("Too few arguments!\nUsage: " + sys.argv[0].split("/")[-1] + " <pattern> <starting_path>")
    sys.exit(1)

def find(pattern, path):
    size = 0.0
    result = []
    subdirs = set()
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                filepath = os.path.join(root, name)
                size += os.path.getsize(filepath) / 1048576.
                result.append(filepath)
                subdirpath = '/'.join(filepath.split('/')[:-1])
                subdirs.add(subdirpath)
    return result, subdirs, size

pat = sys.argv[1]
rootpath = sys.argv[2]
files, subdirs, size = find(pat, rootpath)
print('\n'.join(['"' + subdir + '"' for subdir in subdirs]))
#print('\n'.join(subdirs))

#print('\nTotal Size: {0:.2f} MiB'.format(size))
