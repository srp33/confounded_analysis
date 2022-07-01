import glob
import sys

in_file_pattern = sys.argv[1]
out_file_path = sys.argv[2]

in_file_paths = sorted(glob.glob(in_file_pattern))

with open(out_file_path, "w") as out_file:
    with open(in_file_paths[0]) as in_file:
        for line in in_file:
            out_file.write(line)

    for in_file_path in in_file_paths[1:]:
        with open(in_file_path) as in_file:
            in_file.readline()

            for line in in_file:
                out_file.write(line)
