from ont_fast5_api.fast5_file import Fast5File
from argparse import ArgumentParser
import sys
import glob


def main():
    args = get_args()
    if args.file:
        for filename in args.file:
            extract_signal(filename)
    else:
        for filename in glob.iglob(args.dir + '**/*.fast5', recursive=args.recursive):
            extract_signal(filename)


def extract_signal(fn):
    """
    Extract raw data from the fast5 file.
    Important parameter could be scale
    If True, returns scaled floating point values in pA
    if False, returns raw DAQ values as 16 bit integers
    """
    f = Fast5File(fn, "r")
    print("{}\t{}".format(fn,
                          ','.join([str(i) for i in list(f.get_raw_data(scale=False))])))


def get_args():
    parser = ArgumentParser(description="Extract signal level from fast5 files")
    parser.add_argument("-d", "--dir", help="directory with fast5 file(s) to extract signal from")
    parser.add_argument("-f", "--file", help="fast5 file(s) to extract signal from", nargs='*')
    parser.add_argument("-r", "--recursive",
                        help="recursively go through directories",
                        action="store_true")
    args = parser.parse_args()
    if not args.dir and not args.file:
        sys.exit("ARGUMENT ERROR: Exactly one of --dir or --file is required")
    if args.dir and args.file:
        sys.exit("ARGUMENT ERROR: Exactly one of --dir or --file is required")
    if args.recursive and args.file:
        sys.exit("ARGUMENT ERROR: -r/--recursive is only applicable when selecting a directory")
    return args


if __name__ == '__main__':
    main()
