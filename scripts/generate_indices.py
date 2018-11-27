from ont_fast5_api.fast5_file import Fast5File
from argparse import ArgumentParser
import sys
import glob


def main():
    args = get_args()
    for filename in glob.iglob(args.dir + '**/**/*.fast5', recursive=True):
        extract_read_id(filename)


def extract_read_id(fn):
    f = Fast5File(fn, "r+")
    read_id = f.status.read_info[0].read_id
    f.close()
    print(read_id, fn, sep='\t')

def get_args():
    parser = ArgumentParser(description="Extract signal level from fast5 files")
    parser.add_argument("-d", "--dir", help="directory with fast5 file(s) to extract signal from")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
