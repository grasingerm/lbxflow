import argparse
import os

parser = argparse.ArgumentParser(description='Rename files to ffmpeg')
parser.add_argument('dir', help='directory where files to be renamed are stored')
parser.add_argument('--exe', '-x', default='.png', 
                    help='file extension of images')
parser.add_argument('--numdigits', '-n', default=9, 
                    help='number of digits to format number to')
parser.add_argument('--verbose', '-v', default=False, action='store_true', 
                    help='print diagnostic information')
args = parser.parse_args()

_format_str = "{0:0" + str(args.numdigits) + "d}"

for (i, fname) in enumerate(os.listdir(args.dir)):
    if args.verbose:
        print('Checking ' + fname + ' to see if it needs renamed')
    if fname.endswith(args.exe):
        pieces = fname.split("_")
        pieces[-1] = _format_str.format(i) + args.exe
        new_fname = "_".join(pieces)
        os.rename(os.path.join(args.dir, fname), os.path.join(args.dir, new_fname))
        if args.verbose:
            print('\t' + fname + ' -> ' + new_fname)
