import argparse
import yaml

def load_raw_params(f):
    return yaml.load(f.read())

def load_names(f):
    names = {}
    for (lnum, line) in enumerate(f):
        if line.startswith('#!!'):
            asides = [x.strip() for x in line[3:].split('=')]
            assert len(asides) == 2, """assignment command on line {} is invald
                                     """.format(lnum+1)
            names[asides[0]] = asides[1]
            print("{} = {}".format(asides[0], asides[1]))
    return names

def parse_s(s, names):
    new_s = ""
    p = 0
    while p < len(s):
        start = s.find('&', p)
        if start == -1:
            new_s = new_s + s[p:]
            break
        end = s.find('&', start+1)
        assert end != -1, """there was an opening '&' without a
                             closing '&' in {}[{}]""".format(s, p)
        key = elem[start+1:end]
        assert key in names.keys(), key + " not found in names"
        new_s = new_s + s[p:start] + names[key]
        p = end + 1
    return new_s

def parse_raw_params(params, names):
    if len(names) == 0:
        return
    else:
        for (key, cat) in params.items():
            if key == 'version':
                continue
            elif key == 'comments':
                for (i, comment) in enumerate(params['comments']):
                    params['comments'][i] = parse_s(comment, names)
            for (key, elem) in cat.items():
                cat[key] = parse_s(elem, names)
    return

parser = argparse.ArgumentParser(description="Compile input files for a parametric study")
parser.add_argument('fname', help='file name of pstudy input file')
parser.add_argument('--verbose', '-v', action='store_true', default=False,
                    help='print diagnostic information')
args = parser.parse_args()

with open(args.fname) as f:
    params = load_raw_params(f)
    f.seek(0)
    names = load_names(f)
    print("names: ", names)
    parse_raw_params(params, names)
    for (key, vals) in params.items():
        print(key + ": ")
        for val in vals:
            print(val)
        print()

