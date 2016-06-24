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
    return names

def parse_s(s, names):
    new_s = ''
    p = 0
    while p < len(s):
        start = s.find('&', p)
        if start == -1:
            new_s = new_s + s[p:]
            break
        end = s.find('&', start+1)
        assert end != -1, """there was an opening '&' without a
                             closing '&' in {}[{}]""".format(s, p)
        key = s[start+1:end]
        assert key in names.keys(), key + ' not found in names'
        new_s = new_s + s[p:start] + names[key]
        p = end + 1
    return new_s

def parse_raw_params(params, names, tnl=False):
    if len(names) == 0:
        return
    else:
        for (idx, param) in enumerate(params):
            if type(param) == str:
                params[idx] = parse_s(param, names)
                if tnl: params[idx] = params[idx].rstrip()
            elif type(param) == dict:
                for (k, v) in param.items():
                    param[k] = parse_s(v, names)
                    if tnl: param[k] = param[k].rstrip()
            else:
                raise Exception('unexpected type')
    return

def compile_file(fname, strs, ext='.yaml'):
    with open(fname + ext, 'w') as w:
        w.write(('\n'.join(strs)).replace('~DIR_LS~', 
                        ', '.join(['"'+ s +'"' for s in fname.split('_')])))

def compile_input_files(fname, *cats):
    if type(cats[-1]) != dict:
        compile_file(fname[:-1], cats)
    else:
        for (k, v) in cats[-1].items():
            if len(k) != 0:
                new_fname = k + "_" + fname
            else:
                new_fname = fname
            compile_input_files(new_fname, v, *cats[:-1])

parser = argparse.ArgumentParser(description='Compile input files for a parametric study')
parser.add_argument('fname', help='file name of pstudy input file')
parser.add_argument('--verbose', '-v', action='store_true', default=False,
                    help='Print diagnostic information')
parser.add_argument('--tnl', '-N', action='store_true', default=False,
                    help='Remove trailing new lines')
args = parser.parse_args()

with open(args.fname) as f:
    params = load_raw_params(f)
    f.seek(0)
    names = load_names(f)
    if args.verbose: print('names: ', names)
    parse_raw_params(params, names, tnl=args.tnl)
    ncombos = 1
    if args.verbose:
        for param in params:
            if type(param) == str:
                print(param)
                print()
            elif type(param) == dict:
                ncombos = ncombos * len(param)
                for (key, val) in param.items():
                    print(key + ": " + val)
                    print()
            else:
                raise Exception('unexpected type')
    
    compile_input_files("", *params)
    print("Compilation complete")
    if args.verbose: print("Total files: ", ncombos)
