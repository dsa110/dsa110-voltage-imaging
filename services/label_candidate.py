import argparse
from dsaT3 import labels

def local_run(args):

    if args.candidate is not None:
        if args.label is not None:
            try:
                labels.set_label(args.candidate,args.label)
            except:
                print('Could not set label: ',args.candidate,args.label)

        if args.label=='archive':
            if not args.search:
                print('Searching for voltage files, because you are archiving')
                labels.check_voltages(args.candidate)

    if args.search and args.candidate is not None:
        labels.check_voltages(args.candidate)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--candidate', type=str, default=None, help='Candidate name', required=True)
    parser.add_argument('-l', '--label', type=str, default=None, help='Label ('+str(labels._allowed)+')')
    parser.add_argument('-s', '--search', action='store_true', help='Search for voltage files')
    the_args = parser.parse_args()
    local_run(the_args)
