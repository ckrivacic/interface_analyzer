"""
Usage: plot_relaxed.py [options]

Options:
    --filename=STR  Only plot a single origin
    -x STR  What to plot on x-axis  [default: origin]
    -y STR  What to plot on y-axis  [default: dG_separated]
"""
from matplotlib import pyplot as plt
import glob, os
import pandas as pd
import docopt


def read_and_calculate(pdb_path):
    """
    Parse PDB file for scores
    """
    
    record = {'path': os.path.basename(pdb_path)}
    origin = os.path.basename(pdb_path).split('.')[0].split('_')[:-1]
    origin = '_'.join(origin)
    print(origin)
    record['origin'] = origin

    with open(pdb_path, 'r') as f:
        lines = f.readlines()

    metrics = [
            'dG_cross',
            'dG_cross/dSASAx100',
            'dG_separated',
            'dG_separated/dSASAx100',
            'dSASA_hphobic',
            'dSASA_int',
            'dSASA_polar',
            'delta_unsatHbonds',
            'hbond_E_fraction',
            'hbond_E_fraction',
            'nres_all',
            'nres_int',
            'packstat',
            'per_residue_energy_int',
            'sc_value',
            'side1_normalized',
            'side1_score',
            'side2_normalized',
            'side2_score',
            'complex_normalized'
            ]

    for line in lines:
        #line = line.decode('utf8')

        if line.startswith('pose'):
            record['total_score'] = float(line.split()[-1])

        elif len(line.split()) == 2:
            if line.split()[0] in metrics:
                record[line.split()[0]] = float(line.split()[1])

    return record


def load_pdbs(folder, filename=None):
    pdbfiles = sorted(glob.glob(folder + '/*.pdb'))
    pdbfiles_final = []
    if filename:
        for pdb in pdbfiles:
            if os.path.basename(pdb).startswith(filename):
                pdbfiles_final.append(pdb)
    else:
        pdbfiles_final = pdbfiles

    records = []
    for pdb in pdbfiles_final:
        print('Reading pdb file {}'.format(pdb))
        records.append(read_and_calculate(pdb))

    records = pd.DataFrame(records)
    return records


def plot(records, x='origin', y='dG_separated'):
    xvals = records[x]
    yvals = records[y]

    fig, ax = plt.subplots()

    if x=='origin':
        yvals = []
        labels = list(set(xvals))
        labels_str = [str(l) for l in labels]
        for label in labels:
            record = records[records[x] == label]
            yvals.append(record[y].values)

        #width = 0.45
        #widthmodifier = 0.3
        #widthlist = widthmodifier * len(labels) * width/(len(labels))

        vplt = ax.violinplot(yvals,
                #positions=[ind for ind in range(len(labels))],
                showmedians=True)

        for j, v in enumerate(yvals):
            print(j, v)
            pos = max(v)# + 0.05 * abs(max(v))
            print(pos)
            ax.text(j, pos, labels[j])

        ax.set_xlabel('Parent structure')
        ax.set_ylabel(y)
        ax.set_xticklabels(labels_str)
    
    plt.show()

if __name__=='__main__':
    args = docopt.docopt(__doc__)
    x = args['-x']
    y = args['-y']
    filename = args['--filename']
    folder = 'relaxed_outputs'
    records = load_pdbs(folder)
    plot(records, x=x, y=y)
