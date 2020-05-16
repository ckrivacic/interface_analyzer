"""
Usage: plot_relaxed.py [options]

Options:
    --filename=STR  Only plot a single origin
    -x STR  What to plot on x-axis  [default: origin]
    -y STR  What to plot on y-axis  [default: dG_separated]
    --best=STR  Return the best based on the following metric
"""
from matplotlib import pyplot as plt
import sys
import glob, os
import pandas as pd
import docopt


def read_and_calculate(pdb_path):
    """
    Parse PDB file for scores
    """
    
    record = {'filename': os.path.basename(pdb_path)}
    record['path'] = pdb_path
    origin = os.path.basename(pdb_path).split('.')[0].split('_')[:-1]
    origin = '_'.join(origin)
    #print(origin)
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


def load_pdbs(folder, filename=None, return_best=False):
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
    if return_best:
        from shutil import copyfile
        print('BEST SCORING MODELS BASED ON {}'.format(return_best))
        if '/' in return_best:
            return_best_folder = return_best.replace('/', '_')
        else:
            return_best_folder = return_best
        out = os.path.join('relaxed_outputs',
                'best_{}'.format(return_best_folder))
        if not os.path.exists(out):
            os.mkdir(out)
        for name, group in records.groupby('origin'):
            print(name)
            row = group.loc[group[return_best].idxmin()]
            best = row['path']
            basename = row['filename']
            copyfile(best, os.path.join(out,
                basename))
        sys.exit()
    return records


def plot(records, x='origin', y='dG_separated'):
    xvals = records[x]
    yvals = records[y]

    fig, ax = plt.subplots()

    if x=='origin':
        yvals = []
        origin_yvals = []
        labels = sorted(list(set(xvals)))
        labels_str = [str(l) for l in labels]
        ind1 = []
        ind2 = []
        i = 0
        for label in labels:
            ind1.append(i + 1)
            ind2.append(i + 2)
            i += 2
            origin_path = os.path.join('relaxed_outputs', label +
                    '_origin')
            origin_record = None
            print(label)
            print(origin_path)
            if os.path.exists(origin_path):
                origin_record = load_pdbs(origin_path)
                if not origin_record.empty:
                    origin_yvals.append(origin_record[y].values)
                else:
                    origin_yvals.append([-3])
            else:
                origin_yvals.append([-3])
            record = records[records[x] == label]
            yvals.append(record[y].values)

        #width = 0.45
        #widthmodifier = 0.3
        #widthlist = widthmodifier * len(labels) * width/(len(labels))

        vplt1 = ax.violinplot(yvals,
                positions=ind1,
                #positions=[ind for ind in range(len(labels))],
                showmedians=True)

        #if not origin_record:
        vplt2 = ax.violinplot(origin_yvals, 
                positions=ind2,
                showmedians=True)
        facecolor = 'lightcoral'
        edgecolor='firebrick'
        for pc in vplt2['bodies']:
            pc.set_edgecolor(edgecolor)
            pc.set_facecolor(facecolor)
        for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
            vp = vplt2[partname]
            vp.set_edgecolor(edgecolor)

        suppress_text = False
        if not suppress_text:
            for j, v in enumerate(yvals):
                #print(j, v)
                pos = max(v)# + 0.05 * abs(max(v))
                #print(pos)
                ax.text(2 * j - 1, pos, labels[j].split('_')[0])

        ax.set_xlabel('Parent structure')
        ax.set_ylabel(y)
        #ax.set_xticklabels(labels_str)
    
    plt.show()

if __name__=='__main__':
    args = docopt.docopt(__doc__)
    x = args['-x']
    y = args['-y']
    filename = args['--filename']
    folder = 'relaxed_outputs'
    return_best = args['--best']
    records = load_pdbs(folder, return_best=return_best)
    plot(records, x=x, y=y)
