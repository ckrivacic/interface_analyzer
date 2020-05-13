'''
Please put path to pdbredo as $PDBREDO in your bashrc
Usage:
    submit.py

'''
import sys, os, subprocess, re, glob
import docopt
#task_len = 500

def file_len(fname):
    with open(fname) as f:
        for i,l in enumerate(f):
            pass
    return i + 1

def submit(**params):
    from klab import cluster, process
    import json
    args = docopt.docopt(__doc__)
    #script_path = os.path.expanduser(args['<script>'])
    logdir = 'rosettascripts_logs'
    csv_path = '2020-03-18_Krogan_SARSCoV2_27baits.txt'
    if not os.path.exists(logdir):
        os.makedirs(logdir, exist_ok=True)

    pdb_path = os.path.expanduser('~/sars/figures/pymol_sessions/')
    num_tasks = len(glob.glob(pdb_path + '/*/*/*.pdb')) * 100

    max_runtime = params.get('max_runtime','24:00:00')
    max_memory = params.get('max_memory','4G')

    xml = '/wynton/home/kortemme/krivacic/sars/interface_analyzer/relax.xml'

    python = '/wynton/home/kortemme/krivacic/software/anaconda3/bin/python3'
    script_path = os.path.expanduser('~/sars/interface_analyzer/run_rosettascripts.py')

    qsub_command = 'qsub', '-h', '-cwd',
    qsub_command += '-b',
    qsub_command += 'y',
    qsub_command += '-o', logdir,
    #qsub_command += '-e', error_directory,
    qsub_command += '-j', 'y',
    qsub_command += '-t', '1-{0}'.format(num_tasks),
    qsub_command += '-l', 'h_rt={0}'.format(max_runtime),
    qsub_command += '-l', 'mem_free={0}'.format(max_memory),
    qsub_command += python,
    qsub_command += script_path,
    qsub_command += xml,
    qsub_command += '--parent_dir', pdb_path,
    print(qsub_command)

    status = process.check_output(qsub_command, stderr=subprocess.STDOUT).decode('utf-8')

    status_pattern = re.compile(r'Your job-array (\d+).[0-9:-]+ \(".*"\) has been submitted')
    status_match = status_pattern.match(status)

    if not status_match:
        print(status)
        sys.exit()

    # Figure out the job id, then make a params file for it.
    job_id = status_match.group(1)

    qrls_command = 'qrls', job_id
    process.check_output(qrls_command)
    print(status)

if __name__ == '__main__':
    submit()
