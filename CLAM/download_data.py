import os
import sys
import subprocess


def parser(args):
    """DOCSTRING
    Args
    Returns
    """
    try:
        genome = args.genome
        download_genome(genome)
    except KeyboardInterrupt():
        sys.exit(0)


def download_genome(genome):
    curr_dir = os.path.abspath('.')
    
    admin = (os.getuid() == 0)
    cmd = []
    home = os.environ['HOME']
    if admin:
        profile = '/etc/profile'
    else:
        profile = '{home}/.bashrc'.format(home=home)
    
    if not os.path.isdir('{home}/.clam_data'.format(home=home)):
        os.mkdir('{home}/.clam_data'.format(home=home))
    os.chdir('{home}/.clam_data'.format(home=home))

    if 'CLAM_DAT' not in os.environ or not os.environ['CLAM_DAT'] == '{home}/.clam_data'.format(home=home):
        cmd.append('echo "export CLAM_DAT=\'{clam_data}\'" >> {profile}'.format(
            clam_data=os.path.abspath('.'), profile=profile))
        cmd.append('source {profile}'.format(profile=profile))
        os.environ['CLAM_DAT'] = os.path.abspath('.')

    if not check_genome_data(genome):
        cmd.append('chmod -R 755 {home}/.clam_data'.format(home=home))
        cmd.append(
            'wget https://raw.githubusercontent.com/wkdeng/clam_data/master/{genome}.zip'.format(genome=genome))
        cmd.append('unzip -o {genome}.zip'.format(genome=genome))
        cmd.append('rm  {genome}.zip'.format(genome=genome))
        for item in cmd:
            subprocess.call(item, shell=True, executable='/bin/bash')
        print('Download finished')    
    os.chdir(curr_dir)

def check_genome_data(genome):
    if not os.path.isdir(os.environ['CLAM_DAT'] + '/' + genome):
        return False
    if not os.path.exists(os.environ['CLAM_DAT'] + '/' + genome + '/3UTRs.bed'):
        return False
    if not os.path.exists(os.environ['CLAM_DAT'] + '/' + genome + '/5UTRs.bed'):
        return False
    if not os.path.exists(os.environ['CLAM_DAT'] + '/' + genome + '/cds.bed'):
        return False
    if not os.path.exists(os.environ['CLAM_DAT'] + '/' + genome + '/exons.bed'):
        return False
    if not os.path.exists(os.environ['CLAM_DAT'] + '/' + genome + '/introns.bed'):
        return False
    if not os.path.exists(os.environ['CLAM_DAT'] + '/' + genome + '/proximal200_intron.bed'):
        return False
    if not os.path.exists(os.environ['CLAM_DAT'] + '/' + genome + '/proximal500_intron.bed'):
        return False
    return True

if __name__ == '__main__':
    download_genome('hg38')
