import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem
import os, subprocess

os.environ['OE_LICENSE'] = '/db1/OpenEye/license/oe_license.txt'
if not os.path.isfile(os.environ['OE_LICENSE']):
    os.environ['OE_LICENSE'] = '/db1/OpenEye/oe_license.txt'

BROOD_EXE = '/db1/OpenEye/ubuntu18/bin/brood'
DB_CONFIG = '-db /db1/brood/chembl20'
EXTRA_PARAMS = '-minMolWt 0.0 -minpsa 0 -minHvyAtom 0'

default_DB = (
        '/db1/brood/chembl20',
        '/db1/brood/gostarDB',
        '/db1/brood/surechemblDB/part1',
        '/db1/brood/surechemblDB/part2',
        '/db1/brood/surechemblDB/part3',
        '/db1/brood/surechemblDB/part4',
        '/db1/brood/surechemblDB/part5',
    )
dblist = ','.join(default_DB)

def generateFragments(core_mol, scaffold_mol, logfile):
    frag_to_be_replaced = AllChem.ReplaceCore(core_mol, scaffold_mol)
    frag_to_be_replaced = Chem.MolToSmiles(frag_to_be_replaced).split('.')

    logfile.write('Fragments:\n')
    for i, frag in enumerate(frag_to_be_replaced):
        logfile.write(f' {i:3d}: {frag}\n')

    return frag_to_be_replaced


def showDBInfo(db_list, logfile):
    logfile.write("DB List:\n")
    for i, db in enumerate(db_list):
        logfile.write(f' {i:3d}: {db}\n')


def genFileFromSmilesText(smiles, filename):
    os.system(f"echo '{smiles}' > {filename}")


def runBrood(queryMol, queryFrag, db, prefix, nproc, other_params, logfile):
    cmd = [BROOD_EXE]
    cmd += ['-db', db]
    cmd += ['-queryMol', queryMol, '-queryFrag', queryFrag, '-prefix', prefix]
    cmd += EXTRA_PARAMS.split()
    if nproc > 1:
        cmd += ['-mpi_np', str(nproc)]
    if other_params is not None:
        cmd += other_params.split()

    logfile.write(f'{cmd}\n')
    subprocess.run(cmd)


def genTitleFile(csvfile):
    with open(csvfile) as rfile:
        title_str = rfile.readline()
        for l, line in enumerate(rfile):
            if l > 0:
                break

        if l > 0:
            with open('title.csv', 'w') as titlefile:
                titlefile.write(title_str)
            return False
        else:
            return True


def appendCSVFile(base_csv, new_csv, header=2):
    with open(new_csv) as newfile:
        for i, line in enumerate(newfile):
            if i >= header:
                base_csv.write(line)


def replaceMolFragWithBrood(core, scaffold, db_list=default_DB, nproc=1, other_params=None):
    with open('log.dat', 'w', -1) as logfile, open('res.csv', 'w') as csvfile:

        genFileFromSmilesText(core, 'core.smi')
        core_mol = Chem.MolFromSmiles(core)
        scaffold_mol = Chem.MolFromSmiles(scaffold)

        frag_to_be_replaced = generateFragments(core_mol, scaffold_mol, logfile)        
        showDBInfo(db_list, logfile)

        need_title = True

        for i, frag in enumerate(frag_to_be_replaced):
            fragfile = f'frag{i:03d}.smi'
            genFileFromSmilesText(frag, fragfile)

            for j, db in enumerate(db_list):
                prefix = f'frag{i:03d}_db{j:03d}'
                temp_csv = f'{prefix}.csv'
                runBrood('core.smi', fragfile, db, prefix, nproc, other_params, logfile)
                appendCSVFile(csvfile, temp_csv)
                if need_title:
                    need_title = genTitleFile(temp_csv)

    if not need_title:
        os.system('cat title.csv res.csv > res_wtitle.csv')
        os.system('mv res_wtitle.csv res.csv')
        os.system('rm title.csv')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='brood workflow')
    parser.add_argument('core_smiles', help='Smiles string for core molecule')
    parser.add_argument('scaffold_smiles', help='Smiles string for scaffold')
    parser.add_argument('--db', dest='brood_db', default=dblist, help=f'brood DB (separated by commas, default: {dblist})')
    parser.add_argument('--nproc', dest='nproc', default=1, type=int, help='number of processors per brood job (default: 1)')
    parser.add_argument('--other_params', dest='other_params', default=None, help='Brood parameters to change (default: None)')
    args = parser.parse_args()

    replaceMolFragWithBrood(sys.argv[1], sys.argv[2], args.brood_db.split(','), args.nproc, args.other_params)
