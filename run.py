#!/usr/bin/env python
import sys
import os
import click
from pathlib import Path

script_folder = Path(__file__).parent.resolve()
print(Path(__file__).parent.absolute())
sys.path.append(str((script_folder / 'HADDOCK-antibody-antigen')))
sys.path.append(str((script_folder / 'haddock-tools')))

from restrain_bodies import *
from cli import cli

def make_tbl(atom_lst: list, restraints: list) -> str:
    """
    Makes a list of TBL-formatted restraints.
    """
    result: str = ""
    for r in restraints:
        i, j = r
        atom_i, atom_j = atom_lst[i], atom_lst[j]
        dist_ij = calc_euclidean(atom_i[3], atom_j[3])

        tbl = "assign (segid {0[0]} and resi {0[1]} and name {0[2]}) ".format(atom_i)
        tbl += "(segid {0[0]} and resi {0[1]} and name {0[2]}) ".format(atom_j)
        tbl += "{0:3.3f} 0.0 0.0".format(dist_ij)
        result += tbl + "\n"
    return result

def prepare(antibody: Path, antigen: Path = None, output: Path = Path("./")) -> str:
    # Main logic
    print(f"antibody pdb is {antibody.name}, antigen pdb is {antigen.name} output folder is {str(output)}")
    atom_lst = [read_structure(antibody)]
    restraints = build_restraints(get_bodies(atom_lst))
    return generate_tbl(atom_lst, restraints)

'''
@click.command()
@click.option('--pdb', help='pdb file or a folder with pdb files to run protocol at, for example 4G6K.pdb (file) or my_antibodies (folder)')
@click.option('--output', default="output", help='output folder to store results')
@click.option('--scheme', default="c", help="numbering scheme")
@click.option('--fvonly', default=True, help="use only fv region")
@click.option('--rename', default=True, help="renaming")
@click.option('--splitscfv', default=True, help="splitscfv")
@click.option('--chain', default="A", help="chain to extract active regions from")
@click.option('--delete_intermediate', default=True, help="Delete intermediate files")
def cli(pdb: str, output: str, scheme: str, fvonly: bool, rename: bool, splitscfv: bool, chain: str, delete_intermediate: bool):
'''


@click.command()
@click.option('--antibody', type=click.Path(exists=True), help='pdb file or a folder with pdb files to run protocol at, for example 4G6K.pdb (file) or my_antibodies (folder)')
@click.option('--antigen', type=click.Path(exists=False), help='pdb file with the antigen')
@click.option('--output', default="output", help='output folder to store results')
@click.option('--scheme', default="c", help="numbering scheme")
@click.option('--fvonly', default=True, help="use only fv region")
@click.option('--rename', default=True, help="renaming")
@click.option('--splitscfv', default=True, help="splitscfv")
@click.option('--chain', default="A", help="chain to extract active regions from")
@click.option('--delete_intermediate', default=True, help="Delete intermediate files")
def run(antibody: str, antigen: str, output: str, scheme: str, fvonly: bool, rename: bool, splitscfv: bool, chain: str, delete_intermediate: bool):
    antibody_path = Path(antibody)
    output_path = Path(output)
    tbl = prepare(antibody_path)
    print("TBL IS:")
    print(tbl)
    cli(antibody, output, scheme, fvonly, rename, splitscfv, chain, delete_intermediate)

if __name__ == '__main__':
    run()
