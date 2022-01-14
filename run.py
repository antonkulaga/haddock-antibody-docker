#!/usr/bin/env python
import sys
import os
import click
from pathlib import Path
import importlib

script_folder = Path(__file__).parent.resolve()
print(Path(__file__).parent.absolute())
sys.path.append(str((script_folder / 'HADDOCK-antibody-antigen')))
sys.path.append(str((script_folder / 'haddock-tools')))

from restrain_bodies import *
calc_accessibility = __import__("calc-accessibility")

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


@click.group(chain=True)
def cli():
    return True

@cli.command("restrain")
@click.option('--antibody', type=click.Path(exists=True), help='')
def restrain(antibody: Path) -> str:
    atom_lst = [read_structure(antibody)]
    restraints = build_restraints(get_bodies(atom_lst))
    return generate_tbl(atom_lst, restraints)


def apply_cutoff(access_data, cutoff: float):
    """Apply a cutoff to the sidechain relative accessibility and display on stdout"""
    print(f'Applying cutoff to side_chain_rel - {cutoff}')
    for chain in access_data:
        result_list = []
        for res in access_data[chain]:
            sidechain_relative_accessibility = access_data[chain][res]['side_chain_rel']
            if sidechain_relative_accessibility >= cutoff * 100:
                result_list.append(res)
        result_list = list(set(result_list))
        result_list.sort()
        result_str = ','.join(map(str, result_list))
        print(f'Chain {chain} - {result_str}')

@cli.command("access")
@click.option('--pdb', type=click.Path(exists=True), help='pdb file to compute accessibility for')
@click.option('--cutoff', default=0.15, help="access") #0.4 in the tutorial
def access(pdb: str, cutoff: float):
    access_dic = calc_accessibility.get_accessibility(pdb)
    apply_cutoff(access_dic, cutoff)


@cli.command()
@click.option('--antibody', type=click.Path(exists=True), help='pdb file or a folder with pdb files to run protocol at, for example 4G6K.pdb (file) or my_antibodies (folder)')
@click.option('--antigen', type=click.Path(), help='pdb file with the antigen')
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
    chains_restrain = restrain(antibody_path)
    (output_path / "antibody-unambig.tbl").write_text(chains_restrain)
    cli(antibody, output, scheme, fvonly, rename, splitscfv, chain, delete_intermediate)


if __name__ == '__main__':
    cli()
