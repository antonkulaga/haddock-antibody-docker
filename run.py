#!/usr/bin/env python
import sys
import os
import click
from pathlib import Path
import importlib

import numpy as np
import pandas as pd
from functional import seq
from numpy import int8

script_folder = Path(__file__).parent.resolve()
print(Path(__file__).parent.absolute())
sys.path.append(str((script_folder / 'HADDOCK-antibody-antigen')))
sys.path.append(str((script_folder / 'haddock-tools')))
import cli
from restrain_bodies import read_structure, build_restraints, calc_euclidean, get_bodies, generate_tbl
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
def app():
    return True

def restrain(antibody: Path) -> str:
    atom_lst = read_structure(str(antibody))
    restraints = build_restraints(get_bodies(atom_lst))
    return make_tbl(atom_lst, restraints)


@app.command("restrain")
@click.option('--antibody', type=click.Path(exists=True), help='')
def restrain_command(antibody: str) -> str:
    return restrain(Path(antibody))


def apply_cutoff(access_data, cutoff: float):
    """Apply a cutoff to the sidechain relative accessibility and display on stdout"""
    print(f'Applying cutoff to side_chain_rel - {cutoff}')
    def chain_result(chain) -> list:
        result_list = []
        for res in access_data[chain]:
            sidechain_relative_accessibility = access_data[chain][res]['side_chain_rel']
            if sidechain_relative_accessibility >= cutoff * 100:
                result_list.append(res)
        result = list(set(result_list))
        result.sort()
        return result
    return seq(access_data).map(chain_result).to_list()


def access(pdb: str, cutoff: float):
    access_dic = calc_accessibility.get_accessibility(pdb)
    return apply_cutoff(access_dic, cutoff)[0]


@app.command("access")
@click.option('--pdb', type=click.Path(exists=True), help='pdb file to compute accessibility for')
@click.option('--cutoff', default=0.15, help="access") #0.4 in the tutorial
def access_command(pdb: str, cutoff: float = 0.15):
    return access(pdb, cutoff)


@app.command()
@click.option('--antibody', type=click.Path(exists=True), help='pdb file or a folder with pdb files to run protocol at, for example 4G6K.pdb (file) or my_antibodies (folder)')
@click.option('--antigen', type=click.Path(), help='pdb file with the antigen')
@click.option('--output', default="output", help='output folder to store results')
@click.option('--scheme', default="c", help="numbering scheme")
@click.option('--fvonly', default=True, help="use only fv region")
@click.option('--rename', default=True, help="renaming")
@click.option('--splitscfv', default=True, help="splitscfv")
@click.option('--chain', default="A", help="chain to extract active regions from")
@click.option('--delete_intermediate', default=False, help="Delete intermediate files")
@click.option('--cutoff', default=0.15, help="access") #0.4 in the tutorial
def run(antibody: str, antigen: str, output: str,
        scheme: str, fvonly: bool, rename: bool,
        splitscfv: bool, chain: str,
        delete_intermediate: bool, cutoff: float):
    print(f"calling run with {antibody} antibody and {antigen} antigen where output is {output}")
    antibody_path = Path(antibody)
    output_path = Path(output)
    output_path.mkdir(exist_ok=True)
    pdb_name = antibody_path.name
    tidy_pdb = (output_path / antibody_path.name.replace(".pdb", f"_HADDOCK_tidy.pdb")).resolve()
    cli.process_pdb(antibody_path, output_path, scheme, fvonly, rename, splitscfv, chain, delete_intermediate)
    chains_restrain = restrain(tidy_pdb)
    (output_path / "antibody-unambig.tbl").write_text(chains_restrain)
    active = pd.read_csv(output_path / pdb_name.replace(".pdb", f"_active_sites.txt"), header=None).values.tolist()[0]
    print(active)
    print(type(active))
    print("====")
    ares = access(str(tidy_pdb), cutoff)
    acc_residues = np.array(ares, dtype=int8)
    antibody_final = np.intersect1d(active, acc_residues)
    print(antibody_final)

if __name__ == '__main__':
    app()
