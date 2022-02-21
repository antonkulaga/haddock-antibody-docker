#!/usr/bin/env python
import sys
from typing import Union

import click
from pathlib import Path

import numpy as np
import pandas as pd
from functional import seq
from numpy import int8

from beartype import beartype
from Bio import SeqIO


script_folder = Path(__file__).parent.resolve()
print(script_folder)
sys.path.append(str((script_folder / 'haddock-antibody')))
sys.path.append(str((script_folder / 'haddock-tools')))

import run
from run import process_pdb, tidy_up
from restrain_bodies import read_structure, build_restraints, calc_euclidean, get_bodies, generate_tbl

calc_accessibility = __import__("calc-accessibility")

#rewrote in a proper way
def active_passive_to_ambig(active1: list, passive1: list, active2: list, passive2: list, segid1='A', segid2='B'):
    """Convert active and passive residues to Ambiguous Interaction Restraints

    Parameters
    ----------
    active1 : list
        List of active residue numbers of the first segid

    passive1 : list
        List of passive residue numbers of the first segid

    passive2 : list
        List of passive residue numbers of the second segid

    active2 : list
        List of active residue numbers of the second segid

    active2 : list
        List of passive residue numbers of the second segid

    segid1 : string
        Segid to use for the first model

    segid2 : string
        Segid to use for the second model

    """

    all1 = active1 + passive1
    all2 = active2 + passive2
    result = ""

    for resi1 in active1:
        result += 'assign (resi {:d} and segid {:s})\n'.format(resi1, segid1)
        result += '(\n'
        c = 0
        for resi2 in all2:
            result += '       (resi {:d} and segid {:s})\n'.format(resi2, segid2)
            c += 1
            if c != len(all2):
                result += '        or\n'
        result += ') 2.0 2.0 0.0\n'

    for resi2 in active2:
        result += 'assign (resi {:d} and segid {:s})\n'.format(resi2, segid2)
        result += '(\n'
        c = 0
        for resi1 in all1:
            result += '       (resi {:d} and segid {:s})\n'.format(resi1, segid1)
            c += 1
            if c != len(all1):
                result += '        or\n'

        result += ') 2.0 2.0 0.0\n'
    return result


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


@beartype
def restrain(antibody: Path) -> str:
    """
    generates restrains tbl for the antibody
    :param antibody:
    :return:
    """
    atom_lst = read_structure(str(antibody))
    restraints = build_restraints(get_bodies(atom_lst))
    return make_tbl(atom_lst, restraints)


@app.command("restrain")
@click.option('--antibody', type=click.Path(exists=True), help='antibody pdb')
def restrain_command(antibody: str) -> str:
    return restrain(Path(antibody))


@beartype
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


@app.command("run_param")
@click.option('--antibody', type=click.Path(exists=True), help='pdb')
@click.option('--antigen', type=click.Path(exists=True), help='pdb')
@click.option('--ambig', type=click.Path(exists=True), help='amibg')
@click.option('--unambig', type=click.Path(exists=True), help="unambig")
@click.option('--project', type=click.Path(exists=True), help="project")
@click.option("--haddock_dir", default="/data/sources/haddock2.4", help="folder where haddock is located")
@click.option("--n_comp", default=2, help="N_COMP")
@click.option("--run_number", default=1, help="run_number")
def run_params_command(antibody: str, antigen: str, ambig: str, unambig: str, project: str, haddock_dir: Path = Path("/data/sources/haddock2.4"), n_comp: Union[str, int] = 2, run_number: Union[int, str] =1) -> Path:
    #run_params(result_pdb, antigen_path, ambig, unambig, output_path.absolute(), Path(haddock_dir), n_comp, run_number)
    params =  run_params(Path(antibody), Path(antigen), Path(ambig), Path(unambig), Path(project), Path(haddock_dir), n_comp, run_number)
    path = Path("run.param")
    with path.open("w") as f:
        f.write(params)
    print(f"run.param are: \n {params}")
    return path


@beartype
def run_params(a: Path, b: Path, ambig: Path, unambig: Path, project: Path, haddock_dir: Path = Path("/data/sources/haddock2.4"),
               n_comp: Union[str, int] = 2, run_number: Union[int, str] =1):
    return f'''AMBIG_TBL={str(ambig.resolve())}
HADDOCK_DIR={haddock_dir.resolve()}
N_COMP={n_comp}
PDB_FILE1={str(a.resolve())}
PDB_FILE2={str(b.resolve())}
PROJECT_DIR={project.resolve()}
PROT_SEGID_1=A
PROT_SEGID_2=B
RUN_NUMBER={run_number}
UNAMBIG_TBL={str(unambig.resolve())}
    '''


@beartype
def run_params_in_project(project: Path, a: Path, b: Path, ambig: Path, unambig: Path, haddock_dir: Path = Path("/data/sources/haddock2.4"),
                          n_comp: Union[str, int] = 2, run_number: Union[int, str] =1):
    return run_params(project / a.name, project / b.name, project / ambig.name, project / unambig.name, project, haddock_dir, n_comp, run_number)


def pdb_2_fasta(pdb_file: Path, where: Path):
    print(f"extracting fasta from {pdb_file} to {where}")
    name = pdb_file.stem
    with pdb_file.open("r") as pdb_file:
        sequences = [record for record in SeqIO.parse(pdb_file, 'pdb-atom')]
    with where.open("w+") as fasta:
        for sequence in sequences:
            #for chain in sequence.
            #where.write_text(">" + sequence.id.replace("????:", name +":") + "\n"+str(sequence._seq))
            SeqIO.write(sequence, fasta, "fasta")

#for example start.py extract_fv
@app.command("extract_fv")
@click.argument("folder", type=click.Path(exists=True))
@click.option('--output', default="output", help='output folder to store results')
@click.option('--scheme', default="c", help="numbering scheme")
@click.option('--mode', default="all", help="also extracts fasta sequences")
@click.option('--chain', default="A", help="chain to extract active regions from")
@click.option('--rename', default=True, help="renaming")
def extract_fv(folder: str, output: str, scheme: str, mode: str, chain: str, rename: bool) -> Path:
    print(f"extract_fv for {folder}, the results will be extracted to {output}")
    pdb_path = Path(folder)
    output_path = Path(output)
    output_path.mkdir(exist_ok=True)
    print(f"{str(pdb_path)} is folder, processing all subfolders and pdb files inside of it!")
    for child in pdb_path.iterdir():
        output_subpath = output_path / child.name
        if child.is_dir() and not child.is_symlink() and any(child.iterdir()):
            print("FOLDER")
            output_subpath.mkdir(exist_ok=True)
            #run.process_folder(child, output_subpath, scheme, True, rename, True, chain, True)
            extract_fv(str(child), str(output_subpath), scheme, fasta, chain, rename)
        elif child.is_file() and "pdb" in child.suffix:
            print("FILE")
            try:
                print(f"TRYING {child}")
                pdb = child if mode == "just_fasta" else run.process_pdb(child, output_path, scheme, True, True, rename, chain, True)
                if mode != "no_fasta":
                    pdb_2_fasta(pdb, output_path / (pdb.stem + ".fasta"))
            except BaseException as e:
                print(f"pdb {child.name} FAILED, exception is {e}")
                continue
    return output_path


@app.command("start")
@click.option('--antibody', type=click.Path(exists=True), help='pdb file or a folder with pdb files to run protocol at, for example 4G6K.pdb (file) or my_antibodies (folder)')
@click.option('--antigen', type=click.Path(), help='pdb file with the antigen')
@click.option('--output', default="output", help='output folder to store results')
@click.option('--project', type=click.Path(exists=False), default=None, help = "project path to generate")
@click.option('--scheme', default="c", help="numbering scheme")
@click.option('--fvonly', default=True, help="use only fv region")
@click.option('--rename', default=True, help="renaming")
@click.option('--splitscfv', default=True, help="splitscfv")
@click.option('--chain', default="A", help="chain to extract active regions from")
@click.option('--delete_intermediate', default=False, help="Delete intermediate files")
@click.option('--cutoff', default=0.15, help="access") #0.4 in the tutorial
@click.option('--antigen_cutoff', default=0.15, help="access") #0.4 in the tutorial
@click.option("--haddock_dir", default="/data/sources/haddock2.4", help="folder where haddock is located")
@click.option("--n_comp", default =2, help = "N_COMP")
@click.option("--run_number", default = 1, help = "run_number")
def start(antibody: str, antigen: str, output: str,
          project: str, scheme: str, fvonly: bool, rename: bool,
        splitscfv: bool, chain: str,
        delete_intermediate: bool, cutoff: float, antigen_cutoff: float,
        haddock_dir: str, n_comp: str, run_number: int):
    print(f"calling run with {antibody} antibody and {antigen} antigen where output is {output}")
    antibody_path = Path(antibody)
    output_path = Path(output)
    output_path.mkdir(exist_ok=True)
    pdb_name = antibody_path.name
    result_pdb = run.process_pdb(antibody_path, output_path, scheme, fvonly, rename, splitscfv, chain, delete_intermediate)
    chains_restrain = restrain(result_pdb)
    unambig = (output_path / "antibody-unambig.tbl")
    unambig.write_text(chains_restrain)
    active = pd.read_csv(output_path / pdb_name.replace(".pdb", f"_active.txt"), header=None).values.tolist()[0]
    print(active)
    ares = access(str(result_pdb), cutoff)
    acc_residues = np.array(ares, dtype=int8)
    intersection = np.intersect1d(active, acc_residues)
    active_passive_antibody = output_path / pdb_name.replace(".pdb", f"_antibody_active_passive.txt")
    antibody_active_solv = intersection.tolist()
    with active_passive_antibody.open("w") as antibody_f:
        antibody_f.write(" ".join([str(i) for i in antibody_active_solv]) + "\n")
    antigen_path = (output_path / Path(antigen).name.replace(".pdb", f"_tidy.pdb")).resolve()
    tidy_up(Path(antigen), antigen_path)
    ares_antigen = access(str(result_pdb), antigen_cutoff)
    print("===========")
    ares_antigen_str = ' '.join([str(i) for i in ares_antigen])
    print("passive residues of the antigen are: " + "".join(ares_antigen_str))
    active_passive_antigen = output_path / antigen_path.name.replace(".pdb", f"_antigen_active_passive.txt")
    with active_passive_antigen.open("w") as antigen_f:
        antigen_f.write("\n" + ares_antigen_str)
    passive_active = active_passive_to_ambig(antibody_active_solv, [], [], ares_antigen)
    ambig = (output_path / "antibody-antigen-ambig.tbl")
    with ambig.open("w") as af:
        af.write(passive_active)
    if project is None:
        run_str = run_params(result_pdb, antigen_path, ambig, unambig, output_path.absolute(), Path(haddock_dir), n_comp, run_number)
    else:
        run_str = run_params_in_project(Path(project), result_pdb, antigen_path, ambig, unambig, Path(haddock_dir), n_comp, run_number)
    run_file = (output_path / "run.param")
    print(f"writing run parameters to {str(run_file)}")
    with run_file.open("w") as f:
        f.write(run_str)
    return run_file



if __name__ == '__main__':
    app()
