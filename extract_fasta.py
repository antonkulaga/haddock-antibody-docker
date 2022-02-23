#!/usr/bin/env python

from collections import OrderedDict

import Bio
import click
from Bio import SeqIO, SeqRecord
from beartype import beartype
from pycomfort.files import *


@beartype
def pdb_to_fasta(pdf_path: Path) -> OrderedDict:
    def map_chain(record: Bio.SeqRecord.SeqRecord):
        record.id = record.id.replace("????:", ":")
        record.id = record.id.replace("????:", ":")
        return record.id[1:], record.seq
    with pdf_path.open("r") as pdb_file:
        sequences = seq([s for s in SeqIO.parse(pdb_file, 'pdb-atom')]).map(map_chain)
    return OrderedDict(sequences)


@beartype
def write_fasta_pdb(pdb_path: Path, where: Path) -> Path:
    d = pdb_to_fasta(pdb_path)
    with where.open("w") as f:
        for k,v in d.items():
            f.write(f">:{k}\n")
            f.write(f"{str(v)}\n")
    return where

@beartype
def extract_fasta(pdb_path: Path, where: Path) -> Path:
    if pdb_path.is_dir():
        where.mkdir(exist_ok=True)
        with_ext(pdb_path, "pdb").for_each(lambda child: extract_fasta(child, where / (child.stem + ".fasta")))
        dirs(pdb_path).for_each(lambda child: extract_fasta(child, where / child.name))
        return where
    else:
        return write_fasta_pdb(pdb_path, where)


@click.command("extract_fasta")
@click.argument("pdb",  type=click.Path(exists=True))
@click.argument("where",  type=click.STRING)
def cli(pdb: str, where: str):
    return extract_fasta(Path(pdb), Path(where))

if __name__ == '__main__':
    cli()
