# License: GNU Affero General Public License v3 or later
from __future__ import absolute_import

from builtins import str
import time
import os

from .lib import run_cmd
from .log import return_logger
import subprocess
from sys import argv

logger = return_logger('RRE-finder', False)

def run_rre(infile, rre_path, cores, proj_name):
    cmd = ['python3', 'RRE.py', '-i', infile, '-o', rre_path, '-c', str(cores), '-t', 'fasta', proj_name]
    run_cmd(cmd)


def write_fasta_all(rre_path, operons):
    outf = os.path.join(rre_path, 'all_proteins.fasta')
    if not os.path.isfile(outf):
        fa_text = '\n'.join([gene.write_fasta() for gene in operons.itergenes()])
        with open(outf, 'w') as out:
             out.write(fa_text)
    return outf


def read_rre_data(rre_results):
    with open(rre_results) as inf:
        header = inf.readline()
        data = [line.strip().split('\t') for line in inf.readlines()]
    return data


def parse_rre_resubmit(rre_out, operons):
    rre_data = read_rre_data(rre_out)
    logger.info('Found %s RREs in exploratory mode' % len(rre_data))
    for hit in rre_data:
        gene_name = hit[0]
        hit_id = hit[3]
        loc = '-'.join(hit[10:])
        rre_data = hit[4:10] + [loc] # [Prob, E-val, P-val, Score, SS, Columns, loc]
        operons[gene_name].RRE_hit = True
        operons[gene_name].RRE_data = {hit_id: rre_data}
        operons[gene_name].best_hit = hit_id
    return operons

  
def parse_rre_precision(rre_out, operons):
    rre_data = read_rre_data(rre_out)
    logger.info('Found %s RREs in precision mode' % len(set([r[0] for r in rre_data])))
    for hit in rre_data:
        gene_name = hit[0]
        hit_id = hit[3]
        loc = '-'.join(hit[6:])
        rre_data = ['n/a',hit[4],'n/a',hit[5],'n/a','n/a',loc] # [Prob, E-val, P-val, Score, SS, Columns, loc]
        # if another RRE was already found, set lowest E-value as best hit
        if operons[gene_name].RRE_hit:
            operons[gene_name].RRE_data[hit_id] = rre_data
            if rre_data[1] < operons[gene_name].RRE_data[operons[gene_name].best_hit][1]:
                operons[gene_name].best_hit = hit_id
        else:
            operons[gene_name].RRE_hit = True
            operons[gene_name].RRE_data = {hit_id: rre_data}
            operons[gene_name].best_hit = hit_id
    return operons


def add_nr_rre(operons):
    for operon in operons:
        nr = 0
        for _ in operon.itergenes(RRE_hit=True):
            nr += 1
        operon.nr_RRE = nr
    return operons


def main(path, rre_path, genome_dict, operons, settings):
    rre_project = 'RRE_all'
    rre_resubmit = os.path.join(rre_path, rre_project, rre_project+'_rrefinder_RRE_resubmit_results.txt')
    rre_precision = os.path.join(rre_path, rre_project, rre_project+'_rrefam_results.txt')
    if not os.path.isfile(rre_resubmit) or not os.path.isfile(rre_precision):
        fasta_all = write_fasta_all(rre_path, operons)
        run_rre(fasta_all, rre_path, settings['cores'], rre_project)

    for gene in operons.itergenes():
        gene.RRE_hit = False
    operons = parse_rre_resubmit(rre_resubmit, operons)
    operons = parse_rre_precision(rre_precision, operons)
    operons = add_nr_rre(operons)
    return operons
