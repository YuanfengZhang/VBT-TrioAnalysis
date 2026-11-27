#!/usr/bin/env python3
"""Synthetic end-to-end native XY tests; no biological inputs or network."""
from pathlib import Path
import argparse
import gzip
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('binary', type=Path)
    parser.add_argument('work', type=Path)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    root = args.work.resolve()
    binary = args.binary.resolve()
    par = root / 'par'
    par.mkdir(exist_ok=True)
    for chrom in ('X', 'Y'):
        (par / f'chr{chrom}_PAR.bed').write_text(f'chr{chrom}\t0\t100\nchr{chrom}\t900\t1000\n')
        (par / f'chr{chrom}_nonPAR.bed').write_text(f'chr{chrom}\t100\t400\nchr{chrom}\t500\t900\n')
    reference = root / 'ref.fa'
    reference.write_text(''.join(f'>{c}\n' + 'A' * 1000 + '\n' for c in ('chr1', 'chrX', 'chrY')))
    roles = ('mother', 'father', 'child')
    header = '##fileformat=VCFv4.2\n' + ''.join(f'##contig=<ID={c},length=1000>\n' for c in ('chr1', 'chrX', 'chrY'))
    header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    cases = [
        ('male_x', 'male', False, 'chrX', 301, 'A', 'C', ('0/1', '.', '1'), 'consistent', 'pass'),
        ('male_y', 'male', False, 'chrY', 301, 'A', 'C', ('.', '1', '1'), 'consistent', 'pass'),
        ('female_x', 'female', False, 'chrX', 301, 'A', 'C', ('0/1', '1', '0/1'), 'consistent', 'pass'),
        ('female_x_bad', 'female', False, 'chrX', 301, 'A', 'C', ('1/1', '1', '0/1'), 'violation', 'mendelian_violation'),
        ('female_y', 'female', False, 'chrY', 301, 'A', 'C', ('.', '1', '1'), 'unsupported', 'female_y_nonpar'),
        ('unknown', 'unknown', False, 'chrX', 301, 'A', 'C', ('0/1', '1', '1'), 'not_evaluable', 'unknown_child_sex'),
        ('x_par', 'male', True, 'chrX', 31, 'A', 'C', ('0/1', '0/1', '1/1'), 'consistent', 'pass'),
        ('x_par_unmasked', 'male', False, 'chrX', 31, 'A', 'C', ('0/1', '0/1', '1/1'), 'not_evaluable', 'par_requires_ypar_masked'),
        ('y_par', 'male', True, 'chrY', 31, 'A', 'C', ('0/1', '0/1', '1/1'), 'not_evaluable', 'y_par_excluded'),
        ('pseudo', 'male', False, 'chrX', 301, 'A', 'C', ('0/1', '.', '1/1'), 'consistent', 'pass'),
        ('hetero', 'male', False, 'chrX', 301, 'A', 'C', ('0/1', '.', '0/1'), 'not_evaluable', 'haploid_heterozygous'),
        ('missing', 'male', False, 'chrX', 301, 'A', 'C', ('./.', '.', '1'), 'not_evaluable', 'required_gt_missing'),
        ('ref', 'male', False, 'chrX', 301, 'A', 'C', ('0/0', '1', '0'), 'not_evaluable', 'all_required_reference'),
        ('gap', 'male', False, 'chrX', 451, 'A', 'C', ('0/1', '.', '1'), 'not_evaluable', 'unclassified_region'),
        ('boundary', 'male', False, 'chrX', 399, 'AAA', 'A', ('0/1', '.', '1'), 'not_evaluable', 'boundary_crossing'),
        ('triploid', 'male', False, 'chrX', 301, 'A', 'C', ('0/1', '.', '1/1/1'), 'not_evaluable', 'unexpected_ploidy'),
        ('insert_anchor', 'male', False, 'chrX', 400, 'A', 'AC', ('0/1', '.', '1'), 'consistent', 'pass'),
        ('symbolic', 'male', False, 'chrX', 301, 'A', '<DEL>', ('0/1', '.', '1'), 'unsupported', 'unsupported_variant'),
    ]
    for name, sex, masked, chrom, pos, ref, alt, gts, status, reason in cases:
        out = root / name
        out.mkdir(exist_ok=True)
        for role, gt in zip(roles, gts):
            (out / f'{role}.vcf').write_text(header + f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{role}\n' + f'{chrom}\t{pos}\t{name}\t{ref}\t{alt}\t60\tPASS\t.\tGT\t{gt}\n')
        command = [str(binary), 'mendelian', '--child-sex', sex, '--par-dir', str(par), '-ref', str(reference), '-outDir', str(out), '-out-prefix', 'test']
        for role in roles:
            command.extend([f'-{role}', str(out / f'{role}.vcf')])
        if masked:
            command.append('--ypar-masked')
        result = subprocess.run(command, text=True, capture_output=True)
        (out / 'run.log').write_text(result.stdout + result.stderr)
        assert result.returncode == 0, (name, result.returncode, result.stderr[-2000:])
        lines = [line for line in (out / 'test_child_xy.vcf').read_text().splitlines() if not line.startswith('#')]
        assert len(lines) == 1, (name, lines)
        columns = lines[0].split('\t')
        info = dict(item.split('=', 1) for item in columns[7].split(';'))
        assert (info['XY_STATUS'], info['XY_REASON']) == (status, reason), (name, info)
        assert info['XY_SOURCE_ID'] == 'child:1'
        assert columns[9] == gts[2], (name, columns[9], gts[2])
    print(f'native XY CLI: {len(cases)} cases passed')


if __name__ == '__main__':
    main()
