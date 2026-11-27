#!/usr/bin/env python3
"""Additional native representation, phase, exclusion and baseline regressions."""
from pathlib import Path
import argparse
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('binary', type=Path)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('work', type=Path)
    args = parser.parse_args()
    root = args.work.resolve()
    root.mkdir(parents=True, exist_ok=True)
    binary, baseline = args.binary.resolve(), args.baseline.resolve()
    par = root / 'par'
    par.mkdir(exist_ok=True)
    for chrom in ('X', 'Y'):
        (par / f'chr{chrom}_PAR.bed').write_text(f'chr{chrom}\t0\t100\nchr{chrom}\t900\t1000\n')
        (par / f'chr{chrom}_nonPAR.bed').write_text(f'chr{chrom}\t100\t400\nchr{chrom}\t500\t900\n')
    reference = root / 'ref.fa'
    reference.write_text(''.join(f'>{c}\n' + 'A' * 1000 + '\n' for c in ('chr1', 'chrX', 'chrY')))
    header = '##fileformat=VCFv4.2\n' + ''.join(f'##contig=<ID={c},length=1000>\n' for c in ('chr1', 'chrX', 'chrY'))
    header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    roles = ('mother', 'father', 'child')
    # pos, ref, alt, GT; records absent from irrelevant VCF do not remove chr.
    cases = {
        'mnp': ('male', 'chrX', [[(301, 'AA', 'CG', '0/1')], [], [(301, 'A', 'C', '1'), (302, 'A', 'G', '1')]], True),
        'female_mnp': ('female', 'chrX', [[(301, 'AA', 'CG', '0/1')], [(301, 'AA', 'TT', '1')], [(301, 'A', 'C,T', '1|2'), (302, 'A', 'G,T', '1|2')]], True),
        'phase_no_ps_bad': ('female', 'chrX', [[(301, 'A', 'C', '1|1'), (302, 'A', 'G', '0|0')], [(301, 'A', 'C', '0'), (302, 'A', 'G', '1')], [(301, 'A', 'C', '0|1'), (302, 'A', 'G', '0|1')]], False),
        'phase_no_ps_good': ('female', 'chrX', [[(301, 'A', 'C', '1|1'), (302, 'A', 'G', '0|0')], [(301, 'A', 'C', '0'), (302, 'A', 'G', '1')], [(301, 'A', 'C', '1|0'), (302, 'A', 'G', '0|1')]], True),
        'overlap_trim_end': ('male', 'chrX', [[(301, 'AAA', 'AA', '1/1'), (303, 'A', 'C', '1/1')], [], [(301, 'AAA', 'AC', '1')]], True),
        'overlap_trim_begin': ('male', 'chrX', [[(301, 'A', 'C', '1/1'), (301, 'AAA', 'AA', '1/1')], [], [(301, 'AAA', 'CA', '1')]], True),
        'shifted': ('male', 'chrY', [[], [(301, 'AA', 'A', '1')], [(305, 'AA', 'A', '1')]], True),
    }
    for name, (sex, chrom, samples, consistent) in cases.items():
        out = root / name
        out.mkdir(exist_ok=True)
        command = [str(binary), 'mendelian', '--child-sex', sex, '--par-dir', str(par), '-ref', str(reference), '-outDir', str(out), '-out-prefix', 'test']
        for role, records in zip(roles, samples):
            file = out / f'{role}.vcf'
            # Irrelevant empty role deliberately has no X/Y contig declaration.
            role_header = header if records else header.replace(f'##contig=<ID={chrom},length=1000>\n', '')
            file.write_text(role_header + f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{role}\n' + ''.join(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT\t{gt}\n' for pos, ref, alt, gt in records))
            command.extend([f'-{role}', str(file)])
        result = subprocess.run(command, text=True, capture_output=True)
        (out / 'run.log').write_text(result.stdout + result.stderr)
        assert result.returncode == 0, (name, result.stderr[-2000:])
        rows = [line for line in (out / 'test_child_xy.vcf').read_text().splitlines() if not line.startswith('#')]
        has_violation = any('XY_STATUS=violation' in row for row in rows)
        assert has_violation != consistent, (name, rows)
        assert not any('unsupported' in row for row in rows), (name, rows)
    # One unchanged chr1 fixture drives both old mode and explicit-policy mode.
    inputs = root / 'autosome-inputs'
    inputs.mkdir(exist_ok=True)
    for role in roles:
        rows = [(101, '0/1'), (201, '1/1'), (301, '0/1'), (401, '1/1')]
        (inputs / f'{role}.vcf').write_text(header + f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{role}\n' + ''.join(f'chr1\t{pos}\t.\tA\tC\t60\tPASS\t.\tGT\t{gt}\n' for pos, gt in rows))
    outputs = []
    for name, executable, explicit in [('baseline', baseline, False), ('legacy', binary, False), ('explicit', binary, True)]:
        out = root / name
        out.mkdir(exist_ok=True)
        command = [str(executable), 'mendelian', '-ref', str(reference), '-outDir', str(out), '-out-prefix', 'test', '--autosome-only']
        for role in roles:
            command.extend([f'-{role}', str(inputs / f'{role}.vcf')])
        if explicit:
            command.extend(['--child-sex', 'unknown', '--par-dir', str(par)])
        result = subprocess.run(command, text=True, capture_output=True)
        (out / 'run.log').write_text(result.stdout + result.stderr)
        assert result.returncode == 0, (name, result.stderr[-2000:])
        outputs.append((out / 'test_trio.vcf').read_text())
    assert outputs[0] == outputs[1] == outputs[2], 'Autosome baseline differs'
    print('native representation/phase/overlap: 7 passed; baseline/legacy/explicit autosome byte-identical')


if __name__ == '__main__':
    main()
