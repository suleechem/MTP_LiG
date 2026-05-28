import os
import re
import subprocess
from pathlib import Path


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip()


# Process directories: no letters, not containing "00"
prc = sorted(d for d in os.listdir('.') if not re.search(r'[a-zA-Z]', d) and '00' not in d)

with open('exec_files/to-relax.cfg') as f:
    config = str(sum(1 for line in f if 'BEGIN' in line))

HDR = "{:>12s} \t {:>12s} \t {:>12s} \t {:>12s} \t {:<12s} \t {:>12s} \t {:>12s}"
ROW = "{:>12s} \t {:>12s} \t {:>12s} \t {:>12s} \t {:<12s} \t {:>12s}\t {:>12s}"
print(HDR.format("Process", "Config", "Relaxed", "Preselected", "->Selected", "Computed", "Trained"))
print(HDR.format("", "", "", "Conf-Relax", "", "Tot.OPT.step", "Train+Compt"))

last = len(prc)
graden = []
grade = []

with open("xerrlist.txt", 'w') as fe:
    for i, proc in enumerate(prc):
        if not Path(f'{proc}/B-preselected.cfg').is_file():
            print(f"Skip process {proc}")
            last = len(prc) - 1
            break

        relax_path = f'{proc}/B-relax/relax-output.txt'
        is_last_two = i >= len(prc) - 2

        if is_last_two:
            fe.write("=" * 10 + "process " + proc + "=" * 10 + "\n")

        with open(relax_path) as f:
            lines = f.readlines()

        for a, line in enumerate(lines):
            if 'MV-grade' in line:
                mv = line.split()[7:9]
                pline = lines[a - 1].split() if a > 0 else []
                if is_last_two:
                    val = mv[-1].replace(")", "")
                    suffix = f" {pline[-1].replace('Breaking', '')}" if pline else ""
                    fe.write(f"Grade exceed {val}{suffix}\n")
            elif 'linesearch' in line and (a == 0 or 'linesearch' not in lines[a - 1]):
                if is_last_two:
                    pline = lines[a - 1].split() if a > 0 else []
                    if pline:
                        fe.write(f"Breaking linesearch {pline[-1]}\n")

        # Get top-2 MV-grade line numbers and values in a single pipeline
        grade_lines = sh(
            f"grep -n MV-grade {relax_path} | "
            "sed -e 's/)//g' -e 's/:/ /g' -e 's/(//g' | "
            "sort -g -k 10 | tail -2"
        ).splitlines()
        graden = [ln.split()[0] for ln in grade_lines if ln]
        grade = [p for ln in grade_lines for p in ln.split()[8:10]]

        cfg_files = sorted(Path(proc).glob('*.cfg'))
        lst = len(cfg_files)
        cfg = [str(sum(1 for l in open(f) if 'BEGIN' in l)) for f in cfg_files]
        cfg += [" - "] * (6 - lst)

        if lst == 6:
            row = ROW.format(proc, config, cfg[4], cfg[0], cfg[1], cfg[2], cfg[3])
        elif lst == 3:
            row = ROW.format(proc, config, cfg[2], cfg[0], cfg[1], cfg[4], cfg[3])
        elif lst == 1:
            row = ROW.format(proc, config, cfg[4], cfg[4], cfg[4], cfg[4], "0")
        elif lst == 4:
            row = ROW.format(proc, config, cfg[2], cfg[0], cfg[1], "= Done =", cfg[3])
        else:
            continue
        print(f"{row} {graden}{grade}")

print()

with open("xerrlist.txt") as f:
    xerr = f.readlines()

tline = len(xerr)
proc_line_nums = [i + 1 for i, l in enumerate(xerr) if 'proc' in l]
nline = proc_line_nums[-1] if proc_line_nums else 0
m = tline - nline

print(f'Number of Grade Exceed and Breaking linesearch at process [{last}]: {m}')
print("=" * 150)

if not graden:
    exit()

proc_dir = str(last).zfill(2)
print()
print(f"process: {proc_dir}")
print("-" * 70)
print(f"Line {graden[-1]} of {proc_dir}/B-relax/relax-output.txt file")
print("-" * 70)
subprocess.run(f'head -n {graden[0]} {proc_dir}/B-relax/relax-output.txt | tail -n 3', shell=True)
print()

if len(graden) == 2:
    print("-" * 70)
    print(f"Line {graden[-2]} of {proc_dir}/B-relax/relax-output.txt file")
    print("-" * 70)
    subprocess.run(f'head -n {graden[-1]} {proc_dir}/B-relax/relax-output.txt | tail -n 3', shell=True)
    print()

"""
#0 B-preselected.cfg
#1 C-selected.cfg
#2 D-computed.cfg
#3 E-train.cfg
#4 relaxed.cfg
#5 train.cfg
"""
