#!/usr/bin/env python3
"""Run the three small drop fixtures and compare complete CPU ST/MT traces."""
import argparse
import csv
import json
from pathlib import Path
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('binary', type=Path)
p.add_argument('--cwd', type=Path, default=Path('.build/bin'))
p.add_argument('--output', type=Path, default=Path('.build/sized-validation'))
p.add_argument('--gif', action='store_true')
a = p.parse_args()
binary, cwd, output = a.binary.resolve(), a.cwd.resolve(), a.output.resolve()
reports = {}
for backend in ('cpu-st', 'cpu-mt'):
    for case in ('unit', 'double', 'mixed-face'):
        folder = output / backend / case
        folder.mkdir(parents=True, exist_ok=True)
        command = [str(binary), '--physics=' + backend, '--debug-scenario=sized-drop-' + case,
                   '--debug-steps=600', '--debug-output=' + str(folder)]
        if a.gif:
            command.append('--debug-gif-fps=20')
        with (folder / 'run.log').open('w') as log:
            result = subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        r = json.loads((folder / 'report.json').read_text())
        assert result.returncode == 0 and r['passed'], (backend, case, r['failures'])
        assert r['activeBackend'] == backend
        if backend == 'cpu-mt':
            assert r['sizedWorkerCalls'] > 0
        reports[backend, case] = r
        print(backend, case, r['sized'], flush=True)
for case in ('unit', 'double', 'mixed-face'):
    paths = [output / backend / case / 'sized-trace.csv' for backend in ('cpu-st', 'cpu-mt')]
    traces = [list(csv.DictReader(path.open())) for path in paths]
    assert len(traces[0]) == len(traces[1]) == 601
    for x, y in zip(*traces):
        for field in ('step', 'particles', 'independent', 'mass'):
            assert x[field] == y[field], (case, field, x['step'])
        for field in ('cx', 'cy', 'cz', 'penetration', 'strain', 'speed', 'drift', 'dependency'):
            assert abs(float(x[field]) - float(y[field])) <= 0.00005, (case, field, x['step'])
    particle_paths = [output / backend / case / 'sized-positions.csv' for backend in ('cpu-st', 'cpu-mt')]
    particle_traces = [list(csv.DictReader(path.open())) for path in particle_paths]
    assert len(particle_traces[0]) == len(particle_traces[1])
    for x, y in zip(*particle_traces):
        assert (x['step'], x['particle']) == (y['step'], y['particle'])
        for field in ('x', 'y', 'z', 'vx', 'vy', 'vz'):
            assert abs(float(x[field]) - float(y[field])) <= 0.00005, (case, field, x['step'], x['particle'])
(output / 'summary.json').write_text(json.dumps([
    dict(backend=backend, case=case, passed=r['passed'], workerCalls=r['sizedWorkerCalls'], **r['sized'])
    for (backend, case), r in reports.items()], indent=2) + '\n')
print('All six drops and full CPU trajectory comparisons passed.')
