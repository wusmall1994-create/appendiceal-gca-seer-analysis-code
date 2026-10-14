"""Audit existing official exports; does not execute Joinpoint.
Usage: python Code/audit_fixed_rho_official.py Source_Data/Joinpoint_Error_Sensitivity
"""
import csv, hashlib, json, re, sys
from pathlib import Path

root = Path(sys.argv[1])
def read(path, delimiter=','):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))
mapping = {(r['batch'], r['cohort_id']): r for r in read(root/'Cohort_Map.csv')}
results, audits, hashes = [], [], {}
for plan in read(root/'Run_Plan.csv'):
    prefix = plan['prefix']; batch = prefix.split('_rho')[0]; rho = float(plan['fixed_rho'])
    def file(suffix):
        p = root/'Exports'/(prefix+suffix)
        assert p.is_file(), p
        hashes[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        return p
    for suffix in ['.jpo','.Export.APC.txt','.Export.Model.Selection.txt']:
        file(suffix)
    summary = file('.Export.Run.Summary.txt').read_text(encoding='utf-8-sig')
    assert summary.count('Joinpoint Regression Program Version') == 1
    for pattern in [rf'correlation = {rho:.4f}',r'Maximum\s+4\b',r'Minimum\s+0\b',
        r'Number of Permutations for Each Test\s+4499',r'Seed for Randomly Permuting Data\s+7160',
        r'Number of Resamples\s+5001',r'Empirical Quantile Final Selected Model Only',
        r'Number of Cohorts That Couldn.t Be Executed\s+0',r'Number of Excluded Cohorts\s+0',
        r'Log Transformation\s+Yes',r'Selected Option\s+Standard Error',r'Constant Seed',
        r'Min. Obs. From Joinpoint to End of Data\s+2',r'Min. Obs. Between Two Joinpoints\s+2',
        r'Points Between Adjacent Observed X Values\s+0',r'Permutation Test',
        rf'Full Range\s+{plan["AAPC_start"]} - {plan["AAPC_end"]}']:
        assert re.search(pattern, summary), (prefix, pattern)
    assert plan['input_file'].lower() in summary.lower()
    inputs = read(root/plan['input_file'], '\t')
    lookup = {(r['cohort_id'],r['year']):r for r in inputs}
    assert len(lookup) == len(inputs) == int(plan['input_rows'])
    data = read(file('.Export.Data.txt'), '\t')
    assert len(data) == 5*len(inputs)
    assert len({(r['cohort_id'],r['year'],r['Joinpoints']) for r in data}) == len(data)
    for r in data:
        source = lookup[r['cohort_id'],r['year']]
        assert abs(float(source['rate'])-float(r['Observed rate'])) <= .005000001
        assert abs(float(source['se'])-float(r['Standard Error'])) <= .000050001
    selected = read(file('.Export.Final.Selected.Model.txt'), '\t')
    assert len(selected) == int(plan['cohorts'])
    assert {r['cohort_id'] for r in selected} == {r['cohort_id'] for r in inputs}
    aapcs = read(file('.Export.AAPC.txt'), '\t')
    estimates = read(file('.Export.Model.Estimates.txt'), '\t')
    for r in estimates:
        assert abs(float(r['Auto Corr'])-rho) < 1e-8
    for selected_row in selected:
        cohort, model = selected_row['cohort_id'], selected_row['Model']
        matches = [r for r in aapcs if r['cohort_id']==cohort and r['Joinpoint Model']==model]
        assert len(matches)==1
        a = matches[0]
        assert a['Start Obs']==plan['AAPC_start'] and a['End Obs']==plan['AAPC_end']
        knots = sorted({r['Joinpoint'] for r in estimates if r['cohort_id']==cohort and r['Model']==model and r['Joinpoint']!='NA'})
        assert len(knots)==int(model)
        results.append(dict(batch=batch,rho=rho,label=mapping[batch,cohort]['label'],model=model,knots=';'.join(knots),AAPC=a['AAPC'],lower=a['AAPC C.I. Low'],upper=a['AAPC C.I. High']))
    audits.append(dict(prefix=prefix,rho=rho,cohorts=len(selected),inputs=len(inputs),exported_rows=len(data),status='PASS'))
assert len(results)==48 and len(audits)==9
with (root/'Fixed_Rho_Official_Summary.csv').open('w',encoding='utf-8',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
report=dict(status='PASS',sessions=audits,completed_sessions=9,selected_models=48,
    input_observations=sum(r['inputs'] for r in audits),exported_observations=sum(r['exported_rows'] for r in audits),
    all_selected_AAPC_lower_bounds_positive=all(float(r['lower'])>0 for r in results),
    note='Fixed positive correlations are sensitivity assumptions, not estimated dependence. Archive_actual_rho021 is excluded.',sha256=hashes)
(root/'Fixed_Rho_Official_Audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('PASS: 9 sessions, 48 selected fits, 1383 input and 6915 model-specific rows.')
