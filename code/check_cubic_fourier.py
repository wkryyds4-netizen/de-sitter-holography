"""High-precision independent checks of the limiting cubic frequency density."""
from pathlib import Path
import json,hashlib
import mpmath as mp
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-kernel'
def density(w):
    if not w:return mp.mpf(8)/15
    return mp.pi*w*(w*w+16)*(w*w+64)/(7680*mp.sinh(mp.pi*w/4))
def run(dps):
    mp.mp.dps=dps;rows=[]
    for w in [0,1,4,12]:
        actual=2*mp.quad(lambda t:mp.sech(2*t)**6*mp.cos(w*t),[0,1,2,4,8,mp.inf])
        expected=density(mp.mpf(w))
        rows.append(dict(kind='Fourier',w=w,value=mp.nstr(actual,dps),expected=mp.nstr(expected,dps),
          relative_error=mp.nstr(abs(actual/expected-1),dps)))
    for n,expected in [(0,1),(2,24),(4,1920)]:
        actual=mp.quad(lambda w:w**n*density(w)/mp.pi,[0,1,4,12,32,mp.inf])
        rows.append(dict(kind='moment',n=n,value=mp.nstr(actual,dps),expected=str(expected),
          relative_error=mp.nstr(abs(actual/expected-1),dps)))
    assert max(mp.mpf(r['relative_error']) for r in rows)<mp.mpf('1e-40')
    return rows
def main():
    target=OUT/'fourier-check.json'
    if target.exists():raise FileExistsError(target)
    low=run(50);high=run(75)
    max_change=max(abs(mp.mpf(a['value'])/mp.mpf(b['value'])-1) for a,b in zip(low,high))
    result=dict(status='passed',Js=1,precisions=[50,75],checks_per_precision=7,
      max_relative_precision_change=mp.nstr(max_change,20),low_precision=low,high_precision=high,
      scope='numerical quadrature check of analytic beta/Gamma transform; not interval certification',
      environment=dict(mpmath=mp.__version__),
      input_sha256={'code/check_cubic_fourier.py':hashlib.sha256((ROOT/'code/check_cubic_fourier.py').read_bytes()).hexdigest()})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['max_relative_precision_change'])
if __name__=='__main__':main()
