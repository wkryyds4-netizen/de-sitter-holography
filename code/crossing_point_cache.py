"""Content-addressed raw 6j point cache; only the parent process writes SQLite.

The point identity includes exact decimal input strings, precision settings,
kernel/cache code hashes and the archived primary source. No Gaussian, time,
or quadrature weight is part of a raw point. Frozen studies export every used
record separately so later cache growth cannot make old result hashes stale.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,sqlite3,platform
import mpmath
from observer_crossing_regular import RegularCrossingEvaluator
ROOT=Path(__file__).resolve().parents[1]
SCHEMA="crossing-point-v1"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False)


def digest_text(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def configuration(dps=15,guard_digits=8):
    return dict(schema=SCHEMA,dps=int(dps),guard_digits=int(guard_digits),
        max_working_dps=180,sum_method="r+s+e",mpmath=mpmath.__version__,
        python=platform.python_version(),
        code_sha256={n:sha(ROOT/"code"/n) for n in [
            "crossing_point_cache.py","observer_crossing_regular.py","observer_crossing.py"]},
        source_sha256=sha(ROOT/"sources/kolchmeyer-liu-2411.08090v1.html"))


def point_identity(config,nu,labels):
    if len(labels)!=4:raise ValueError("Four ordered labels required")
    return dict(configuration=config,nu=str(nu),labels=[str(x) for x in labels])


def point_key(identity):
    return digest_text(canonical(identity))


def point_evaluate(identity):
    cfg=identity["configuration"]
    ev=RegularCrossingEvaluator(cfg["dps"],cfg["guard_digits"],
        cfg["max_working_dps"],cfg["sum_method"])
    labels=identity["labels"]
    r=ev.evaluate(identity["nu"],*labels)
    c=ev.ctx
    with c.workdps(r.working_dps):
        cond=r.hh_cross/ev.rho_hh(labels[-1])
        def real_text(x):return c.nstr(x,n=r.working_dps)
        def complex_text(x):return dict(real=real_text(c.re(x)),imag=real_text(c.im(x)))
        output=dict(identity=identity,key=point_key(identity),
            conditional=complex_text(cond),hh_cross=complex_text(r.hh_cross),
            symmetric_s=complex_text(r.symmetric_s),
            method=r.method,chart_labels=[real_text(x) for x in r.chart_labels],
            chart_rotation=r.chart_rotation,requested_dps=r.requested_dps,
            working_dps=r.working_dps,max_observed_formula_dps=r.max_observed_formula_dps,
            cancellation_ratio=real_text(r.cancellation_ratio),
            reality_relative_error=real_text(r.reality_relative_error),
            residue_relative_error=real_text(r.residue_relative_error),
            refinement_count=r.refinement_count,hyper_calls=r.hyper_calls,
            elapsed_seconds=r.elapsed_seconds,
            generated_utc=datetime.now(timezone.utc).isoformat())
    return output


class PointCache:
    def __init__(self,path):
        self.path=Path(path)
        self.path.parent.mkdir(parents=True,exist_ok=True)
        self.db=sqlite3.connect(self.path)
        self.db.execute("CREATE TABLE IF NOT EXISTS points (key TEXT PRIMARY KEY, payload TEXT NOT NULL, payload_sha256 TEXT NOT NULL)")
        self.db.commit()
        self.pending=0

    def get(self,identity):
        key=point_key(identity)
        row=self.db.execute("SELECT payload,payload_sha256 FROM points WHERE key=?",(key,)).fetchone()
        if row is None:return None
        payload,expected=row
        if digest_text(payload)!=expected:raise ValueError("Corrupt point-cache payload: "+key)
        record=json.loads(payload)
        if record.get("identity")!=identity or record.get("key")!=key:
            raise ValueError("Point-cache identity mismatch: "+key)
        return record

    def put(self,record):
        key=record["key"]
        if point_key(record["identity"])!=key:raise ValueError("Refuse mismatched point key")
        payload=canonical(record);checksum=digest_text(payload)
        old=self.db.execute("SELECT payload_sha256 FROM points WHERE key=?",(key,)).fetchone()
        if old is not None:
            # No replacement even when a recomputation only changes timing.
            existing=self.get(record["identity"])
            if existing!=record:raise ValueError("Refuse to overwrite an existing point")
            return
        self.db.execute("INSERT INTO points VALUES (?,?,?)",(key,payload,checksum))
        self.pending+=1
        if self.pending>=32:self.commit()

    def commit(self):
        self.db.commit();self.pending=0

    def close(self):
        self.commit();self.db.close()

    def __enter__(self):return self
    def __exit__(self,*args):self.close()