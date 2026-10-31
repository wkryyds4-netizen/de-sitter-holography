"""Archive versioned primary papers for local verification; no credentials required."""
from pathlib import Path
import hashlib
import json
import time
from datetime import datetime, timezone
import requests

ROOT = Path(__file__).resolve().parent.parent
SOURCES = {
    "giulini-marolf-gr-qc9902045v2": "https://arxiv.org/html/gr-qc/9902045v2",
    "held-maxfield-2509.05412v1": "https://arxiv.org/html/2509.05412v1",
    "han-hu-2610.14175v1": "https://arxiv.org/html/2610.14175v1",
    "cui-kolchmeyer-2607.13665v1": "https://arxiv.org/html/2607.13665v1",
    "harlow-zhao-2607.14215v2": "https://arxiv.org/html/2607.14215v2",
    "collier-et-al-2501.01486v1": "https://arxiv.org/html/2501.01486v1",
    "ali-et-al-2003.11046v2": "https://arxiv.org/html/2003.11046v2",
    "witten-2303.02837v2": "https://arxiv.org/html/2303.02837v2",
    "chen-stanford-tang-yang-2607.14042v1": "https://arxiv.org/html/2607.14042v1",
    "narovlansky-2506.02109v2": "https://arxiv.org/html/2506.02109v2",
    "pinamonti-0806.0803v2": "https://arxiv.org/pdf/0806.0803v2",
    "neeb-2011-aif2660": "https://www.numdam.org/item/10.5802/aif.2660.pdf",
    "kolchmeyer-liu-2411.08090v1": "https://arxiv.org/html/2411.08090v1",
    "chandra-liska-post-2610.30247v1": "https://arxiv.org/html/2610.30247v1",
    "bros-et-al-0901.4223v1": "https://arxiv.org/html/0901.4223v1",
    "bros-et-al-0901.4223v1.pdf": "https://arxiv.org/pdf/0901.4223v1",
    "dlmf-14.24": "https://dlmf.nist.gov/14.24",
    "marolf-yu-2610.24073v1": "https://arxiv.org/html/2610.24073v1",
    "goldberg-schlag-math0306108v2": "https://arxiv.org/html/math/0306108v2",
    "germain-pusateri-rousset-2018.pdf": "https://ems.press/content/serial-article-files/15932",
    "dlmf-15.8": "https://dlmf.nist.gov/15.8",
    "dlmf-15.2": "https://dlmf.nist.gov/15.2",
    "sun-2111.04591v2": "https://arxiv.org/html/2111.04591v2",
    "dlmf-14.12": "https://dlmf.nist.gov/14.12",
    "dlmf-14.20": "https://dlmf.nist.gov/14.20",
    "dlmf-14.3": "https://dlmf.nist.gov/14.3",
    "dlmf-5.11": "https://dlmf.nist.gov/5.11",
}

def main():
    folder = ROOT / "sources"
    folder.mkdir(exist_ok=True)
    manifest = []
    failures = []
    for name, url in SOURCES.items():
        is_pdf = name.endswith(".pdf") or "/pdf/" in url or url.endswith(".pdf")
        dest = folder / (name if name.endswith((".html", ".pdf")) else name + (".pdf" if is_pdf else ".html"))
        cached = dest.exists()
        if not cached:
            try:
                response = requests.get(url, timeout=45)
                response.raise_for_status()
                valid = response.content.lstrip().startswith(b"%PDF") if is_pdf else b"<html" in response.content.lower()
                if not valid:
                    raise ValueError(f"Unexpected file type: {url}")
                dest.write_bytes(response.content)
                time.sleep(1)
            except (requests.RequestException, ValueError) as exc:
                failures.append({"name":name,"url":url,"error":str(exc)})
                print(f"Unavailable: {name}: {exc}", flush=True)
                continue
        data = dest.read_bytes()
        manifest.append({"name": name, "url": url,
                         "file": dest.relative_to(ROOT).as_posix(),
                         "sha256": hashlib.sha256(data).hexdigest(),
                         "bytes": len(data), "cache_reused": cached})
        print(f"{name}: {len(data)} bytes", flush=True)
    (folder / "manifest.json").write_text(json.dumps({
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "sources": manifest, "unavailable": failures,
    }, indent=2) + "\n", encoding="utf-8")
    if failures:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
