### 📌 High-Performance Batched PubChem Similarity SDF Downloader[cite: 9]

[![View Script](https://img.shields.io/badge/VIEW_SCRIPT-PUBCHEM__SIMILARITY__DOWNLOADER.PY-blue?style=for-the-badge&logo=python)](#)
[![Status](https://img.shields.io/badge/STATUS-EXECUTABLE-brightgreen?style=for-the-badge)](#)

<details>
<summary>🔍 <b>Click to View Script Details & Usage Header</b></summary>

```python
# ==============================================================================
# SCRIPT NAME : pubchem_similarity_downloader.py
# AUTHOR      : K. R. Janahan
# DESCRIPTION : Fast, batched downloader that resolves PubChem CIDs for compound
#               names, fetches all 2D similar structures above a given percentage 
#               threshold, and retrieves 3D/2D SDF records using bulk requests.
# DEPENDENCIES: Requests, Python3
# USAGE       : python3 pubchem_similarity_downloader.py
# ==============================================================================
```
</details>

<details>
<summary>💻 <b>Click to View / Copy Full Source Code</b></summary>

```python

"""
PubChem Similarity SDF Downloader (fast, batched)
===================================================
Add any compound name + desired similarity % below and run the script.

For every entry in COMPOUNDS this script will:
  1. Resolve the PubChem CID from the compound name
  2. Find ALL compounds at or above the given 2D similarity threshold
     (no cap on result count)
  3. Download SDFs for every match — in BATCHES of CHUNK_SIZE CIDs per
     request (3D conformer if available, otherwise falls back to 2D),
     instead of one request per compound. This is what makes it fast:
     PubChem itself recommends sending CID lists in bulk rather than
     making one request per identifier.

Files land in "<CompoundName>_<threshold>_similarity/".
Already-downloaded files are skipped, so you can stop and re-run safely.

Usage:
    1. Edit the COMPOUNDS dict below.
    2. python pubchem_similarity_downloader.py
"""

import os
import re
import time
import requests
import sys

# --- JANAHAN OS // UNIVERSAL STARTUP ANIMATION BLOCK ---
try:
    import sys, os
    # 1. Automatically inject your custom tool folder into the path
    custom_path = os.path.expanduser("~/intro_script")
    if custom_path not in sys.path:
        sys.path.insert(0, custom_path)

    # 2. Dynamically detect the name of the script running
    try:
        script_name = os.path.splitext(os.path.basename(__file__))[0].replace("_", " ").upper()
    except NameError:
        script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0].replace("_", " ").upper() if sys.argv and sys.argv[0] else "AUTOMATION PIPELINE"

    # 3. Trigger animation with dynamic naming
    from bio_intro import play_intro
    play_intro(f"PROTOCOL: {script_name} // EXECUTION MATRIX")

except Exception as e:
    # DEBUG FLAG: Print the exact reason why it fails instead of silently passing!
    print(f"\n[JANAHAN OS NOTICE] Startup animation skipped: {e}\n")
# -------------------------------------------------------

# =========================================================
# CONFIG — just add "Compound Name": similarity_percentage
# =========================================================
COMPOUNDS = {
    "Imipenem": 90,
    "Meropenem": 95,
    "Ertapenem": 90,
    "Doripenem": 95,
    "Panipenem": 90,
    "Tebipenem": 95,
}

# How many CIDs to send per batch request. PubChem handles a few hundred
# fine via POST (no URL length limit since CIDs go in the request body).
CHUNK_SIZE = 200

# Small pause between BATCH requests (not per-file) to stay well within
# PubChem's request-rate policy.
REQUEST_DELAY = 0.2

# Retries per batch request
MAX_RETRIES = 5

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

session = requests.Session()
session.headers.update({"User-Agent": "Mozilla/5.0"})


def chunked(items, size):
    for i in range(0, len(items), size):
        yield items[i:i + size]


def get_cid_from_name(name):
    """Resolve a PubChem CID from a compound name."""
    url = f"{PUG}/compound/name/{requests.utils.quote(name)}/cids/JSON"
    resp = session.get(url, timeout=30)
    resp.raise_for_status()
    return resp.json()["IdentifierList"]["CID"][0]


def get_similar_cids(cid, threshold):
    """Return every CID at or above the given 2D similarity threshold."""
    url = (
        f"{PUG}/compound/fastsimilarity_2d/cid/{cid}/cids/JSON"
        f"?Threshold={threshold}"
    )
    resp = session.get(url, timeout=60)
    resp.raise_for_status()
    return resp.json()["IdentifierList"]["CID"]


def fetch_titles(cids):
    """Batch-fetch {cid: title} for a list of CIDs."""
    titles = {}
    for chunk in chunked(cids, CHUNK_SIZE):
        for attempt in range(MAX_RETRIES):
            try:
                resp = session.post(
                    f"{PUG}/compound/cid/property/Title/JSON",
                    data={"cid": ",".join(map(str, chunk))},
                    timeout=60,
                )
                if resp.status_code == 200:
                    for prop in resp.json()["PropertyTable"]["Properties"]:
                        titles[prop["CID"]] = prop.get("Title") or f"CID{prop['CID']}"
                break
            except (requests.RequestException, KeyError, ValueError):
                time.sleep(1)
        time.sleep(REQUEST_DELAY)
    return titles


def parse_sdf_records(text):
    """Split a multi-record SDF blob into {cid: sdf_text}."""
    records = {}
    for block in text.split("$$$$\n"):
        if not block.strip():
            continue
        m = re.search(r"<PUBCHEM_COMPOUND_CID>\s*\n(\d+)", block)
        if m:
            records[int(m.group(1))] = block.rstrip("\n") + "\n$$$$\n"
    return records


def fetch_sdf_records(cids, record_type):
    """Batch-fetch {cid: sdf_text} for a list of CIDs at the given record type."""
    records = {}
    for chunk in chunked(cids, CHUNK_SIZE):
        for attempt in range(MAX_RETRIES):
            try:
                resp = session.post(
                    f"{PUG}/compound/cid/SDF",
                    data={"cid": ",".join(map(str, chunk)), "record_type": record_type},
                    timeout=120,
                )
                if resp.status_code == 200 and resp.text.strip():
                    records.update(parse_sdf_records(resp.text))
                break
            except requests.RequestException:
                time.sleep(1)
        time.sleep(REQUEST_DELAY)
    return records


def clean_filename(name):
    for ch in '/\\:*?"<>|':
        name = name.replace(ch, "_")
    return name.replace(" ", "_")


def process_compound(compound_name, threshold):
    print(f"\n========== {compound_name} ({threshold}% similarity) ==========")

    try:
        cid = get_cid_from_name(compound_name)
        print(f"Resolved {compound_name} -> CID {cid}")
    except Exception as e:
        print(f"Could not resolve CID for {compound_name}: {e}")
        return

    folder_name = f"{clean_filename(compound_name)}_{threshold}_similarity"
    os.makedirs(folder_name, exist_ok=True)

    try:
        similar_cids = get_similar_cids(cid, threshold)
        print(f"Found {len(similar_cids)} compounds at >= {threshold}% similarity")
    except Exception as e:
        print(f"Error retrieving similar compounds for {compound_name}: {e}")
        return

    if not similar_cids:
        return

    print(f"Fetching titles in batches of {CHUNK_SIZE}...")
    titles = fetch_titles(similar_cids)

    print(f"Fetching 3D SDFs in batches of {CHUNK_SIZE}...")
    sdf_3d = fetch_sdf_records(similar_cids, "3d")

    missing = [c for c in similar_cids if c not in sdf_3d]
    sdf_2d = {}
    if missing:
        print(f"{len(missing)} compounds had no 3D conformer, fetching 2D for those...")
        sdf_2d = fetch_sdf_records(missing, "2d")

    total = len(similar_cids)
    downloaded, skipped, failed = 0, 0, 0
    for i, cid_ in enumerate(similar_cids, start=1):
        title = clean_filename(titles.get(cid_, f"CID{cid_}"))
        filename = f"{title}_CID{cid_}.sdf"
        filepath = os.path.join(folder_name, filename)

        if os.path.exists(filepath):
            skipped += 1
            continue

        block = sdf_3d.get(cid_) or sdf_2d.get(cid_)
        if block:
            with open(filepath, "w", encoding="utf-8") as f:
                f.write(block)
            downloaded += 1
            record_type = "3D" if cid_ in sdf_3d else "2D"
            if i % 25 == 0 or i == total:
                print(f"[{i}/{total}] ... {downloaded} saved so far ({record_type} last)")
        else:
            failed += 1

    print(
        f"{compound_name}: {downloaded} downloaded, {skipped} already had, "
        f"{failed} failed (no SDF available)"
    )


def main():
    for compound_name, threshold in COMPOUNDS.items():
        process_compound(compound_name, threshold)
    print("\nFinished downloading all similarity sets.")


if __name__ == "__main__":
    main()

```
</details>
