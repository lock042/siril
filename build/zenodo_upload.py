#!/usr/bin/env python3
"""Publish a Siril source tarball as a new version of the Zenodo record.

Usage: zenodo_upload.py <tarball> <version> [--date YYYY-MM-DD] [--draft]
Needs ZENODO_TOKEN (scopes: deposit:write, deposit:actions). Set
ZENODO_URL=https://sandbox.zenodo.org to test (with a sandbox token and
ZENODO_RECORD pointing to a sandbox record).
"""

import argparse
import datetime
import json
import os
import sys
import urllib.error
import urllib.request

import yaml

ZENODO_URL = os.environ.get("ZENODO_URL", "https://zenodo.org").rstrip("/")
# any published version of the Siril concept record (10.5281/zenodo.13880497)
ZENODO_RECORD = os.environ.get("ZENODO_RECORD", "13880498")
CFF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "CITATION.cff")


def request(method, url, data=None, content_type="application/json", auth=True):
    headers = {"Content-Type": content_type}
    if auth:
        headers["Authorization"] = "Bearer " + os.environ["ZENODO_TOKEN"]
    if data is not None and content_type == "application/json":
        data = json.dumps(data).encode()
    req = urllib.request.Request(url, data=data, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req) as resp:
            body = resp.read()
    except urllib.error.HTTPError as e:
        sys.exit(f"{method} {url} failed: {e.code} {e.read().decode(errors='replace')}")
    return json.loads(body) if body else None


def metadata(version, date):
    with open(CFF, encoding="utf-8") as f:
        cff = yaml.safe_load(f)
    creators = []
    for a in cff["authors"]:
        c = {"name": f"{a['family-names']}, {a['given-names']}"}
        if "orcid" in a:
            c["orcid"] = a["orcid"].rsplit("/", 1)[-1]
        if "affiliation" in a:
            c["affiliation"] = a["affiliation"]
        creators.append(c)
    tag_url = f"{cff['repository-code']}/-/tree/{version}"
    article = cff["preferred-citation"]
    return {
        "upload_type": "software",
        "title": cff["title"],
        "version": version,
        "publication_date": date,
        "creators": creators,
        "description": f"<p>{cff['abstract'].strip()}</p>"
                       f"<p>Source code of Siril {version}: "
                       f"<a href=\"{tag_url}\">{tag_url}</a></p>",
        "license": cff["license"].lower(),
        "keywords": cff["keywords"],
        "related_identifiers": [
            {"identifier": article["doi"], "relation": "isDescribedBy",
             "scheme": "doi", "resource_type": "publication-article"},
            {"identifier": tag_url, "relation": "isSupplementTo",
             "scheme": "url", "resource_type": "software"},
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tarball")
    parser.add_argument("version")
    parser.add_argument("--date", default=datetime.date.today().isoformat())
    parser.add_argument("--draft", action="store_true", help="do not publish")
    args = parser.parse_args()

    meta = metadata(args.version, args.date)
    api = f"{ZENODO_URL}/api/deposit/depositions"

    latest = request("GET", f"{ZENODO_URL}/api/records/{ZENODO_RECORD}/versions/latest",
                     auth=False)
    if latest.get("metadata", {}).get("version") == args.version:
        sys.exit(f"version {args.version} is already the latest on Zenodo")

    # returns the pending draft instead if a previous run did not publish
    res = request("POST", f"{api}/{latest['id']}/actions/newversion")
    draft_url = res.get("links", {}).get("latest_draft", f"{api}/{res['id']}")
    draft = request("GET", draft_url)

    for f in draft.get("files", []):
        request("DELETE", f"{api}/{draft['id']}/files/{f['id']}")
    with open(args.tarball, "rb") as f:
        request("PUT", f"{draft['links']['bucket']}/{os.path.basename(args.tarball)}",
                data=f.read(), content_type="application/octet-stream")
    request("PUT", f"{api}/{draft['id']}", data={"metadata": meta})

    if args.draft:
        print(f"draft ready: {draft['links']['html']}")
        return
    pub = request("POST", f"{api}/{draft['id']}/actions/publish")
    print(f"published {args.version}: https://doi.org/{pub['doi']}")


if __name__ == "__main__":
    main()
