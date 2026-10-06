#!/usr/bin/env python


"""Copyright (c) 2005-2026, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# Flag PRs merged into develop since the last release that are not yet cited on
# the "changes since last release" page. Findings go in the AUDIT-TABLE region of
# tracking issue #641; --dry-run prints them instead.
#
# Workflow secrets: GH_TOKEN reads Chaste/Chaste and edits #641; WEBSITE_TOKEN
# reads the private Chaste.github.io page (local runs fall back to gh's own auth);
# ORG_READ_TOKEN (needs read:org) lists org members - without it, every human
# participant is tagged rather than just members.

import argparse
import json
import os
import re
import subprocess
import sys
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone

REPO = "Chaste/Chaste"
ORG, NAME = REPO.split("/")
BRANCH = "develop"
WEBSITE_REPO = "Chaste/Chaste.github.io"
CHANGES_PATH = "site/content/docs/release-notes/changes-since-last-release.md"
TRACKING_ISSUE = 641
EXEMPT_LABEL = "not for release notes"

# @-mentioned when a PR has no org-member owner.
BACKUP_TAGGEES = ["fcooper8472", "mirams"]
# Never @-mentioned: service accounts, and bots GitHub reports as ordinary users.
NON_HUMAN = {"ChastePi", "chastedevelopment", "Copilot", "copilot-swe-agent"}
# Fewer members than this coming back means the org fetch failed.
MIN_MEMBERS = 15

TABLE_START, TABLE_END = "<!-- AUDIT-TABLE:START -->", "<!-- AUDIT-TABLE:END -->"


def run(cmd, stdin=None, env=None):
    """Run a command, exit on failure, return stdout."""
    p = subprocess.run(cmd, capture_output=True, text=True, input=stdin,
                       env={**os.environ, **env} if env else None)
    if p.returncode:
        sys.exit(f"command failed: {' '.join(cmd)}\n{p.stderr.strip()}")
    return p.stdout


def gh_json(cmd):
    return json.loads(run(cmd))


def issue(subcommand, *opts, **kw):
    """Run `gh issue <subcommand> #641 ...`."""
    return run(["gh", "issue", subcommand, str(TRACKING_ISSUE), "--repo", REPO, *opts], **kw)


# --- gather ---

def last_release():
    r = gh_json(["gh", "api", f"repos/{REPO}/releases/latest"])
    return r["name"] or r["tag_name"], r["published_at"]


def merged_prs(since):
    prs = gh_json(["gh", "pr", "list", "--repo", REPO, "--state", "merged",
                   "--base", BRANCH, "--limit", "1000",
                   "--search", f"merged:>={since[:10]}",
                   "--json", "number,title,mergedAt,author,mergedBy,labels"])
    return [p for p in prs if p["mergedAt"] > since]


def org_members():
    """Org member logins, or None if unavailable (ORG_READ_TOKEN missing/unscoped)."""
    token = os.environ.get("ORG_READ_TOKEN")
    p = subprocess.run(["gh", "api", f"orgs/{ORG}/members", "--paginate", "--jq", ".[].login"],
                       capture_output=True, text=True,
                       env={**os.environ, "GH_TOKEN": token} if token else None)
    members = set(p.stdout.split()) if not p.returncode else set()
    if len(members) < MIN_MEMBERS:
        why = (p.stderr.strip().splitlines() or [f"{len(members)} returned"])[-1]
        print(f"warning: {ORG} member list unavailable ({why}); tagging all human "
              f"participants - does ORG_READ_TOKEN carry the read:org scope?")
        return None
    return members


def pr_details(numbers):
    """number -> (closing issue numbers, participant logins)."""
    out = {}
    for i in range(0, len(numbers), 40):
        batch = numbers[i:i + 40]
        fields = " ".join(
            f"p{n}: pullRequest(number: {n}) {{ "
            f"closingIssuesReferences(first: 20) {{ nodes {{ number }} }} "
            f"participants(first: 100) {{ nodes {{ login }} }} }}"
            for n in batch)
        query = f'query {{ repository(owner: "{ORG}", name: "{NAME}") {{ {fields} }} }}'
        repo = gh_json(["gh", "api", "graphql", "-f", f"query={query}"])["data"]["repository"]
        for n in batch:
            node = repo[f"p{n}"]
            out[n] = ({x["number"] for x in node["closingIssuesReferences"]["nodes"]},
                      {x["login"] for x in node["participants"]["nodes"]})
    return out


def changes_page():
    path = f"repos/{WEBSITE_REPO}/contents/{CHANGES_PATH}"
    token = os.environ.get("WEBSITE_TOKEN")
    if not token:  # local runs: use gh's own credentials
        return run(["gh", "api", path, "-H", "Accept: application/vnd.github.raw"])
    req = urllib.request.Request(f"https://api.github.com/{path}", headers={
        "Authorization": f"token {token}",
        "Accept": "application/vnd.github.raw",
        "User-Agent": "chaste-release-notes-audit"})
    with urllib.request.urlopen(req) as resp:
        return resp.read().decode()


# --- matching ---

def title_issues(title):
    """Issue numbers in a PR title: an explicit '#NNN', or a leading '[#]NNN'
    followed by space/colon (so a version like '3.25 support' is not a match)."""
    nums = {int(n) for n in re.findall(r"#(\d+)\b", title)}
    lead = re.match(r"\s*#?(\d+)(?=[\s:]|$)", title)
    return nums | ({int(lead.group(1))} if lead else set())


def referenced(page):
    """Issue/PR numbers cited on the changes page: '#NNN' tokens and Chaste
    issue/pull URLs only, so bare years and version numbers are ignored."""
    return {int(n) for n in re.findall(r"Chaste/Chaste/(?:issues|pull)/(\d+)", page)} \
        | {int(n) for n in re.findall(r"(?<![\w/])#(\d+)\b", page)}


def recorded_reason(pr, refs):
    """A short string saying why the PR counts as recorded, or None."""
    if any(label["name"] == EXEMPT_LABEL for label in pr["labels"]):
        return "labelled"
    if pr["number"] in refs:
        return f"PR #{pr['number']} cited"
    hit = (pr["closes"] | pr["title_issues"]) & refs
    return f"issue #{min(hit)} cited" if hit else None


def is_bot(user):
    login = (user or {}).get("login", "")
    return (not user or bool(user.get("is_bot")) or login in NON_HUMAN
            or login.endswith("[bot]") or login.startswith("app/"))


def taggees(pr, members):
    """Who to @-mention: the PR author, else the merger, else org-member
    participants, else the backup list. members is None when the member list is
    unavailable, in which case everyone counts as in-org."""
    for user in (pr["author"], pr["mergedBy"]):
        if not is_bot(user) and (members is None or user["login"] in members):
            return [user["login"]]
    humans = {p for p in pr["participants"]
              if p not in NON_HUMAN and not p.endswith("[bot]")}
    return sorted(humans if members is None else humans & members) or list(BACKUP_TAGGEES)


# --- report + reconcile ---

def _row(pr, members):
    title = pr["title"].replace("|", "\\|")
    tags = " ".join("@" + t for t in taggees(pr, members))
    return f"| #{pr['number']} | {title} | {tags} |"


def section_md(release, since, errant, checked, members):
    """The Markdown for the AUDIT-TABLE region of issue #641."""
    head = f"Merged into `{BRANCH}` since **{release}** ({since[:10]}). _Checked {checked}._"
    if not errant:
        return f"{head}\n\n✅ Nothing outstanding."
    rows = "\n".join(_row(pr, members) for pr in errant)
    return f"{head}\n\n| PR | Title | Owner |\n|----|-------|-------|\n{rows}"


def prs_in_body(body):
    return {int(n) for n in re.findall(r"\|\s*#(\d+)\s*\|", body or "")}


def splice(body, section):
    """Replace the AUDIT-TABLE region, appending the markers if they are absent."""
    block = f"{TABLE_START}\n{section}\n{TABLE_END}"
    pattern = re.compile(re.escape(TABLE_START) + ".*?" + re.escape(TABLE_END), re.DOTALL)
    if pattern.search(body):
        return pattern.sub(lambda _: block, body)
    print("warning: AUDIT-TABLE markers missing from issue #641; appending them")
    return f"{body.rstrip()}\n\n{block}\n"


def reconcile(section, errant, checked):
    before = json.loads(issue("view", "--json", "state,body"))
    was_open = before["state"] == "OPEN"
    issue("edit", "--body-file", "-", stdin=splice(before["body"], section))

    if not errant:
        if was_open:
            issue("comment", "--body", f"✅ All merged PRs recorded as of {checked}. Closing.")
            issue("close")
        return

    if not was_open:
        issue("reopen")
    new = sorted({pr["number"] for pr in errant} - prs_in_body(before["body"]))
    if new and was_open:
        issue("comment", "--body",
              f"Newly unrecorded as of {checked}: {', '.join(f'#{n}' for n in new)}")


def main():
    parser = argparse.ArgumentParser(
        description="Audit merged PRs against the changes-since-last-release page.")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the findings instead of editing issue #641")
    dry_run = parser.parse_args().dry_run

    checked = datetime.now(timezone.utc).date().isoformat()
    members = org_members()
    release, since = last_release()
    prs = merged_prs(since)
    print(f"{len(prs)} PR(s) merged into {BRANCH} since {release} ({since[:10]})")

    details = pr_details([pr["number"] for pr in prs])
    for pr in prs:
        pr["closes"], pr["participants"] = details[pr["number"]]
        pr["title_issues"] = title_issues(pr["title"])
    refs = referenced(changes_page())

    errant, by_issue = [], defaultdict(list)
    for pr in sorted(prs, key=lambda p: p["number"]):
        reason = recorded_reason(pr, refs)
        if not reason:
            errant.append(pr)
        tags = " ".join("@" + t for t in taggees(pr, members))
        print(f"  #{pr['number']}: {reason or 'NOT RECORDED -> ' + tags}")
        for n in (pr["closes"] | pr["title_issues"]) & refs:
            by_issue[n].append(pr["number"])

    for n, nums in sorted(by_issue.items()):
        if len(nums) > 1:
            print(f"note: changes-page #{n} stands in for PRs "
                  f"{', '.join(f'#{x}' for x in nums)}")

    section = section_md(release, since, errant, checked, members)
    if dry_run:
        print(f"\n--- issue #{TRACKING_ISSUE} section ---\n\n{section}")
    else:
        reconcile(section, errant, checked)


if __name__ == "__main__":
    main()
