#!/usr/bin/env python3
"""Verify the immutable handoff metadata; never runs Julia or modifies files."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote


def digest(data):
    return hashlib.sha256(data).hexdigest()


def require(condition, message):
    if not condition:
        raise SystemExit('HANDOFF CHECK FAILED: ' + message)


def git(root, *args):
    return subprocess.run(['git', '-C', str(root), *args], check=True,
                          stdout=subprocess.PIPE).stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--allow-source-changes', action='store_true',
                        help='Permit later source edits while verifying the frozen packet')
    options = parser.parse_args()
    packet = Path(__file__).resolve().parent
    root = packet.parents[1]
    manifest = json.loads((packet / 'MANIFEST.json').read_text())
    for name, expected in manifest['artifacts'].items():
        path = packet / name
        require(path.is_file(), 'missing packet file: ' + name)
        data = path.read_bytes()
        require(len(data) == expected['bytes'] and digest(data) == expected['sha256'],
                'packet hash/length: ' + name)
    actual = {str(p.relative_to(packet)) for p in packet.rglob('*') if p.is_file()}
    require(actual == set(manifest['artifacts']) | {'MANIFEST.json'},
            'packet inventory contains missing or unlisted files')
    for name, expected in manifest['documents'].items():
        require(digest((root / name).read_bytes()) == expected,
                'handoff document changed: ' + name)
    baseline = manifest['code_commit']
    paths = ['src', 'ext', 'test', 'Project.toml', '.github/workflows']
    frozen = (packet / 'CODE_TREE_BASELINE.txt').read_bytes()
    require(digest(frozen) == manifest['code_tree_listing_sha256'], 'code listing hash')
    require(git(root, 'ls-tree', '-r', baseline, '--', *paths) == frozen,
            'baseline Git objects differ from packet')
    current_matches = git(root, 'ls-tree', '-r', 'HEAD', '--', *paths) == frozen
    # Include uncommitted and untracked source in the comparison.
    source_status = git(root, 'status', '--porcelain', '--untracked-files=all', '--', *paths)
    if not options.allow_source_changes:
        require(current_matches and not source_status, 'source changed from frozen baseline')
    crosswalk = json.loads((packet / 'ROADMAP_CROSSWALK.json').read_text())
    rows = crosswalk['items']
    require(len(rows) == crosswalk['count'] == 79, '79-row crosswalk')
    require(sum(r['legacy_checked'] for r in rows) == crosswalk['legacy_checked'] == 23,
            'legacy checked count')
    require(crosswalk['open'] == 56, 'open count')
    require([r['id'] for r in rows] == [f'M{i:03d}' for i in range(1, 80)],
            'milestone ID sequence')
    source = git(root, 'show', baseline + ':' + crosswalk['source_file']).decode().splitlines()
    require(sum(bool(re.match(r'^\s*- \[[ xX]\]', line)) for line in source) == 79,
            'original milestone count')
    for row in rows:
        index = row['source_line'] - 1
        match = re.match(r'^\s*- \[([ xX])\]\s*(.*)$', source[index])
        require(match is not None, 'source milestone location: ' + row['id'])
        parts = [match[2]]
        index += 1
        while index < len(source) and source[index].startswith('  ') and source[index].strip():
            parts.append(source[index].strip())
            index += 1
        text = ' '.join(parts)
        require(text == row['source_text'] and digest(text.encode()) == row['source_text_sha256'],
                'milestone text: ' + row['id'])
        require((match[1].lower() == 'x') == row['legacy_checked'], 'legacy status: ' + row['id'])
    packages = json.loads((packet / 'WORK_PACKAGES.json').read_text())['packages']
    by_id = {p['id']: p for p in packages}
    require(len(packages) == len(by_id) == 44, '44 unique work packages')
    require(set(by_id) == {f'P{i:02d}' for i in range(44)}, 'package ID sequence')
    used = {p for row in rows for p in row['work_packages'] if p != 'ALL'}
    require(used == set(by_id), 'every package covered by roadmap')
    seen, visiting = set(), set()

    def visit(pid):
        require(pid in by_id, 'unknown dependency: ' + pid)
        require(pid not in visiting, 'cyclic acceptance dependency: ' + pid)
        if pid in seen:
            return
        visiting.add(pid)
        for dependency in by_id[pid]['acceptance_dependencies']:
            visit(dependency)
        visiting.remove(pid)
        seen.add(pid)

    for package in packages:
        visit(package['id'])
        for name in package['canonical_owners']:
            require((root / name).exists(), 'missing canonical owner: ' + name)
    plan = (root / 'docs/OPENCODE_IMPLEMENTATION_HANDOFF.md').read_text()
    for pid in by_id:
        require(plan.count('<a id="' + pid.lower() + '"></a>') == 1,
                'plan anchor: ' + pid)
    links = 0
    for name in manifest['documents']:
        path = root / name
        for target in re.findall(r'\]\(([^)]+)\)', path.read_text()):
            target = target.strip('<>')
            if re.match(r'^[a-zA-Z][\w+.-]*:', target):
                continue
            local, _, fragment = target.partition('#')
            local = unquote(local)
            destination = (path.parent / local) if local else path
            require(destination.exists(), 'broken link in ' + name + ': ' + target)
            if fragment.startswith('p') and re.fullmatch(r'p\d\d', fragment):
                require('<a id="' + fragment + '"></a>' in destination.read_text(),
                        'missing package anchor: ' + target)
            links += 1
    print(json.dumps({'result': 'PASS (metadata only; no Julia execution)',
                      'packet_files': len(actual), 'milestones': len(rows),
                      'legacy_checked': 23, 'open': 56, 'work_packages': len(packages),
                      'dependency_graph': 'acyclic', 'document_links_checked': links,
                      'source_equals_baseline': current_matches and not bool(source_status)}, indent=2))


if __name__ == '__main__':
    main()
