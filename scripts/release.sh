#!/usr/bin/env bash
# Copyright (c) 2015-2026 Claudio Satriano
# SPDX-License-Identifier: GPL-3.0-or-later

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: scripts/release.sh VERSION [--dry-run]

VERSION input:
    - You can pass "1.7" or "v1.7" (both are accepted)
    - Recommended input is "1.7"
    - The created git tag is always "vVERSION" (example: "v1.7")

Create a release commit+tag and start the next development cycle:
1) Replace CHANGELOG "[unreleased]" heading with "[VERSION] - YYYY-MM-DD",
   and replace the "[unreleased]:" footer link with the new version link
2) Update version in .zenodo.json and CITATION.cff
3) Commit as "Release vVERSION" and create annotated tag
4) Re-add "[unreleased]" section and footer link to CHANGELOG
5) Commit as "Start next development cycle"

Options:
  --dry-run  Validate and print the actions without changing files
  -h, --help Show this help

Note:
    Run help as "scripts/release.sh -h" (not just "-h").

Examples:
    scripts/release.sh 1.7
    scripts/release.sh v1.7
    scripts/release.sh 1.7.1
    scripts/release.sh v1.7 --dry-run
EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -eq 0 ]]; then
    usage
    exit 0
fi

VERSION="$1"
shift

DRY_RUN=false
for arg in "$@"; do
    case "$arg" in
        --dry-run)
            DRY_RUN=true
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $arg"
            ;;
    esac
done

VERSION="${VERSION#v}"
[[ "$VERSION" =~ ^[0-9]+(\.[0-9]+){1,2}$ ]] || die "Invalid version '$VERSION'. Expected N.N or N.N.N (optional 'v' prefix is allowed)"

require_cmd git
require_cmd awk
require_cmd mktemp

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || die "Not inside a git repository"
cd "$REPO_ROOT"

CHANGELOG="CHANGELOG.md"
ZENODO_JSON=".zenodo.json"
CITATION_CFF="CITATION.cff"
[[ -f "$CHANGELOG" ]] || die "Missing $CHANGELOG"
[[ -f "$ZENODO_JSON" ]] || die "Missing $ZENODO_JSON"
[[ -f "$CITATION_CFF" ]] || die "Missing $CITATION_CFF"

if ! git diff-index --quiet HEAD --; then
    die "Working tree is not clean. Commit or stash changes before releasing."
fi

git fetch --tags --quiet

TAG="v$VERSION"
if git show-ref --tags --verify --quiet "refs/tags/$TAG"; then
    die "Tag $TAG already exists"
fi

if ! grep -q '^## \[unreleased\]$' "$CHANGELOG"; then
    die "Could not find '## [unreleased]' in $CHANGELOG"
fi

PREV_VERSION=$(grep -m1 '^## \[[0-9]' "$CHANGELOG" | sed 's/^## \[\([^]]*\)\].*/\1/')
[[ -n "$PREV_VERSION" ]] || die "Could not determine previous version from $CHANGELOG"

REPO_URL="$(git remote get-url origin | sed 's/\.git$//; s|^git@github\.com:|https://github.com/|')"

RELEASE_DATE="$(date +%Y-%m-%d)"
COMMIT_MSG="Release $TAG"

echo "Release configuration"
echo "  version: $VERSION"
echo "  tag: $TAG"
echo "  date: $RELEASE_DATE"
echo "  dry-run: $DRY_RUN"

if $DRY_RUN; then
    echo
    echo "Dry run checks passed. Planned actions:"
    echo "  1) Replace CHANGELOG '[unreleased]' heading with version and date, replace '[unreleased]:' footer link with new version link"
    echo "  2) Update version in $ZENODO_JSON"
    echo "  3) Update version and date-released in $CITATION_CFF"
    echo "  4) git add $CHANGELOG $ZENODO_JSON $CITATION_CFF"
    echo "  5) git commit -m \"$COMMIT_MSG\""
    echo "  6) git tag -a $TAG -m \"$COMMIT_MSG\""
    echo "  7) Re-add '[unreleased]' section and footer link to $CHANGELOG"
    echo "  8) git add $CHANGELOG"
    echo "  9) git commit -m \"Start next development cycle\""
    BRANCH="$(git rev-parse --abbrev-ref HEAD)"
    echo "  10) git push origin $BRANCH $TAG"
    exit 0
fi

TMP_FILE="$(mktemp)"
trap 'rm -f "$TMP_FILE"' EXIT

awk -v ver="$VERSION" -v d="$RELEASE_DATE" -v prev="$PREV_VERSION" \
    -v repo="$REPO_URL" '
BEGIN { heading_done=0; link_done=0 }
!heading_done && /^## \[unreleased\]$/ {
    print "## [" ver "] - " d
    heading_done=1
    next
}
!link_done && /^\[unreleased\]:/ {
    print "[" ver "]: " repo "/compare/v" prev "...v" ver
    link_done=1
    next
}
{ print }
END {
    if (!heading_done) { exit 2 }
    if (!link_done) { exit 3 }
}
' "$CHANGELOG" > "$TMP_FILE" || die "Failed to update $CHANGELOG"

mv "$TMP_FILE" "$CHANGELOG"
trap - EXIT

TMP_FILE="$(mktemp)"
trap 'rm -f "$TMP_FILE"' EXIT

awk -v ver="$VERSION" '
BEGIN { done=0 }
{
    if (!done && $0 ~ /^[[:space:]]*"version"[[:space:]]*:[[:space:]]*"[^"]+",?[[:space:]]*$/) {
        sub(/"version"[[:space:]]*:[[:space:]]*"[^"]+"/, "\"version\": \"" ver "\"")
        done=1
    }
    print $0
}
END {
    if (!done) {
        exit 2
    }
}
' "$ZENODO_JSON" > "$TMP_FILE" || die "Failed to update version in $ZENODO_JSON"

mv "$TMP_FILE" "$ZENODO_JSON"
trap - EXIT

TMP_FILE="$(mktemp)"
trap 'rm -f "$TMP_FILE"' EXIT

awk -v ver="$VERSION" -v d="$RELEASE_DATE" '
BEGIN { version_done=0; date_done=0 }
{
    if (!version_done && $0 ~ /^[[:space:]]*version:[[:space:]]*"[^"]+"[[:space:]]*$/) {
        print "version: \"" ver "\""
        version_done=1
    } else if (!date_done && $0 ~ /^[[:space:]]*date-released:[[:space:]]*"[^"]+"[[:space:]]*$/) {
        print "date-released: \"" d "\""
        date_done=1
    } else {
        print $0
    }
}
END {
    if (!version_done || !date_done) {
        exit 2
    }
}
' "$CITATION_CFF" > "$TMP_FILE" || die "Failed to update version/date-released in $CITATION_CFF"

mv "$TMP_FILE" "$CITATION_CFF"
trap - EXIT

git add "$CHANGELOG" "$ZENODO_JSON" "$CITATION_CFF"
git commit -m "$COMMIT_MSG"
git tag -a "$TAG" -m "$COMMIT_MSG"

TMP_FILE="$(mktemp)"
trap 'rm -f "$TMP_FILE"' EXIT

awk -v ver="$VERSION" -v repo="$REPO_URL" '
!heading_done && /^## \[[0-9]/ {
    print "## [unreleased]"
    print ""
    heading_done=1
}
!link_done && /^\[[0-9]/ {
    print "[unreleased]: " repo "/compare/v" ver "...HEAD"
    link_done=1
}
{ print }
' "$CHANGELOG" > "$TMP_FILE" || die "Failed to re-add [unreleased] section to $CHANGELOG"

mv "$TMP_FILE" "$CHANGELOG"
trap - EXIT

git add "$CHANGELOG"
git commit -m "Start next development cycle"

echo
echo "Done. Release commit+tag and next-cycle commit created locally."
echo "Push manually with:"
echo "  git push origin $(git rev-parse --abbrev-ref HEAD) $TAG"
