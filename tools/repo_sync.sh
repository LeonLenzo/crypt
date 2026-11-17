#!/bin/bash
# Keep the local clone and the Setonix clone from drifting, in the only two places it matters.
#
#   tools/repo_sync.sh require-current
#   tools/repo_sync.sh publish "<commit subject>" <path> [<path> ...]
#
# WHY THIS EXISTS. One directory tree has two authors. Code and decisions are written locally;
# the bytes (CDS, reads, indexes) are written on Setonix and never enter git. The trouble is
# the third class: the MANIFESTS that describe those bytes are small, tracked, and only the
# machine holding the bytes can write them truthfully. In one session that collided three
# times -- host_cds_manifest.tsv, reads_manifest.tsv and cds_fetch_manifest.tsv all existed on
# Setonix, were committed from a local copy, and then blocked the next merge -- and the
# .gitignore already records two earlier cases (predict_sample.tsv, ref_candidates.tsv) that
# failed the other way, as jobs dying on an input that never arrived.
#
# The fix is ONE AUTHOR PER FILE. A manifest is committed by the machine that produced it, so
# it reaches the other side through git rather than through scp, and no file has two writers
# for git to arbitrate between.
#
# `require-current` is the other half: a job that runs on a stale checkout produces results
# attributable to code nobody can identify afterwards. It is cheap to check and the failure it
# prevents is expensive.

set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO" || exit 1
BRANCH="${REPO_SYNC_BRANCH:-master}"

# Two jobs can finish within seconds of each other, and git takes a lock per index write.
git_locked() {
    local n=0
    while [ -f .git/index.lock ] && [ $n -lt 30 ]; do sleep 2; n=$((n + 1)); done
    git "$@"
}

require_current() {
    local head before
    before="$(git rev-parse --short HEAD)"
    if ! git fetch -q origin "$BRANCH" 2>/dev/null; then
        echo "repo_sync: WARNING cannot reach origin; running on $before unverified"
        return 0
    fi
    if git merge-base --is-ancestor "origin/$BRANCH" HEAD; then
        echo "repo_sync: current at $before"
        return 0
    fi
    # Behind, or diverged. A fast-forward is the only move a job may make on its own: anything
    # else means somebody has unpushed work here and a job must not resolve that by guessing.
    if git_locked merge --ff-only "origin/$BRANCH" >/dev/null 2>&1; then
        head="$(git rev-parse --short HEAD)"
        echo "repo_sync: updated $before -> $head"
        return 0
    fi
    echo "repo_sync: ABORT. This checkout is not a fast-forward of origin/$BRANCH."
    echo "  HEAD      $(git log --oneline -1 HEAD)"
    echo "  origin    $(git log --oneline -1 "origin/$BRANCH")"
    echo "  Resolve by hand. A job will not merge diverged history or overwrite local commits."
    return 1
}

publish() {
    local msg="$1"; shift
    [ $# -gt 0 ] || { echo "repo_sync: publish needs at least one path"; return 1; }

    local existing=()
    local p
    for p in "$@"; do
        if [ -e "$p" ]; then existing+=("$p"); else echo "repo_sync: no such path, skipped: $p"; fi
    done
    [ ${#existing[@]} -gt 0 ] || { echo "repo_sync: nothing to publish"; return 0; }

    git_locked add -- "${existing[@]}" || return 1
    if git diff --cached --quiet -- "${existing[@]}"; then
        echo "repo_sync: no change in ${existing[*]}, nothing to publish"
        git_locked reset -q -- "${existing[@]}" 2>/dev/null
        return 0
    fi

    # Only the named paths are committed. The tree carries other modifications (tables a job
    # regenerated, a manifest another job is mid-write) and they are not this commit's business.
    local body
    body="Written on $(hostname) by SLURM job ${SLURM_JOB_ID:-none}, published by the job that
produced it so the file has one author and reaches the other machine through git."
    if ! git_locked commit -q --only -m "$msg" -m "$body" -- "${existing[@]}"; then
        echo "repo_sync: commit failed"
        return 1
    fi
    echo "repo_sync: committed $(git rev-parse --short HEAD)  ${existing[*]}"

    if git push -q origin "HEAD:$BRANCH" 2>/dev/null; then
        echo "repo_sync: pushed"
        return 0
    fi
    # Rejected: origin moved while the job ran. --autostash because the tree is usually dirty
    # with other jobs' output, and rebase refuses to start otherwise.
    echo "repo_sync: push rejected, rebasing onto origin/$BRANCH"
    if git_locked pull -q --rebase --autostash origin "$BRANCH" && git push -q origin "HEAD:$BRANCH"; then
        echo "repo_sync: pushed after rebase"
        return 0
    fi
    echo "repo_sync: COULD NOT PUSH. The commit is here at $(git rev-parse --short HEAD) and"
    echo "  the file is on disk, but the other machine will not see it until this is resolved."
    git status --short | head -5
    return 1
}

case "${1:-}" in
    require-current) require_current ;;
    publish) shift; publish "$@" ;;
    *) echo "usage: repo_sync.sh require-current | publish \"<subject>\" <path>..."; exit 2 ;;
esac
