---
name: vibe
description: Use when doing scientific research, literature review, hypothesis validation, gap finding, or any research task requiring rigorous evidence tracking and adversarial review. Activates the Vibe Science research engine with OTAE loops, Reviewer 2, quality gates, and serendipity tracking.
---

# Vibe Science — AI-Native Research Engine

> Think first, analyse second. Find what has NOT been done, then verify it with data.

## What This Does

Vibe Science turns Claude into a disciplined research agent. Every claim passes through adversarial review (Reviewer 2), quality gates (34 total), and confounder testing. Only what survives gets published.

**Core loop:** Literature gap → hypothesis → targeted analysis → validation → claim.

## Commands

| Command | Purpose |
|---------|---------|
| `/vibe-science:init` | Initialize a new research session |
| `/vibe-science:loop` | Execute one OTAE research cycle |
| `/vibe-science:search` | Literature search across Scopus, PubMed, OpenAlex |
| `/vibe-science:reviewer2` | Invoke adversarial review |
| `/vibe-science:start` | Full session initialization with SOLO/TEAM choice |

## How It Works

1. **OBSERVE** — Read current state, identify what's changed
2. **THINK** — Plan the highest-value next action
3. **ACT** — Execute ONE action (search, analyze, compute)
4. **EVALUATE** — Extract claims, score confidence, check gates

Every cycle crystallizes state to files. Reviewer 2 challenges all major findings. Nothing advances without evidence.

## Key Principles

- **NO DATA = NO GO** — No thesis without evidence from data
- **Every claim has a confidence score** — Computed formula, not vibes
- **Reviewer 2 is co-pilot** — Can VETO, REDIRECT, and FORCE re-investigation
- **Gates block progress** — 34 quality gates, 8 schema-enforced
- **Serendipity is tracked** — Unexpected findings scored and preserved
- **Everything persists to files** — If it's not in a file, it doesn't exist

## Full Methodology

The complete methodology (OTAE-Tree, R2 Ensemble, Evidence Engine, 21 protocols, 34 gates) is in the root `SKILL.md` (~1,300 lines). Load it progressively:

1. **Start here** — This file gives you the commands and principles
2. **For deep methodology** — Read `SKILL.md` at the plugin root
3. **For specific protocols** — Read from `protocols/` directory
4. **For gate definitions** — Read from `gates/gates.md`
5. **For schemas** — Read from `schemas/` directory

## Plugin Enforcement

When installed as a Claude Code plugin, Vibe Science adds code-level enforcement:

- **SessionStart hook** — Injects research context (~700 tokens), creates session in SQLite
- **PostToolUse hook** — Gate enforcement (exit code 2 = BLOCK), auto-logging, permission checks
- **Stop hook** — Blocks session end if unreviewed claims exist (LAW 4)
- **SQLite persistence** — Cross-session memory, R2 calibration, audit trail

## Important Rules

- **Web searches MUST be INLINE** — Use WebSearch/WebFetch directly, NOT via background sub-agents (Task tool). Sub-agents cannot inherit web permissions.
- **Scientific skills invoked directly** — Use the Skill tool for PubMed, GEO, OpenAlex, etc.
- **Template-first pattern** — Always read templates before creating files. Downstream hooks depend on template structure.
