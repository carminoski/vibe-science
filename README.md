<p align="center">
  <img src="logos/logo-v6.0.svg" alt="Vibe Science" width="700">
</p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.18665031"><img src="https://zenodo.org/badge/1148022920.svg" alt="DOI"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-Apache_2.0-blue.svg" alt="License"></a>
  <img src="https://img.shields.io/badge/version-6.0.0-purple.svg" alt="Version">
  <img src="https://img.shields.io/badge/node-%3E%3D18.0.0-brightgreen.svg" alt="Node">
  <img src="https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey.svg" alt="Platform">
</p>

# Vibe Science

> An AI-native research engine that loops until discovery — with adversarial review, quality gates, serendipity tracking, and plugin-enforced integrity.

Vibe Science turns Claude Code into a disciplined research agent. Instead of letting the AI rush to conclusions, it forces every claim through adversarial review ("Reviewer 2"), 34 quality gates, and confounder testing. Only what survives gets published.

**Field-tested over 21 CRISPR research sprints** — caught a claim with p < 10^-100 whose sign reversed under propensity matching. Without Vibe Science, it would have been published.

---

## What Problem Does This Solve?

AI agents are dangerous in science. Not because they hallucinate — that's the easy problem. The dangerous problem is that they find **real patterns** in **real data** and construct **plausible narratives** around them, without ever asking: *"What if this is an artifact?"*

Over 21 CRISPR sprints, we watched the agent celebrate a result with OR=2.30 and p < 10^-100. After propensity matching, the sign reversed. Without a structural adversary, this would have been published as a finding.

**The solution:** embed a "Reviewer 2" whose ONLY job is to destroy claims. Only what survives both builder and destroyer advances.

---

## Why a Plugin, Not Just a Skill?

Vibe Science started as a **skill** (a system prompt). Versions v3.5 through v5.5 were prompt-only: the agent was *told* to run quality gates, use Reviewer 2, and track serendipity. It worked — sometimes. But four subsystems were bypassable because they relied on the agent's voluntary compliance:

| Subsystem | As Skill (v5.5) | As Plugin (v6.0) |
|-----------|-----------------|------------------|
| **Reviewer 2** | "Please review this claim" | Hook blocks session end if claims are unreviewed |
| **Quality Gates** | "Check gate DQ4 before proceeding" | Exit code 2 = tool action blocked until gate passes |
| **Research Logging** | "Write to PROGRESS.md" | Auto-logged to SQLite after every tool use |
| **Memory Recall** | "Read STATE.md at session start" | Hook injects ~700 tokens of context automatically |

The plugin wraps the skill in **code-level enforcement**: 5 lifecycle hooks, a gate engine, a permission engine, and SQLite persistence across sessions. The skill (SKILL.md, ~1,300 lines) remains the methodology brain. The plugin is the body that makes compliance structural, not voluntary.

### Evolution

| Version | Codename | Key Innovation |
|---------|----------|----------------|
| **v3.5** | TERTIUM DATUR | Foundation: OTAE loop, R2 double-pass, field-tested (21 sprints) |
| **v4.0** | ARBOR VITAE | Tree search over hypotheses |
| **v4.5** | ARBOR VITAE (Pruned) | Phase 0 brainstorm, R2 6 modes |
| **v5.0** | IUDEX | Seeded Fault Injection, blind-first pass, R3 judge |
| **v5.5** | ORO | 7 data quality gates, R2 INLINE mode |
| **v6.0** | NEXUS | **Plugin enforcement**, SQLite, 5 hooks, permissions |

---

## Requirements

| Requirement | Version | Why | Check |
|-------------|---------|-----|-------|
| **Node.js** | >= 18.0.0 | Runtime for hooks and plugin scripts | `node --version` |
| **Claude Code** | >= 1.0.33 | Plugin host | `claude --version` |
| **Git** | any | Clone the repo | `git --version` |
| **C++ Build Tools** | — | Required by `better-sqlite3` (native SQLite binding) | See below |

**C++ Build Tools by platform:**
- **Windows:** Visual Studio Build Tools with "Desktop development with C++" workload, or `npm install -g windows-build-tools`
- **macOS:** `xcode-select --install`
- **Linux:** `sudo apt install build-essential` (Debian/Ubuntu) or equivalent

**Optional:**
- **Bun** — faster alternative to Node for hook execution
- **Python** — only needed if your research involves Python-based analysis

---

## Quick Start

```bash
# 1. Clone and install
git clone https://github.com/th3vib3coder/vibe-science.git
cd vibe-science
npm install

# 2. Launch Claude Code with the plugin
claude --plugin-dir .

# 3. Start a research session
/vibe
```

On first startup, the SessionStart hook auto-creates `~/.vibe-science/`, initializes the SQLite database, and injects research context.

---

## Installation Methods

### Marketplace (Recommended)

```bash
/plugin marketplace add th3vib3coder/vibe-science
/plugin install vibe-science@vibe-science
# Restart Claude Code
```

### `--plugin-dir` (Quick Test)

```bash
claude --plugin-dir /path/to/vibe-science
```

### Global Settings (Permanent)

Add to `~/.claude/settings.json` (or `%USERPROFILE%\.claude\settings.json` on Windows):

```json
{
  "plugins": ["/absolute/path/to/vibe-science"]
}
```

### Project-Level

Add to your project's `.claude/settings.json` to load only in that project.

---

## What Does It Do

When active, every Claude Code session becomes a structured research session:

```
Research question → Brainstorm (Phase 0)
                         ↓
                    OTAE Loop (repeats):
                      OBSERVE  → Read current state
                      THINK    → Plan next action
                      ACT      → Execute ONE action
                      EVALUATE → Extract claims, score confidence
                         ↓
                    Reviewer 2 (adversarial review)
                         ↓
                    Only surviving claims advance
```

**What you'll notice:**
- Every claim gets a confidence score (0-1) with a mathematical formula
- Reviewer 2 assumes every claim is wrong and demands evidence
- 34 quality gates block progress — you can't skip steps
- State files are created automatically (STATE.md, PROGRESS.md, CLAIM-LEDGER.md)
- Serendipity is tracked — unexpected findings get scored and preserved
- Everything persists to SQLite — cross-session memory, R2 calibration, audit trail

---

## Architecture at a Glance

| Layer | What It Does | How |
|-------|-------------|-----|
| **Skill** (prompt-level) | Guides *what the agent thinks* | OTAE loop, R2 Ensemble, 11 Laws, 21 protocols |
| **Plugin** (code-level) | Controls *what the agent can do* | 5 hooks, gate engine, permissions, SQLite |

```
  SKILL (the brain)                PLUGIN (the body)
  ─────────────────                ──────────────────
  OTAE loop                        5 lifecycle hooks
  R2 Ensemble (7 modes)            Gate Engine (34 gates)
  11 Constitutional Laws           Permission Engine (6 roles)
  21 protocols                     Research Spine (auto-log)
  Evidence Engine                  SQLite (11 tables)
  Serendipity Engine               Silent Observer + R2 Calibration

  Guides REASONING                 Enforces BEHAVIOR
```

For deep technical details, see [ARCHITECTURE.md](ARCHITECTURE.md) and [SKILL.md](SKILL.md).

---

## Repository Structure

```
vibe-science/
├── .claude-plugin/          ← Plugin manifests
│   ├── plugin.json
│   └── marketplace.json
├── skills/vibe/SKILL.md     ← Skill entry point (auto-discovered)
├── commands/                ← Slash commands (auto-discovered)
│   ├── start.md, init.md, loop.md, search.md, reviewer2.md
├── agents/reviewer2.md      ← R2 subagent (auto-discovered)
├── hooks/hooks.json         ← 4 lifecycle hooks (auto-discovered)
├── plugin/                  ← Enforcement engine (~6,600 LOC)
│   ├── scripts/             ← Hook implementations (JS)
│   ├── lib/                 ← Engine modules
│   └── db/                  ← Schema + config
├── SKILL.md                 ← Full methodology (~1,300 lines)
├── CLAUDE.md                ← Project constitution
├── protocols/               ← 21 methodology protocols
├── gates/                   ← 34 quality gate specs
├── schemas/                 ← 9 JSON validation schemas
├── assets/                  ← Fault taxonomy, rubrics
└── archive/                 ← Historical versions (v3.5 → v5.5)
```

---

## Troubleshooting

| Problem | Cause | Fix |
|---------|-------|-----|
| Plugin not found | Wrong path | Verify `.claude-plugin/plugin.json` exists at the path |
| `npm install` fails (Windows) | `better-sqlite3` needs C++ | Install VS Build Tools with C++ workload |
| `npm install` fails (macOS) | Missing Xcode tools | `xcode-select --install` |
| Hooks don't fire | Not loaded as plugin | Use `--plugin-dir`, marketplace, or settings.json |
| SQLite errors | DB corruption | Delete `~/.vibe-science/db/` and restart |
| Embedding worker fails | Missing ONNX runtime | Non-critical — falls back to keyword search |

---

## Using Without Claude Code

You can use the methodology with any LLM: upload `SKILL.md` as a system prompt, plus `protocols/` and `gates/gates.md`. Note: without the plugin, gates are prompt-enforced only (voluntary compliance).

---

## Citation

> Vibe Science Contributors (2026). *Vibe Science: an AI-native research engine with adversarial review and serendipity tracking.* GitHub: [th3vib3coder/vibe-science](https://github.com/th3vib3coder/vibe-science) · DOI: [10.5281/zenodo.18665031](https://doi.org/10.5281/zenodo.18665031)

```bibtex
@software{vibe_science_2026,
  title     = {Vibe Science: AI-native research with adversarial review and serendipity tracking},
  author    = {{Vibe Science Contributors}},
  year      = {2026},
  version   = {6.0.0},
  url       = {https://github.com/th3vib3coder/vibe-science},
  doi       = {10.5281/zenodo.18665031},
  license   = {Apache-2.0}
}
```

## License

Apache 2.0 — see [LICENSE](LICENSE).

## Authors

**Carmine Russo, Elisa Bertelli (MD)**

---

*Built with Claude Code · Powered by Claude Opus · Made with adversarial love*
