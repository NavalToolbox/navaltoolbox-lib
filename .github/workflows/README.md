# GitHub Actions Workflows

Ce répertoire contient les workflows CI/CD pour `navaltoolbox-lib`.

## Workflows

### 🔄 ci.yml - Intégration Continue

**Déclencheurs**: `push` et `pull_request` sur `main`

| Job | Description | Plateforme |
|-----|-------------|------------|
| `lint-rust` | `cargo fmt --check` + `cargo clippy` | Ubuntu |
| `test-rust` | `cargo test` | Ubuntu |
| `build-wheels` | Build wheels Python | Linux, macOS, Windows |
| `test-python` | Tests pytest | Python 3.9-3.14 |

### 📚 docs.yml - Documentation

**Déclencheurs**: `push` et `pull_request` sur `main`

| Job | Description |
|-----|-------------|
| `build` | Build Sphinx avec PyData theme |
| `deploy` | Déploiement sur GitHub Pages (main seulement) |

### 📦 publish.yml - Publication

**Déclencheurs**: 
- Tags `v*.*.*`
- `workflow_dispatch` (manuel)

| Job | Destination | Déclencheur |
|-----|-------------|-------------|
| `build-wheels` | - | Toujours |
| `build-sdist` | - | Toujours |
| `publish-testpypi` | TestPyPI | Manuel |
| `publish-pypi` | PyPI | Tag `v*.*.*` ou manuel |
| `publish-cratesio` | crates.io | Tag `v*.*.*` ou manuel |

## Secrets requis

| Secret | Description | Utilisé par |
|--------|-------------|-------------|
| `CARGO_REGISTRY_TOKEN` | Token crates.io | `publish-cratesio` |

> **Note**: La publication sur PyPI utilise le **Trusted Publishing** (OIDC), aucun token requis.

## Plateformes supportées

- 🐧 Linux (ubuntu-latest)
- 🍎 macOS (macos-latest)
- 🪟 Windows (windows-latest)

## Versions Python

- Python 3.9
- Python 3.10
- Python 3.11
- Python 3.12
- Python 3.13
- Python 3.14 (prerelease)

## Usage manuel

```bash
# Déclencher la publication sur TestPyPI
gh workflow run publish.yml -f environment=testpypi

# Déclencher la publication sur PyPI
gh workflow run publish.yml -f environment=pypi

# Déclencher la publication sur crates.io
gh workflow run publish.yml -f environment=cratesio
```
