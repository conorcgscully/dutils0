#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${ROOT_DIR}/.env.local"

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "Missing ${ENV_FILE}."
  echo "Copy .env.local.example to .env.local and fill in your paths."
  exit 1
fi

# shellcheck disable=SC1090
source "${ENV_FILE}"

: "${MICROMAMBA_BIN:?MICROMAMBA_BIN is required in .env.local}"
: "${MICROMAMBA_ENV_PREFIX:?MICROMAMBA_ENV_PREFIX is required in .env.local}"
: "${NCYCLE_CORE_PATH:?NCYCLE_CORE_PATH is required in .env.local}"

"${MICROMAMBA_BIN}" run -p "${MICROMAMBA_ENV_PREFIX}" python -m pip install -e "${NCYCLE_CORE_PATH}"
"${MICROMAMBA_BIN}" run -p "${MICROMAMBA_ENV_PREFIX}" python -m pip install -e "${ROOT_DIR}"
"${MICROMAMBA_BIN}" run -p "${MICROMAMBA_ENV_PREFIX}" python -m pytest -q "${ROOT_DIR}/dutils0/tests"
