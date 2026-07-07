#!/bin/bash
# Install "Chrome for Testing" (CfT) into the user's XDG data directory, no root
# required. CfT is only needed to run the headless-browser GUI end-to-end tests
# (tests/testthat/test-gui-e2e.R), which drive the Shiny app via chromote +
# shinytest2. Those tests skip automatically when no Chrome/Chromium is found, so
# this is optional and only relevant for local development.
#
# Usage:
#   bash misc/scripts/install-chrome-for-testing.sh
#
# The script prints the CHROMOTE_CHROME export line to add to your shell so that
# chromote (and therefore shinytest2) picks up the installed binary. On Ubuntu
# 22.04 the required shared libraries ship with the base system; if Chrome fails
# to start, install its runtime dependencies (needs root), e.g.:
#   sudo apt-get install -y libnss3 libatk-bridge2.0-0 libgbm1 libasound2 \
#                           libxkbcommon0 libgtk-3-0
set -euo pipefail

DEST="${CFT_DIR:-$HOME/.local/share/chrome-for-testing}"
JSON_URL="https://googlechromelabs.github.io/chrome-for-testing/last-known-good-versions-with-downloads.json"

command -v curl >/dev/null || { echo "curl is required" >&2; exit 1; }
command -v unzip >/dev/null || { echo "unzip is required" >&2; exit 1; }

echo "Resolving latest stable Chrome-for-Testing (linux64) ..."
ZIP_URL="$(curl -fsSL "$JSON_URL" \
  | grep -o '"platform":"linux64","url":"[^"]*chrome-linux64.zip"' \
  | grep -m1 -o 'https[^"]*')"
[ -n "$ZIP_URL" ] || { echo "Could not resolve download URL" >&2; exit 1; }
echo "  $ZIP_URL"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
echo "Downloading ..."
curl -fsSL -o "$TMP/cft.zip" "$ZIP_URL"

echo "Installing into $DEST ..."
mkdir -p "$DEST"
rm -rf "$DEST/chrome-linux64"
unzip -q "$TMP/cft.zip" -d "$DEST"

CHROME="$DEST/chrome-linux64/chrome"
"$CHROME" --version

cat <<EOF

Chrome for Testing installed. To let chromote / shinytest2 find it, export:

    export CHROMOTE_CHROME="$CHROME"

Add that to your shell profile to make it permanent. Then run the GUI e2e tests:

    R -e 'devtools::test(filter = "gui-e2e")'

(requires the 'chromote' and 'shinytest2' Suggests packages and a JDK for rcdk).
EOF
