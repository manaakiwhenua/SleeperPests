#!/usr/bin/env sh
set -eu
URL='https://raw.githubusercontent.com/manaakiwhenua/SleeperPests/6e4b9f032bec9001bf7241a9d59625894c776f9b/INApest%20INApestMeta%20core%20function%20code/INApestAnalytical.R'
OUT="$(dirname "$0")/R/INApestAnalytical.R"
curl -L "$URL" -o "$OUT"
EXPECTED='49238e62a9f4c99445072d08356ca33081c60b37f1a4b268b328794990d163c5'
ACTUAL=$(sha256sum "$OUT" | awk '{print $1}')
[ "$ACTUAL" = "$EXPECTED" ] || { echo "SHA-256 mismatch: $ACTUAL" >&2; exit 1; }
echo "PASS INApestAnalytical.R $ACTUAL"
