#!/usr/bin/env bash
# Minimal end-to-end MCP handshake against immunopipe's stdio transport.
# Sends initialize + tools/list over stdio and prints the tool names.
# Expected: a JSON-RPC response listing immunopipe's tools, then "SMOKE OK".
#
# Usage:
#   bash scripts/mcp_smoke.sh                                   # `immunopipe` from PATH
#   IMMUNOPIPE_BIN=/path/to/.venv/bin/immunopipe bash scripts/mcp_smoke.sh
set -uo pipefail

IMMUNOPIPE_BIN="${IMMUNOPIPE_BIN:-immunopipe}"
SMOKE_STDOUT="${SMOKE_STDOUT:-/tmp/mcp_smoke_stdout.json}"
SMOKE_STDERR="${SMOKE_STDERR:-/tmp/mcp_smoke_stderr.log}"

# `timeout` is coreutils on Linux, `gtimeout` with brew on macOS; degrade gracefully.
TIMEOUT_BIN="$(command -v timeout || command -v gtimeout || true)"

PYTHON_BIN="${PYTHON:-}"
if [ -z "$PYTHON_BIN" ]; then
    if command -v python3 >/dev/null 2>&1; then PYTHON_BIN=python3; else PYTHON_BIN=python; fi
fi

echo "--- binary: $IMMUNOPIPE_BIN ---"
if ! command -v "$IMMUNOPIPE_BIN" >/dev/null 2>&1 && [ ! -x "$IMMUNOPIPE_BIN" ]; then
    echo "SMOKE FAILED — '$IMMUNOPIPE_BIN' is not an executable command."
    echo "This is a client-configuration problem, not a server problem: point \"command\" at the"
    echo "console script (\`command -v immunopipe\`) or at an absolute interpreter path that can"
    echo "import immunopipe."
    exit 127
fi

serve() {
    if [ -n "$TIMEOUT_BIN" ]; then
        "$TIMEOUT_BIN" 60 "$IMMUNOPIPE_BIN" mcp --transport stdio
    else
        "$IMMUNOPIPE_BIN" mcp --transport stdio
    fi
}

{
  printf '%s\n' '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"smoke","version":"0.0.1"}}}'
  printf '%s\n' '{"jsonrpc":"2.0","method":"notifications/initialized"}'
  printf '%s\n' '{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}'
} | serve 2>"$SMOKE_STDERR" > "$SMOKE_STDOUT"

rc=$?
echo "--- exit: $rc ---"

if [ -z "$PYTHON_BIN" ]; then
    if grep -q '"tools"' "$SMOKE_STDOUT" 2>/dev/null; then echo "SMOKE OK"; else
        echo "SMOKE FAILED (no python to parse; stdout was:)"; head -c 800 "$SMOKE_STDOUT"; exit 1
    fi
    exit 0
fi

"$PYTHON_BIN" - "$SMOKE_STDOUT" <<'PYEOF'
import json
import sys

path = sys.argv[1]
try:
    raw = open(path).read()
except OSError as exc:
    sys.exit("SMOKE FAILED — no stdout captured: %s" % exc)

frames, bad = [], []
for line in raw.splitlines():
    if not line.strip():
        continue
    try:
        frames.append(json.loads(line))
    except ValueError:
        bad.append(line)

if bad:
    sys.exit("SMOKE FAILED — stdout must carry only JSON-RPC frames, got non-JSON line(s): %r"
             % bad[:3])

init = next((f for f in frames if f.get("id") == 1), None)
tools = next((f for f in frames if f.get("id") == 2), None)
if init is None or tools is None:
    sys.exit("SMOKE FAILED — expected responses with id=1 (initialize) and id=2 (tools/list), got: %r"
             % [f.get("id") for f in frames])

info = init.get("result", {}).get("serverInfo", {})
print("server: %s %s (protocol %s)" % (info.get("name"), info.get("version"),
                                       init.get("result", {}).get("protocolVersion")))
names = [t["name"] for t in tools.get("result", {}).get("tools", [])]
if not names:
    sys.exit("SMOKE FAILED — tools/list returned no tools")
print("tools advertised (%d):" % len(names))
for name in names:
    print("  -", name)

# Notifications must not be answered; anything with a null id is a protocol error.
stray = [f for f in frames if f.get("id") not in (1, 2)]
if stray:
    sys.exit("SMOKE FAILED — unexpected frame(s) beyond the two responses: %r" % stray)
PYEOF

parse_rc=$?
if [ "$parse_rc" -eq 0 ]; then
    echo "SMOKE OK"
else
    echo "SMOKE FAILED — stdout was:"
    head -c 800 "$SMOKE_STDOUT"
    echo "--- stderr ---"
    tail -20 "$SMOKE_STDERR"
    exit 1
fi
