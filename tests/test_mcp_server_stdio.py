"""Protocol-level regression test for the immunopipe MCP server's stdio transport.

This guards the R1.3 defect: the server answered the `notifications/initialized`
notification with a JSON-RPC error whose `id` was null. That frame is not a valid
MCP `JSONRPCMessage` (RequestId must be a string or integer), so a stock MCP client
(e.g. the official Python SDK) reports a parse failure on every connection.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent

INITIALIZE = (
    '{"jsonrpc":"2.0","id":1,"method":"initialize",'
    '"params":{"protocolVersion":"2024-11-05","capabilities":{},'
    '"clientInfo":{"name":"pytest","version":"0.0.1"}}}'
)
INITIALIZED = '{"jsonrpc":"2.0","method":"notifications/initialized"}'
TOOLS_LIST = '{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}'

EXPECTED_TOOLS = {
    "list_pipeline_options",
    "list_processes",
    "list_gbatch_options",
    "generate_full_config",
    "configure_immunopipe",
}


def _run_stdio_session(payload: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "immunopipe", "mcp", "--transport", "stdio"],
        input=payload,
        capture_output=True,
        text=True,
        timeout=120,
        cwd=str(REPO_ROOT),
    )


def _frames(stdout: str) -> list:
    """Every non-empty stdout line must be a JSON-RPC message."""
    frames = []
    for line in stdout.splitlines():
        if not line.strip():
            continue
        frames.append(json.loads(line))  # raises if stdout is polluted
    return frames


def test_stdio_handshake_answers_only_requests():
    proc = _run_stdio_session("\n".join([INITIALIZE, INITIALIZED, TOOLS_LIST]) + "\n")

    assert proc.returncode == 0, proc.stderr
    frames = _frames(proc.stdout)

    # Exactly two responses — the two requests. The notification must not be answered.
    assert [frame.get("id") for frame in frames] == [1, 2], frames
    assert not [frame for frame in frames if frame.get("id") is None], frames
    assert not [frame for frame in frames if "error" in frame], frames

    init, tools = frames
    assert init["result"]["protocolVersion"] == "2024-11-05"
    assert init["result"]["serverInfo"]["name"] == "immunopipe-mcp"
    assert EXPECTED_TOOLS <= {tool["name"] for tool in tools["result"]["tools"]}


def test_stdio_stdout_carries_no_logging():
    """Logging goes to stderr; stdout must stay pure JSON-RPC."""
    proc = _run_stdio_session(INITIALIZE + "\n")

    assert proc.returncode == 0, proc.stderr
    frames = _frames(proc.stdout)
    assert len(frames) == 1
    assert "starting" not in proc.stdout
    assert "Immunopipe MCP server starting" in proc.stderr
