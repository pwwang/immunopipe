# AI-powered chat on reports

Every HTML report immunopipe generates includes a chat panel that answers questions
about the page it is on. It is a separate mechanism from the
[MCP server](mcp-server.md): the MCP server helps an agent *write* a configuration, while
the chat helps you *read* the results of a run that has finished.

## Opening it

Open any report - from your own run, or any report in the [gallery](gallery.md) - and use
the chat panel. There is nothing to install or connect first.

## What it can see

The panel reads the page it is displayed on, plus any file you attach to a message. So on
a cell-type annotation page it can be asked which markers define a cluster; on an
enrichment page, what a term means. It does not read your input files unless you attach
them.

## The model is yours

You supply the model and the API key. They are entered in the panel and stay in your
browser; nothing is sent anywhere until you provide them, and immunopipe ships no key of
its own. Any chat model will do - the panel needs no function or tool calling, so models
that support neither also work.

## Where the panel comes from

It is provided by the report renderer,
[pipen-report](https://github.com/pwwang/pipen-report) 1.2 and later, which integrates the
[page-chat](https://github.com/pwwang/page-chat.js) library.
