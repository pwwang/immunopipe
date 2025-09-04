#!/usr/bin/env python3
"""Example script demonstrating immunopipe MCP server functionality."""

import asyncio
# import json
# from pathlib import Path

from .server import McpServer, SimpleMCPClient


async def demonstrate_mcp_tools():
    """Demonstrate the MCP server tools."""
    print("🚀 Starting Immunopipe MCP Server Demonstration\n")

    # Create server and client
    server = McpServer()
    client = SimpleMCPClient(server)

    try:
        # Initialize connection
        print("1. Initializing connection...")
        init_response = await client.initialize()
        if "result" in init_response:
            print(f"   ✓ Connected to: {init_response['result']['name']}")
            print(f"   ✓ Version: {init_response['result']['version']}")

        # List available tools
        print("\n2. Listing available tools...")
        tools = await client.list_tools()
        print(f"   ✓ Found {len(tools)} tools:")
        for i, tool in enumerate(tools[:5]):  # Show first 5
            print(f"     {i + 1}. {tool['name']}: {tool['description']}")
        if len(tools) > 5:
            print(f"     ... and {len(tools) - 5} more tools")

        # Generate a basic template
        print("\n3. Generating basic configuration template...")
        result = await client.call_tool("generate_basic_template")
        if result.get("result"):
            template = result["result"]["content"]
            print("   ✓ Basic template generated:")
            print("   " + "\n   ".join(template.split("\n")[:8]) + "\n   ...")

        # List pipeline options
        print("\n4. Discovering pipeline options...")
        result = await client.call_tool("list_pipeline_options")
        if result.get("result", {}).get("content"):
            options = result["result"]["content"]
            print(f"   ✓ Found {len(options)} pipeline options:")
            for name, details in list(options.items())[:3]:
                print(f"     - {name}: {details.get('description', 'No description')}")

        # List processes
        print("\n5. Discovering available processes...")
        result = await client.call_tool("list_processes")
        if result.get("result", {}).get("content"):
            processes = result["result"]["content"]
            print(f"   ✓ Found {len(processes)} processes:")
            for name, details in list(processes.items())[:3]:
                print(f"     - {name}: {details.get('description', 'No description')}")

        # Generate TCR analysis template
        print("\n6. Generating TCR analysis template...")
        result = await client.call_tool("generate_tcr_template")
        if result.get("result"):
            template = result["result"]["content"]
            print("   ✓ TCR template generated:")
            print("   " + "\n   ".join(template.split("\n")[:6]) + "\n   ...")

        # Get help for configuration generation
        print("\n7. Getting help for TCR analysis configuration...")
        result = await client.call_tool(
            "help_generate_config",
            {"description": "I want to analyze single-cell TCR data with clustering"},
        )
        if result.get("result", {}).get("content", {}).get("suggestions"):
            suggestions = result["result"]["content"]["suggestions"]
            print("   ✓ Configuration suggestions:")
            for suggestion in suggestions:
                print(f"     - {suggestion}")

        # Suggest processes for TCR analysis
        print("\n8. Suggesting processes for TCR analysis...")
        result = await client.call_tool("suggest_processes", {"analysis_type": "tcr"})
        if result.get("result", {}).get("content", {}).get("suggested_processes"):
            processes = result["result"]["content"]["suggested_processes"]
            print(f"   ✓ Suggested {len(processes)} processes for TCR analysis:")
            print(f"     {', '.join(processes)}")

        # Generate a complete configuration
        print("\n9. Generating complete configuration...")
        pipeline_opts = {"name": "demo_analysis", "outdir": "./demo_output"}
        process_configs = {
            "SampleInfo": {"in": {"infile": ["sample_info.txt"]}},
            "TOrBCellSelection": {"envs": {"cell_type": "T"}},
        }
        gbatch_opts = {"project": "demo-project", "region": "us-central1"}

        result = await client.call_tool(
            "generate_full_config",
            {
                "pipeline_options": pipeline_opts,
                "processes": process_configs,
                "gbatch_options": gbatch_opts,
                "description": "Demo configuration file",
            },
        )

        if result.get("result"):
            config = result["result"]["content"]
            print("   ✓ Complete configuration generated:")
            print("   " + "\n   ".join(config.split("\n")[:10]) + "\n   ...")

        print("\n🎉 MCP Server demonstration completed successfully!")
        print(
            "    The server provides hierarchical tools for composing "
            "immunopipe configurations"
        )
        print("    using natural language descriptions and structured parameters.")

    except Exception as e:
        print(f"❌ Error during demonstration: {e}")
        return False

    return True


async def interactive_demo():
    """Interactive demonstration of MCP tools."""
    print("🎯 Interactive Immunopipe MCP Demo")
    print("Available commands:")
    print("  1. basic - Generate basic template")
    print("  2. tcr - Generate TCR template")
    print("  3. help <description> - Get configuration help")
    print("  4. processes <type> - Suggest processes")
    print("  5. options - List pipeline options")
    print("  6. quit - Exit")

    server = McpServer()
    client = SimpleMCPClient(server)

    await client.initialize()

    while True:
        try:
            user_input = input("\n> ").strip()

            if user_input == "quit":
                break
            elif user_input == "basic":
                result = await client.call_tool("generate_basic_template")
                if result.get("result"):
                    print("\nBasic Template:")
                    print(result["result"]["content"])
            elif user_input == "tcr":
                result = await client.call_tool("generate_tcr_template")
                if result.get("result"):
                    print("\nTCR Template:")
                    print(result["result"]["content"])
            elif user_input.startswith("help "):
                description = user_input[5:]
                result = await client.call_tool(
                    "help_generate_config", {"description": description}
                )
                if result.get("result", {}).get("content", {}).get("suggestions"):
                    print("\nSuggestions:")
                    for suggestion in result["result"]["content"]["suggestions"]:
                        print(f"  - {suggestion}")
            elif user_input.startswith("processes "):
                analysis_type = user_input[10:]
                result = await client.call_tool(
                    "suggest_processes", {"analysis_type": analysis_type}
                )
                if (
                    result.get("result", {})
                    .get("content", {})
                    .get("suggested_processes")
                ):
                    processes = result["result"]["content"]["suggested_processes"]
                    print(f"\nSuggested processes for {analysis_type}:")
                    print(f"  {', '.join(processes)}")
            elif user_input == "options":
                result = await client.call_tool("list_pipeline_options")
                if result.get("result", {}).get("content"):
                    options = result["result"]["content"]
                    print(f"\nPipeline options ({len(options)} total):")
                    for name, details in list(options.items())[:10]:
                        print(
                            f"  {name}: {details.get('description', 'No description')}"
                        )
            else:
                print("Unknown command. Type 'quit' to exit.")

        except KeyboardInterrupt:
            break
        except Exception as e:
            print(f"Error: {e}")

    print("\nGoodbye!")


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "interactive":
        asyncio.run(interactive_demo())
    else:
        asyncio.run(demonstrate_mcp_tools())
