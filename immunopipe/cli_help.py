from argx import ArgumentParser
from pipen_args import parser as pipen_parser

from . import processes
from .pipeline import Immunopipe


def main(args: list[str]) -> None:
    """Show help message for specific processes.

    Args:
        args: The command line arguments.
    """
    pipe = Immunopipe()
    pipe.build_proc_relationships()
    parser = ArgumentParser(
        description="Show help message for specific processes.",
        prog="immunopipe help",
    )
    help_message = ["The process to show help message for."]
    help_message.append("---")
    help_message.append("Processes:")
    for proc in pipe.procs:
        doc = proc.__doc__ or proc.__bases__[0].__doc__
        help_message.append(f"* {proc.name}: {doc.strip().splitlines()[0]}")

    parser.add_argument(
        "topic",
        nargs="?",
        choices=[proc.name for proc in pipe.procs],
        help="\n".join(help_message),
    )
    parsed_args = parser.parse_args(args)
    topic = parsed_args.topic
    if topic is None:
        print(parser.format_help())
    else:
        proc = getattr(processes, topic)
        # Don't link to other processes
        proc.nexts = []
        # Redefine the start process
        Immunopipe.starts = [proc]
        # Clear previous parser actions
        pipen_parser._actions = []
        # Clear previous parser action groups
        pipen_parser._action_groups = []
        # Set the commands
        pipen_parser.set_cli_args(["immunopipe", "--help+"])
        pipe = Immunopipe(plugin_opts={"args_flatten": False})
        pipe.build_proc_relationships()
        pipe.run()
