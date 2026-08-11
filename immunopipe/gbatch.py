"""Run the entire pipeline on Google Cloud Batch using pipen-cli-gbatch."""

from __future__ import annotations

import sys
from typing import Any

from simpleconf import Config
from panpath import PanPath
from xqute.utils import logger
from pipen.defaults import CONFIG_FILES
from pipen_args.plugin import ArgsPlugin
from pipen_cli_gbatch import (
    CliGbatchDaemonPipeline,
    CliGbatchPlugin,
    __version__ as cli_gbatch_version,
    __file__ as cli_gbatch_file,
)
from .pipeline import Immunopipe, parser  # type: ignore

sys.excepthook = sys.__excepthook__


class ImmunopipeGbatchDaemon(CliGbatchDaemonPipeline):

    def _run_version(self):
        """Print version information for pipen-cli-gbatch and pipen."""
        from .version import __version__
        from pipen import __version__ as pipen_version

        print(f"Immunopipe version: v{__version__}")
        print(f"pipen version: v{pipen_version}")
        print(f"pipen-cli-gbatch version: v{cli_gbatch_version}")

    def _show_versions(self):
        """Log the version information for debugging purposes."""
        from .version import __version__

        logger.info(f"Immunopipe version: v{__version__}")
        super()._show_versions()

    async def handle_workdir(self):
        await super().handle_workdir()

        if "workdir" in self._command_args:
            del self._command_args["workdir"]

        mounted_workdir = await self._get_arg_from_command("workdir")

        # Copy configuration file over
        cf_at = [cmd.startswith("@") for cmd in self.command]
        if any(cf_at):
            command_name = await self.command_name()
            command_workdir = await self.command_workdir()
            cf_index = cf_at.index(True)
            cf_path = PanPath(self.command[cf_index][1:])
            cf_dest = PanPath(command_workdir).joinpath(
                self.daemon_name,
                cf_path.name,
            )
            await cf_path.a_copy(cf_dest)

            self.command[cf_index] = (
                f"@{mounted_workdir}/{command_name}/{self.daemon_name}/{cf_path.name}"
            )


async def main(argv):
    """The main function to run the pipeline on Google Cloud Batch."""
    from .version import __version__

    act_group = parser._extra_parser.add_mutually_exclusive_group(required=False)
    act_group.title = "Actions (extra options)"
    act_group.description = (
        "Specify one of the following actions. "
        "If none is specified, the pipeline will be run and waited for completion."
    )

    cli_gbatch_arg_file = PanPath(cli_gbatch_file).parent / "daemon_args.toml"
    cli_gbatch_args = await Config.a_load_one(cli_gbatch_arg_file, loader="toml")

    action_options = cli_gbatch_args.mutually_exclusive_groups[0].arguments
    for opt in action_options:
        opt_args = opt.pop("flags")
        parser.add_extra_argument(*opt_args, **opt, group="Actions")

    cli_gbatch_arguments = cli_gbatch_args.arguments + [
        arg
        for group in cli_gbatch_args.groups
        for arg in group.arguments
        if arg.flags != ["command"]
    ]
    for arg in cli_gbatch_arguments:
        arg_args = arg.pop("flags")
        if (
            "--name" in arg_args
            or "--plain" in arg_args
            or "--workdir" in arg_args
            or "--jobname-prefix" in arg_args
            # or "--cwd" in arg_args
            or "--entrypoint" in arg_args
            or "--commands" in arg_args
        ):
            continue

        for i, aa in enumerate(arg_args):
            if aa.startswith("--"):
                arg_args[i] = aa.replace("--", "--gbatch.")
            elif aa.startswith("-"):
                arg_args[i] = aa.replace("-", "--gbatch.")
            else:
                arg_args[i] = f"gbatch.{aa}"

        parser.add_extra_argument(
            *arg_args,
            **arg,
            group="Options for pipen-cli-gbatch",
        )

    parser.prog = "immunopipe gbatch"
    parser.description = (
        "Run the entire Immunopipe pipeline on "
        "Google Cloud Batch using pipen-cli-gbatch."
    )
    parser.set_cli_args(argv)

    cli_gbatch_config = parser.parse_extra_args(fromfile_parse=True, fromfile_keep=True)
    for key, value in vars(cli_gbatch_config.gbatch).items():
        setattr(cli_gbatch_config, key, value)
        delattr(cli_gbatch_config, f"gbatch.{key}")
    del cli_gbatch_config.gbatch

    def is_valid(val: Any) -> bool:
        """Check if a value is valid (not None, not empty string, not empty list).
        """
        if val is None:
            return False
        if isinstance(val, bool):
            return True
        return bool(val)

    defaults = await CliGbatchPlugin._get_defaults_from_config(
        CONFIG_FILES,
        cli_gbatch_config.profile,
    )
    default_scheduler_opts = defaults.pop("scheduler_opts", {})
    # update parsed with the defaults
    for key, val in default_scheduler_opts.items():
        if (
            key == "mount"
            and val
            and getattr(cli_gbatch_config, key, None)
        ):
            if not isinstance(val, (tuple, list)):
                val = [val]
            val = list(val)

            kp_mount = getattr(cli_gbatch_config, key)
            if not isinstance(kp_mount, (tuple, list)):
                val.append(kp_mount)
            else:
                val.extend(kp_mount)
            setattr(cli_gbatch_config, key, val)
            continue

        if (
            key == "command"
            or val is None
            or is_valid(getattr(cli_gbatch_config, key, None))
        ):
            continue

        setattr(cli_gbatch_config, key, val)

    cli_gbatch_config.name = ".ImmunopipeGbatch"
    cli_gbatch_config.plain = False
    cli_gbatch_config.workdir = None  # will infer from command
    cli_gbatch_config.jobname_prefix = "immunopipe-gbatch"
    # cli_gbatch_config.cwd = None
    cli_gbatch_config.entrypoint = "/usr/local/bin/_entrypoint.sh"
    cli_gbatch_config.commands = ["{lang}", "{script}"]

    pipe = Immunopipe()
    # Let on_init() hook handle argument parsing, e.g. print help message
    if len(argv) == 0:
        # Print a simple usage message
        print(
            "Usage: immunopipe gbatch [options to run pipeline on Google Cloud Batch]"
        )
        print("Try 'immunopipe gbatch --help' for more information.\n")
        sys.exit(0)
    elif not cli_gbatch_config.version and not cli_gbatch_config.view_logs:
        # Set parser._cli_args
        await ArgsPlugin.on_init.impl(pipe)

    if not cli_gbatch_config.project:
        print("\033[1;4mError\033[0m: --gbatch.project is required.\n")
        sys.exit(1)
    if not cli_gbatch_config.location:
        print("\033[1;4mError\033[0m: --gbatch.location is required.\n")
        sys.exit(1)

    setattr(cli_gbatch_config, "_other_opts", defaults)
    command = ["immunopipe", *parser._cli_args]
    daemon = ImmunopipeGbatchDaemon(cli_gbatch_config, command)
    daemon.envs["IMMUNOPIPE_HOST_VERSION"] = __version__
    await daemon.run()
