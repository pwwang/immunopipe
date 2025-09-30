"""Run the entire pipeline on Google Cloud Batch using pipen-cli-gbatch."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

from simpleconf import Config
from yunpath import AnyPath
from xqute.schedulers.gbatch_scheduler import DEFAULT_MOUNTED_ROOT
from pipen.defaults import CONFIG_FILES
from pipen_args.plugin import ArgsPlugin
from pipen_cli_gbatch import (
    CliGbatchDaemon,
    CliGbatchPlugin,
    __version__ as cli_gbatch_version,
    __file__ as cli_gbatch_file,
)
from .pipeline import Immunopipe, parser

sys.excepthook = sys.__excepthook__


class ImmunopipeGbatchDaemon(CliGbatchDaemon):

    def _run_version(self):
        """Print version information for pipen-cli-gbatch and pipen."""
        from .version import __version__
        from pipen import __version__ as pipen_version

        print(f"Immunopipe version: v{__version__}")
        print(f"pipen-cli-gbatch version: v{cli_gbatch_version}")
        print(f"pipen version: v{pipen_version}")

    def _infer_name(self):
        super()._infer_name()
        immunopipe_name = self._get_arg_from_command("name")
        self.config["workdir"] = f"{self.config['workdir']}/{immunopipe_name}"

        # Copy configuration file over
        cf_at = [cmd.startswith("@") for cmd in self.command]
        if any(cf_at):
            cf_index = cf_at.index(True)
            cf_path = AnyPath(self.command[cf_index][1:])
            cf_dest = AnyPath(self.config["workdir"]).joinpath(
                self.config["name"],
                cf_path.name,
            )
            cf_dest.write_bytes(cf_path.read_bytes())
            self.command[cf_index] = (
                f"@{DEFAULT_MOUNTED_ROOT}/xqute_workdir/{cf_path.name}"
            )


async def main(argv):
    """The main function to run the pipeline on Google Cloud Batch."""
    act_group = parser._extra_parser.add_mutually_exclusive_group(required=False)
    act_group.title = "Actions (extra options)"
    act_group.description = (
        "Specify one of the following actions. "
        "If none is specified, the pipeline will be run and waited for completion."
    )

    cli_gbatch_arg_file = Path(cli_gbatch_file).parent / "daemon_args.toml"
    cli_gbatch_args = Config.load_one(cli_gbatch_arg_file, loader="toml")

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
                arg_args[i] = aa.replace("--", "--cli-gbatch.")
            elif aa.startswith("-"):
                arg_args[i] = aa.replace("-", "--cli-gbatch.")
            else:
                arg_args[i] = f"cli-gbatch.{aa}"

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
    for key, value in vars(cli_gbatch_config.cli_gbatch).items():
        setattr(cli_gbatch_config, key, value)
        delattr(cli_gbatch_config, f"cli_gbatch.{key}")
    del cli_gbatch_config.cli_gbatch

    def is_valid(val: Any) -> bool:
        """Check if a value is valid (not None, not empty string, not empty list).
        """
        if val is None:
            return False
        if isinstance(val, bool):
            return True
        return bool(val)

    defaults = CliGbatchPlugin._get_defaults_from_config(
        CONFIG_FILES,
        cli_gbatch_config.profile,
    )
    # update parsed with the defaults
    for key, val in defaults.items():
        if (
            key == "mount"
            and val
            and getattr(cli_gbatch_config, key, None)
        ):
            if not isinstance(val, (tuple, list)):
                val = [val]
            val = list(val)

            kp_mount = getattr(cli_gbatch_config, key)
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


    mount_as_cwd = getattr(cli_gbatch_config, "mount_as_cwd", None)
    cwd = getattr(cli_gbatch_config, "cwd", None)
    delattr(cli_gbatch_config, "mount_as_cwd")
    if mount_as_cwd and cwd:
        print(
            "\033[1;4mError\033[0m: --mount-as-cwd and --cwd "
            "cannot be used together.\n"
        )
        sys.exit(1)

    mount = getattr(cli_gbatch_config, "mount", None) or []
    if mount_as_cwd:
        mount.append(f"{mount_as_cwd}:/mnt/disks/.cwd")
        setattr(cli_gbatch_config, "mount", mount)
        setattr(cli_gbatch_config, "cwd", "/mnt/disks/.cwd")

    cli_gbatch_config.name = ".ImmunopipeCliGbatch"
    cli_gbatch_config.plain = False
    cli_gbatch_config.workdir = None  # will infer from command
    cli_gbatch_config.jobname_prefix = "immunopipe-cli-gbatch"
    # cli_gbatch_config.cwd = None
    cli_gbatch_config.entrypoint = "/usr/local/bin/_entrypoint.sh"
    cli_gbatch_config.commands = ["{lang}", "{script}"]

    pipe = Immunopipe()
    # Let on_init() hook handle argument parsing
    if len(argv) == 0:
        pipe.run()
    elif not cli_gbatch_config.version and not cli_gbatch_config.view_logs:
        await ArgsPlugin.on_init.impl(pipe)

    if not cli_gbatch_config.project:
        print("\033[1;4mError\033[0m: --cli-gbatch.project is required.\n")
        sys.exit(1)
    if not cli_gbatch_config.location:
        print("\033[1;4mError\033[0m: --cli-gbatch.location is required.\n")
        sys.exit(1)

    command = ["immunopipe", *parser._cli_args]
    daemon = ImmunopipeGbatchDaemon(cli_gbatch_config, command)
    await daemon.run()
