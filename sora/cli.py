"""Command-line interface for inspecting and managing a SORA installation."""

import argparse
import json
import os
import subprocess
import sys
from importlib import metadata
from typing import Any


PACKAGE_NAME = "sora-astro"
CONFIG_LEVELS = {"basic": 1, "advanced": 2, "develop": 3}
PUBLIC_CONFIG_MODES = ("basic", "advanced")


def get_version() -> str:
    """Return the installed SORA distribution version.

    Returns
    -------
    `str`
        Installed version, or ``"unknown"`` when the package metadata is not
        available.
    """
    try:
        return metadata.version(PACKAGE_NAME)
    except metadata.PackageNotFoundError:
        return "unknown"


def is_editable_install() -> bool:
    """Check whether SORA was installed in editable mode.

    The check uses the ``direct_url.json`` metadata defined for direct URL
    installations. Missing or malformed metadata is ignored.

    Returns
    -------
    `bool`
        `True` if any installed SORA distribution is marked as editable;
        otherwise `False`.
    """
    for distribution in metadata.distributions(name=PACKAGE_NAME):
        direct_url = distribution.read_text("direct_url.json")
        if direct_url is None:
            continue

        try:
            install_info = json.loads(direct_url)
        except (json.JSONDecodeError, TypeError):
            continue

        if not isinstance(install_info, dict):
            continue

        dir_info = install_info.get("dir_info")
        if isinstance(dir_info, dict) and dir_info.get("editable") is True:
            return True

    return False


def _new_config():
    """Create a configuration instance for a CLI operation.

    Returns
    -------
    `sora.config.Config`
        Newly loaded SORA configuration.

    Notes
    -----
    The import is intentionally delayed so commands such as ``sora version``
    do not need to import the configuration subsystem.
    """
    from sora.config import Config

    return Config()


def _questionary():
    """Return the lazily imported interactive-prompt module.

    Returns
    -------
    `module`
        The :mod:`questionary` module used by interactive commands.
    """
    import questionary

    return questionary


def _parse_config_value(raw_value: str, current_value: Any) -> Any:
    """Parse a CLI value according to the current configuration type.

    Parameters
    ----------
    raw_value : `str`
        Text supplied on the command line or through an interactive prompt.
    current_value : `object`
        Existing setting whose type determines how ``raw_value`` is parsed.

    Returns
    -------
    `object`
        Parsed value. Strings are returned verbatim, while container and
        untyped values are decoded as YAML.

    Raises
    ------
    ValueError
        If the text cannot be converted to the expected type or a decoded
        container has the wrong type.
    """
    value = raw_value.strip()

    if isinstance(current_value, bool):
        boolean_values = {
            "true": True,
            "yes": True,
            "y": True,
            "1": True,
            "false": False,
            "no": False,
            "n": False,
            "0": False,
        }
        try:
            return boolean_values[value.lower()]
        except KeyError:
            raise ValueError(
                "expected a boolean: true/false, yes/no, y/n, or 1/0"
            ) from None

    if isinstance(current_value, int):
        try:
            return int(value)
        except ValueError:
            raise ValueError("expected an integer") from None

    if isinstance(current_value, float):
        try:
            return float(value)
        except ValueError:
            raise ValueError("expected a number") from None

    if isinstance(current_value, (list, tuple, dict)):
        # YAML supports structured values while remaining convenient to enter
        # as a single command-line argument.
        import yaml

        try:
            parsed = yaml.safe_load(value)
        except yaml.YAMLError as exc:
            raise ValueError(f"invalid YAML value: {exc}") from None

        if isinstance(current_value, tuple) and isinstance(parsed, list):
            return tuple(parsed)
        if not isinstance(parsed, type(current_value)):
            expected = type(current_value).__name__
            raise ValueError(f"expected a YAML {expected}")
        return parsed

    if current_value is None:
        import yaml

        return yaml.safe_load(value)

    return raw_value


def _split_config_key(dotted_key: str) -> tuple[str, str]:
    """Split and validate a configuration key in ``SECTION.KEY`` notation.

    Parameters
    ----------
    dotted_key : `str`
        Configuration key to split.

    Returns
    -------
    section_name : `str`
        Name before the separator.
    key : `str`
        Field name after the separator.

    Raises
    ------
    ValueError
        If the key does not contain exactly one separator or either component
        is empty.
    """
    try:
        section_name, key = dotted_key.split(".", 1)
    except ValueError:
        raise ValueError("configuration keys must use SECTION.KEY") from None

    if not section_name or not key or "." in key:
        raise ValueError("configuration keys must use SECTION.KEY")

    return section_name, key


def _resolve_config_key(config, dotted_key: str):
    """Resolve a ``SECTION.KEY`` name to its configuration section.

    Parameters
    ----------
    config : `sora.config.Config`
        Configuration containing the registered sections.
    dotted_key : `str`
        Configuration key in ``SECTION.KEY`` notation.

    Returns
    -------
    section_name : `str`
        Name of the resolved section.
    key : `str`
        Field name within the section.
    section : `sora.config.BaseConfigSection`
        Section object containing the field.

    Raises
    ------
    ValueError
        If ``dotted_key`` does not use exactly one separator.
    KeyError
        If the section or field is not registered.
    """
    section_name, key = _split_config_key(dotted_key)

    try:
        section = config.sections[section_name]
    except KeyError:
        raise KeyError(f"unknown configuration section: {section_name}") from None

    if key not in section.FIELDS:
        raise KeyError(f"unknown configuration key: {dotted_key}")

    return section_name, key, section


def _format_config_value(value: Any) -> str:
    """Format a configuration value for terminal output.

    Parameters
    ----------
    value : `object`
        Value to display.

    Returns
    -------
    `str`
        Human-readable scalar text or a YAML representation for containers.
    """
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    if isinstance(value, (dict, list, tuple)):
        import yaml

        return yaml.safe_dump(
            value,
            sort_keys=False,
            default_flow_style=False,
        ).rstrip()
    return str(value)


def _typed_choices(values: list[Any]) -> list[dict[str, Any]]:
    """Build prompt choices that preserve their original value types.

    Parameters
    ----------
    values : `list`
        Values accepted by a configuration prompt.

    Returns
    -------
    `list` of `dict`
        Questionary choices with formatted labels and unmodified values.
    """
    return [
        {"name": _format_config_value(value), "value": value}
        for value in values
    ]


def interactive_config(mode: str = "basic") -> int:
    """Interactively edit configuration and persist all changes at once.

    Parameters
    ----------
    mode : {``"basic"``, ``"advanced"``, ``"develop"``}, optional
        Maximum prompt-detail level to expose.

    Returns
    -------
    `int`
        Zero after saving, discarding, or aborting the operation.

    Raises
    ------
    KeyError
        If ``mode`` is not a registered configuration level.
    """
    level = CONFIG_LEVELS[mode]
    config = _new_config()
    schema = config.get_prompt_schema()
    section_names = list(schema)
    prompts = _questionary()

    print(f"Configuring in '{mode}' mode.")
    section_choice = prompts.select(
        "Which section would you like to edit?",
        choices=["All sections", *section_names, "Exit"],
    ).ask()

    if section_choice in (None, "Exit"):
        print("Configuration aborted.")
        return 0

    sections_to_edit = (
        section_names if section_choice == "All sections" else [section_choice]
    )
    changed_keys: list[str] = []

    for section_name in sections_to_edit:
        section = config.sections[section_name]
        properties = schema[section_name]
        visible_properties = {
            key: prompt
            for key, prompt in properties.items()
            if prompt.get("level", 1) <= level
        }
        if not visible_properties:
            continue

        print(f"── {section_name.upper()} ──")
        for key, prompt in visible_properties.items():
            while True:
                current = getattr(section, key)
                keep_label = f"Keep current ({current})"

                if "choices" in prompt:
                    answer = prompts.select(
                        prompt["question"],
                        choices=[keep_label, *_typed_choices(prompt["choices"])],
                        default=keep_label,
                    ).ask()
                    if answer is None or answer == keep_label:
                        break
                    new_value = answer
                else:
                    answer = prompts.text(
                        f"{prompt['question']} (current: {current})",
                        default="",
                    ).ask()
                    if answer is None or not answer.strip():
                        break
                    try:
                        new_value = _parse_config_value(answer, current)
                    except (TypeError, ValueError) as exc:
                        print(
                            f"Invalid value for {section_name}.{key}: {exc}",
                            file=sys.stderr,
                        )
                        continue

                try:
                    section.set_local(key, new_value)
                except (TypeError, ValueError) as exc:
                    print(
                        f"Invalid value for {section_name}.{key}: {exc}",
                        file=sys.stderr,
                    )
                    continue

                changed_keys.append(f"{section_name}.{key}")
                print(f"→ {section_name}.{key} set to {new_value}")
                break

    if not changed_keys:
        print("No configuration changes were made.")
        return 0

    save_changes = prompts.confirm(
        f"Save {len(changed_keys)} configuration change(s)?",
        default=True,
    ).ask()
    if not save_changes:
        print("Configuration changes discarded.")
        return 0

    config.save_local()
    print(f"Configuration saved to {config.config_path}")
    return 0


def cmd_config_edit(args) -> int:
    """Handle the interactive ``config edit`` command.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed arguments containing the requested configuration ``mode``.

    Returns
    -------
    `int`
        Exit status returned by :func:`interactive_config`.
    """
    return interactive_config(mode=args.mode)


def cmd_config_show(args) -> int:
    """Print the effective configuration or only user overrides.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed arguments containing the ``local`` selection flag.

    Returns
    -------
    `int`
        Zero after printing the configuration.
    """
    config = _new_config()
    document = config.to_local_dict() if args.local else config.to_dict()
    print(_format_config_value(document))
    return 0


def cmd_config_path(args) -> int:
    """Print one or all SORA configuration-related paths.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed arguments whose ``kind`` selects the user configuration,
        bundled defaults, data directory, or all paths.

    Returns
    -------
    `int`
        Zero after printing the requested path information.
    """
    config = _new_config()
    paths = {
        "user": config.config_path,
        "default": config.default_path,
        "data": config.data_path,
    }
    if args.kind == "all":
        for name, path in paths.items():
            print(f"{name}: {path}")
    else:
        print(paths[args.kind])
    return 0


def cmd_config_validate(args) -> int:
    """Load and report the validity of the current configuration.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed command arguments. This command defines no command-specific
        arguments.

    Returns
    -------
    `int`
        Zero when configuration loading and validation succeed.

    Notes
    -----
    Validation occurs while the configuration object is constructed; any
    validation exception is handled by :func:`main`.
    """
    config = _new_config()
    print(f"Configuration is valid: {config.config_path}")
    return 0


def cmd_config_get(args) -> int:
    """Print one effective configuration value.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed arguments containing a ``SECTION.KEY`` name in ``key``.

    Returns
    -------
    `int`
        Zero after printing the value.
    """
    config = _new_config()
    _, key, section = _resolve_config_key(config, args.key)
    print(_format_config_value(getattr(section, key)))
    return 0


def cmd_config_set(args) -> int:
    """Parse and persist one user configuration override.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed arguments containing a ``SECTION.KEY`` name and its new value.

    Returns
    -------
    `int`
        Zero after updating the value.
    """
    config = _new_config()
    section_name, key, section = _resolve_config_key(config, args.key)
    current = getattr(section, key)
    value = _parse_config_value(args.value, current)
    config.update(section_name, key, value)
    print(f"{args.key} set to {_format_config_value(value)}")
    return 0


def cmd_config_reset(args) -> int:
    """Remove one user override, optionally after confirmation.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed arguments containing a ``SECTION.KEY`` name and the ``yes``
        flag used to bypass confirmation.

    Returns
    -------
    `int`
        Zero when the override is removed, absent from a valid configuration,
        or retained by choice; one if recovery cannot restore validity.
    """
    try:
        config = _new_config()
    except (KeyError, TypeError, ValueError) as load_error:
        return _reset_invalid_config(args, load_error)

    section_name, key, section = _resolve_config_key(config, args.key)
    if key not in section.to_dict_local():
        print(f"No user override is set for {args.key}.")
        return 0

    if not args.yes:
        confirmed = _questionary().confirm(
            f"Reset {args.key} to its default value?",
            default=False,
        ).ask()
        if not confirmed:
            print("Reset canceled.")
            return 0

    config.remove_override(section_name, key)
    print(f"{args.key} reset to its default value.")
    return 0


def _reset_invalid_config(args, load_error: Exception) -> int:
    """Remove an override directly when normal configuration loading fails.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed arguments containing a ``SECTION.KEY`` name and confirmation
        flag.
    load_error : `Exception`
        Error raised while applying the current user configuration.

    Returns
    -------
    `int`
        Zero if the reset restores a valid configuration or is canceled; one
        if the requested override is absent or another error remains.
    """
    from sora.config import Config

    section_name, key = _split_config_key(args.key)
    if not args.yes:
        confirmed = _questionary().confirm(
            f"Configuration is invalid. Reset {args.key} directly in the "
            "user configuration file?",
            default=False,
        ).ask()
        if not confirmed:
            print("Reset canceled.")
            return 0

    removed = Config.remove_override_from_file(section_name, key)
    if not removed:
        print(f"No user override is set for {args.key}.")
        print(f"Configuration remains invalid: {load_error}", file=sys.stderr)
        return 1

    try:
        _new_config()
    except (KeyError, OSError, TypeError, ValueError) as exc:
        print(
            f"{args.key} was removed, but configuration remains invalid: {exc}",
            file=sys.stderr,
        )
        return 1

    print(f"{args.key} reset to its default value.")
    return 0


def cmd_version(args) -> int:
    """Print the installed SORA version.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed command arguments. This command defines no command-specific
        arguments.

    Returns
    -------
    `int`
        Zero after printing the version.
    """
    print(get_version())
    return 0


def cmd_info(args) -> int:
    """Print installation, Python, configuration, and data information.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed command arguments. This command defines no command-specific
        arguments.

    Returns
    -------
    `int`
        Zero after printing the diagnostic information.
    """
    config = _new_config()
    print(f"version: {get_version()}")
    print(f"python: {sys.version.split()[0]}")
    print(f"editable: {'yes' if is_editable_install() else 'no'}")
    print(f"config: {config.config_path}")
    print(f"data: {config.data_path}")
    return 0


def cmd_update(args) -> int:
    """Upgrade a non-editable SORA installation with the active Python.

    Parameters
    ----------
    args : `argparse.Namespace`
        Parsed command arguments. This command defines no command-specific
        arguments.

    Returns
    -------
    `int`
        Zero when the update succeeds or is skipped for an editable install;
        otherwise the exit status returned by pip.
    """
    if is_editable_install():
        print("SORA is installed in editable mode; skipping update.")
        return 0

    if os.environ.get("CONDA_PREFIX"):
        print(
            "Warning: SORA is running in a Conda environment; "
            "the update will use pip in that environment.",
            file=sys.stderr,
        )

    print("Updating SORA package via pip…")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", "--upgrade", PACKAGE_NAME]
    )
    if result.returncode == 0:
        print("SORA successfully updated.")
    else:
        print("Update failed (see pip output above).", file=sys.stderr)
    return result.returncode


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level parser and register all SORA subcommands.

    Returns
    -------
    `argparse.ArgumentParser`
        Fully configured parser whose namespaces include a command handler.
    """
    parser = argparse.ArgumentParser(prog="sora", description="SORA CLI")
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {get_version()}",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    config_parser = subparsers.add_parser(
        "config",
        help="Inspect and edit SORA configuration",
    )
    config_parser.set_defaults(handler=cmd_config_edit, mode="basic")
    config_subparsers = config_parser.add_subparsers(dest="config_command")

    edit_parser = config_subparsers.add_parser(
        "edit",
        help="Interactively edit configuration",
    )
    edit_parser.add_argument(
        "--mode",
        choices=PUBLIC_CONFIG_MODES,
        default="basic",
        help="Configuration detail level",
    )
    edit_parser.set_defaults(handler=cmd_config_edit)

    show_parser = config_subparsers.add_parser(
        "show",
        help="Show effective configuration",
    )
    show_parser.add_argument(
        "--local",
        action="store_true",
        help="Show only user overrides",
    )
    show_parser.set_defaults(handler=cmd_config_show)

    path_parser = config_subparsers.add_parser(
        "path",
        help="Show configuration and data paths",
    )
    path_parser.add_argument(
        "kind",
        nargs="?",
        choices=["all", "user", "default", "data"],
        default="all",
    )
    path_parser.set_defaults(handler=cmd_config_path)

    validate_parser = config_subparsers.add_parser(
        "validate",
        help="Validate the current configuration",
    )
    validate_parser.set_defaults(handler=cmd_config_validate)

    get_parser = config_subparsers.add_parser(
        "get",
        help="Read one configuration value",
    )
    get_parser.add_argument("key", metavar="SECTION.KEY")
    get_parser.set_defaults(handler=cmd_config_get)

    set_parser = config_subparsers.add_parser(
        "set",
        help="Set one configuration value",
    )
    set_parser.add_argument("key", metavar="SECTION.KEY")
    set_parser.add_argument("value")
    set_parser.set_defaults(handler=cmd_config_set)

    reset_parser = config_subparsers.add_parser(
        "reset",
        help="Reset one configuration value",
    )
    reset_parser.add_argument("key", metavar="SECTION.KEY")
    reset_parser.add_argument(
        "-y",
        "--yes",
        action="store_true",
        help="Reset without confirmation",
    )
    reset_parser.set_defaults(handler=cmd_config_reset)

    dev_parser = subparsers.add_parser(
        "dev",
        help="Developer utilities",
    )
    dev_subparsers = dev_parser.add_subparsers(
        dest="dev_command",
        required=True,
    )
    dev_config_parser = dev_subparsers.add_parser(
        "config",
        help="Edit all configuration fields",
        description="Interactively edit every registered configuration field.",
    )
    dev_config_parser.set_defaults(handler=cmd_config_edit, mode="develop")

    version_parser = subparsers.add_parser(
        "version",
        help="Show installed SORA version",
    )
    version_parser.set_defaults(handler=cmd_version)

    info_parser = subparsers.add_parser(
        "info",
        help="Show installation and configuration information",
    )
    info_parser.set_defaults(handler=cmd_info)

    update_parser = subparsers.add_parser(
        "update",
        help="Upgrade SORA to the latest version via pip",
    )
    update_parser.set_defaults(handler=cmd_update)

    return parser


def main(argv=None) -> int:
    """Run the SORA command-line interface.

    Parameters
    ----------
    argv : sequence of `str`, optional
        Arguments to parse. If omitted, :data:`sys.argv` is used.

    Returns
    -------
    `int`
        Selected command's exit status, ``130`` for interruption, or ``1``
        for a handled command error.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except KeyboardInterrupt:
        print("Operation canceled.", file=sys.stderr)
        return 130
    except (EOFError, KeyError, OSError, TypeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
