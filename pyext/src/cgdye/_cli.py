"""The slice of click's decorator API that cgdye's command lines use, on argparse.

IMP.bff carries no dependency beyond what IMP itself brings, and click was the
last external left in cgdye — 110 ``option`` decorators across five modules.
Rewriting those by hand into argparse would have been a 119-declaration diff
with a transcription error in it somewhere; implementing the surface they
actually use is smaller, mechanical, and verifiable by running every ``--help``.

The surface is deliberately closed, and measured rather than guessed:
``option``, ``argument``, ``command``, ``group``, ``echo``, ``Path``,
``Choice``, ``ClickException``, ``UsageError``, with the keywords
``default``, ``help``, ``show_default``, ``required``, ``multiple``,
``is_flag`` and ``type``. Anything else raises at decoration time instead of
being quietly ignored — a CLI that silently drops an option is worse than one
that will not import.

Decorators run bottom-up, so parameters accumulate in reverse and are reversed
once when the command is built; that is why the order in ``--help`` matches the
source.
"""

from __future__ import annotations

import argparse
import os
import sys

__all__ = [
    "option", "argument", "command", "group", "echo",
    "Path", "Choice", "ClickException", "UsageError",
]

_SUPPORTED = {
    "default", "help", "show_default", "required", "multiple",
    "is_flag", "type", "nargs", "metavar",
}


class ClickException(Exception):
    """An error worth showing the user without a traceback."""

    exit_code = 1

    def show(self) -> None:
        print(f"Error: {self}", file=sys.stderr)


class UsageError(ClickException):
    """The command line itself was wrong."""

    exit_code = 2


class Choice:
    """A fixed set of allowed values."""

    def __init__(self, choices, case_sensitive: bool = True):
        self.choices = list(choices)
        self.case_sensitive = case_sensitive

    def __call__(self, value):
        candidates = self.choices
        if not self.case_sensitive:
            lowered = {c.lower(): c for c in self.choices}
            if str(value).lower() in lowered:
                return lowered[str(value).lower()]
        if value not in candidates:
            raise argparse.ArgumentTypeError(
                f"{value!r} is not one of {', '.join(map(str, candidates))}")
        return value


class Path:
    """A filesystem path, optionally checked for existence and kind."""

    def __init__(self, exists: bool = False, file_okay: bool = True,
                 dir_okay: bool = True, writable: bool = False,
                 readable: bool = True, resolve_path: bool = False,
                 path_type=None):
        self.exists = exists
        self.file_okay = file_okay
        self.dir_okay = dir_okay
        self.resolve_path = resolve_path
        self.path_type = path_type

    def __call__(self, value):
        text = str(value)
        if self.exists and not os.path.exists(text):
            raise argparse.ArgumentTypeError(f"path does not exist: {text}")
        if os.path.exists(text):
            if not self.file_okay and os.path.isfile(text):
                raise argparse.ArgumentTypeError(f"expected a directory: {text}")
            if not self.dir_okay and os.path.isdir(text):
                raise argparse.ArgumentTypeError(f"expected a file: {text}")
        if self.resolve_path:
            text = os.path.realpath(text)
        return self.path_type(text) if self.path_type else text


def echo(message="", err: bool = False, nl: bool = True) -> None:
    """Print, with click's signature."""
    print(message, file=sys.stderr if err else sys.stdout, end="\n" if nl else "")


def _check(attrs: dict) -> None:
    unknown = set(attrs) - _SUPPORTED
    if unknown:
        raise TypeError(
            f"cgdye._cli does not implement click option(s) {sorted(unknown)}; "
            f"add them here rather than reaching for click")


def _pending(func) -> list:
    if not hasattr(func, "__cli_params__"):
        func.__cli_params__ = []
    return func.__cli_params__


def option(*decls, **attrs):
    """Record an optional argument, in click's spelling."""
    _check(attrs)

    def decorator(func):
        _pending(func).append(("option", decls, attrs))
        return func
    return decorator


def argument(*decls, **attrs):
    """Record a positional argument, in click's spelling."""
    _check(attrs)

    def decorator(func):
        _pending(func).append(("argument", decls, attrs))
        return func
    return decorator


def _dest(decls, kind):
    """Longest long-form flag becomes the destination, as click does."""
    if kind == "argument":
        return decls[0].lstrip("-").replace("-", "_")
    longs = [d for d in decls if d.startswith("--")]
    name = (longs[0] if longs else decls[0]).lstrip("-")
    return name.replace("-", "_")


class Command:
    """A callable command built from the recorded parameters."""

    def __init__(self, func, name=None, help=None):
        self.callback = func
        self.name = name or func.__name__.replace("_", "-")
        self.help = help or (func.__doc__ or "").strip()
        self.params = list(reversed(getattr(func, "__cli_params__", [])))

    def build_parser(self, parser=None):
        parser = parser or argparse.ArgumentParser(
            prog=self.name, description=self.help)
        for kind, decls, attrs in self.params:
            attrs = dict(attrs)
            dest = _dest(decls, kind)
            help_text = attrs.pop("help", None)
            show_default = attrs.pop("show_default", False)
            default = attrs.pop("default", None)
            if show_default and default is not None and help_text:
                help_text = f"{help_text}  [default: {default}]"
            kw = {"help": help_text}
            if attrs.pop("is_flag", False):
                kw["action"] = "store_true"
                kw["default"] = bool(default)
            else:
                if attrs.pop("multiple", False):
                    kw["action"] = "append"
                    kw["default"] = list(default) if default else []
                else:
                    kw["default"] = default
                if "type" in attrs and attrs["type"] is not None:
                    kw["type"] = attrs.pop("type")
                for passthrough in ("nargs", "metavar"):
                    if passthrough in attrs:
                        kw[passthrough] = attrs.pop(passthrough)
            required = attrs.pop("required", False)
            if kind == "argument":
                if not required and kw.get("default") is not None:
                    kw["nargs"] = kw.get("nargs", "?")
                kw.pop("help", None) if kw.get("help") is None else None
                parser.add_argument(dest, **{k: v for k, v in kw.items() if v is not None or k == "default"})
            else:
                flags = [d for d in decls if d.startswith("-")] or [f"--{dest}"]
                kw["required"] = required
                kw["dest"] = dest
                parser.add_argument(*flags, **kw)
        return parser

    def main(self, argv=None, standalone_mode: bool = True):
        parser = self.build_parser()
        args = parser.parse_args(argv)
        try:
            return self.callback(**vars(args))
        except ClickException as exc:
            if not standalone_mode:
                raise
            exc.show()
            sys.exit(exc.exit_code)

    # click commands stay callable as plain functions
    def __call__(self, *args, **kwargs):
        if not args and not kwargs:
            return self.main()
        return self.callback(*args, **kwargs)


class Group(Command):
    """A command with subcommands."""

    def __init__(self, func=None, name=None, help=None):
        if func is None:
            func = lambda **_: None                      # noqa: E731
            func.__name__ = name or "cli"
            func.__doc__ = help or ""
        super().__init__(func, name=name, help=help)
        self.commands: dict[str, Command] = {}

    def add_command(self, cmd: Command, name: str | None = None) -> None:
        self.commands[name or cmd.name] = cmd

    def command(self, name=None, **attrs):
        def decorator(func):
            cmd = Command(func, name=name, help=attrs.get("help"))
            self.add_command(cmd)
            return cmd
        return decorator

    def main(self, argv=None, standalone_mode: bool = True):
        parser = argparse.ArgumentParser(prog=self.name, description=self.help)
        subparsers = parser.add_subparsers(dest="_command", required=True)
        for name, cmd in self.commands.items():
            sub = subparsers.add_parser(name, description=cmd.help, help=cmd.help)
            cmd.build_parser(sub)
        args = vars(parser.parse_args(argv))
        chosen = self.commands[args.pop("_command")]
        try:
            return chosen.callback(**args)
        except ClickException as exc:
            if not standalone_mode:
                raise
            exc.show()
            sys.exit(exc.exit_code)


def command(name=None, **attrs):
    """Turn a function into a command, in click's spelling."""
    def decorator(func):
        return Command(func, name=name, help=attrs.get("help"))
    if callable(name):
        func, name = name, None
        return Command(func)
    return decorator


def group(name=None, **attrs):
    """Turn a function into a command group, in click's spelling."""
    def decorator(func):
        return Group(func, name=name, help=attrs.get("help"))
    if callable(name):
        return Group(name)
    return decorator
