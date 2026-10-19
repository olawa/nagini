"""Run a command-line script from a Snakemake ``script:`` rule.

The BCR-ABL1 and Dorado scripts are argparse command-line tools, which is how
their tests and anyone at a prompt use them.  Each also has a
``snakemake_argv(snakemake)`` that turns the rule's input, output and params
into the same argument list, so a rule can use ``script:`` and the two ways of
running a script cannot drift apart: they parse the same arguments.

``run`` does what the shell redirect ``> {log} 2>&1`` did.  Standard library
only.
"""

from __future__ import annotations

import contextlib
import sys
import traceback
from typing import Callable, Optional, Sequence


def run(main: Callable[[Sequence[str]], Optional[int]], argv: Sequence[str], log=None) -> int:
    """Call ``main(argv)`` with stdout and stderr going to the rule's log.

    Returns the exit status.  A ``sys.exit("message")`` inside the script ends
    up in the log, not only on the console, and counts as failure.
    """

    log_path = log[0] if log and not isinstance(log, str) else log
    with contextlib.ExitStack() as stack:
        stream = stack.enter_context(open(log_path, "w")) if log_path else sys.stdout
        if log_path:
            stack.enter_context(contextlib.redirect_stdout(stream))
            stack.enter_context(contextlib.redirect_stderr(stream))
        try:
            main([str(arg) for arg in argv])
            status = 0
        except SystemExit as stop:
            if stop.code is None:
                status = 0
            elif isinstance(stop.code, int):
                status = stop.code
            else:
                print(stop.code)
                status = 1
        except Exception:  # noqa: BLE001 - the traceback belongs in the log
            traceback.print_exc()
            status = 1
        stream.flush()
    return status
