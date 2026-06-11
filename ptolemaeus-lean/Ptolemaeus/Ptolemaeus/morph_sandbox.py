"""Lokaler Stub für Morph.py: führt Code im Prozess aus und fängt stdout.

Für die echte Morph-Cloud-Umgebung siehe https://github.com/morph-labs/morph-python-sdk
(Paket ``morphcloud``) — API weicht ab (kein ``MorphSandbox.execute_code`` in der Doku).
"""

from __future__ import annotations

import asyncio
import io
from contextlib import redirect_stdout
from types import TracebackType
from typing import Any


class MorphSandbox:
    def __init__(self) -> None:
        self._globals: dict[str, Any] = {"__builtins__": __builtins__}

    @classmethod
    async def create(cls) -> MorphSandbox:
        return cls()

    async def __aenter__(self) -> MorphSandbox:
        return self

    async def __aexit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        return None

    async def execute_code(self, code: str) -> dict[str, Any]:
        buf = io.StringIO()

        def run() -> None:
            with redirect_stdout(buf):
                exec(compile(code, "<sandbox>", "exec"), self._globals, self._globals)

        await asyncio.to_thread(run)
        return {"output": buf.getvalue()}
