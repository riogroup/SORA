from __future__ import annotations

import hashlib
import json
import os
import tempfile
import time
from pathlib import Path
from typing import Any, Mapping
from urllib.parse import urlparse

import requests
from tqdm import tqdm

DEFAULT_TIMEOUT = (10, 120)


def validate_http_url(url: str) -> None:
    """Validate a URL used as a remote SORA data source."""
    if not isinstance(url, str) or not url.strip():
        raise TypeError('Download URL must be a non-empty string')
    parsed = urlparse(url)
    if parsed.scheme not in {'http', 'https'} or not parsed.netloc:
        raise ValueError(f'Invalid HTTP or HTTPS URL: {url!r}')


def filename_from_url(url: str) -> str:
    """Return a safe filename extracted from an HTTP URL."""
    validate_http_url(url)
    filename = Path(urlparse(url).path).name
    if not filename or filename in {'.', '..'}:
        raise ValueError(f'URL does not contain a valid filename: {url!r}')
    return filename


def fetch_text(
        url: str,
        *,
        session: requests.Session | None = None,
        timeout=DEFAULT_TIMEOUT,
        retries: int = 3,
) -> str:
    """Fetch text with bounded retries and request timeouts."""
    validate_http_url(url)
    if retries < 1:
        raise ValueError('retries must be at least 1')

    client = session or requests.Session()
    owns_session = session is None
    try:
        for attempt in range(retries):
            try:
                with client.get(url, timeout=timeout) as response:
                    response.raise_for_status()
                    return response.text.strip()
            except requests.RequestException:
                if attempt == retries - 1:
                    raise
                time.sleep(2 ** attempt)
    finally:
        if owns_session:
            client.close()

    raise RuntimeError(f'Failed to fetch {url}')


def file_matches(
        path: Path | str,
        *,
        expected_size: int | None = None,
        expected_sha256: str | None = None,
) -> bool:
    """Return whether an existing file satisfies configured integrity data."""
    path = Path(path)
    if not path.is_file():
        return False
    if path.stat().st_size == 0:
        return False
    if expected_size is not None and path.stat().st_size != expected_size:
        return False
    if expected_sha256 is not None:
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(chunk)
        if digest.hexdigest().lower() != expected_sha256.lower():
            return False
    return True


def download_file(
        url: str,
        destination: Path | str,
        *,
        session: requests.Session | None = None,
        timeout=DEFAULT_TIMEOUT,
        retries: int = 3,
        expected_size: int | None = None,
        expected_sha256: str | None = None,
        description: str | None = None,
        progress: bool = True,
) -> Path:
    """Download a file atomically and validate available integrity metadata."""
    validate_http_url(url)
    if retries < 1:
        raise ValueError('retries must be at least 1')
    if expected_size is not None and expected_size < 0:
        raise ValueError('expected_size must be non-negative')
    if expected_sha256 is not None:
        normalized_hash = expected_sha256.lower()
        if len(normalized_hash) != 64 or any(
                char not in '0123456789abcdef' for char in normalized_hash
        ):
            raise ValueError('expected_sha256 must contain 64 hexadecimal characters')

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    client = session or requests.Session()
    owns_session = session is None

    try:
        for attempt in range(retries):
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=f'.{destination.name}.',
                suffix='.part',
                dir=destination.parent,
            )
            temporary_path = Path(temporary_name)
            try:
                digest = hashlib.sha256()
                bytes_written = 0
                with client.get(url, stream=True, timeout=timeout) as response:
                    response.raise_for_status()
                    header_size = response.headers.get('content-length')
                    response_size = int(header_size) if header_size else None
                    stream = os.fdopen(descriptor, 'wb')
                    descriptor = -1
                    with stream, tqdm(
                            desc=description or destination.name,
                            total=response_size,
                            unit='B',
                            unit_scale=True,
                            unit_divisor=1024,
                            disable=not progress,
                    ) as bar:
                        for chunk in response.iter_content(chunk_size=1024 * 1024):
                            if not chunk:
                                continue
                            stream.write(chunk)
                            digest.update(chunk)
                            bytes_written += len(chunk)
                            bar.update(len(chunk))
                        stream.flush()
                        os.fsync(stream.fileno())

                if response_size is not None and bytes_written != response_size:
                    raise ValueError(
                        f'Incomplete download for {url}: expected '
                        f'{response_size} bytes, received {bytes_written}'
                    )
                if expected_size is not None and bytes_written != expected_size:
                    raise ValueError(
                        f'Unexpected size for {url}: expected '
                        f'{expected_size} bytes, received {bytes_written}'
                    )
                if (
                        expected_sha256 is not None
                        and digest.hexdigest().lower() != expected_sha256.lower()
                ):
                    raise ValueError(f'Checksum mismatch for {url}')

                os.replace(temporary_path, destination)
                return destination
            except (OSError, ValueError, requests.RequestException):
                if descriptor >= 0:
                    os.close(descriptor)
                temporary_path.unlink(missing_ok=True)
                if attempt == retries - 1:
                    raise
                time.sleep(2 ** attempt)
    finally:
        if owns_session:
            client.close()

    raise RuntimeError(f'Failed to download {url}')


def write_json_atomic(path: Path | str, data: Mapping[str, Any]) -> None:
    """Write JSON without exposing a partially written destination file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f'.{path.name}.',
        suffix='.tmp',
        dir=path.parent,
        text=True,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as stream:
            json.dump(data, stream, indent=2)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_path, path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


def write_bytes_atomic(path: Path | str, data: bytes) -> None:
    """Write bytes without exposing a partially written destination file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f'.{path.name}.',
        suffix='.tmp',
        dir=path.parent,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, 'wb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_path, path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise
