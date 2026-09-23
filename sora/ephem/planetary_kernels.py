from __future__ import annotations

import json
from datetime import datetime, timedelta
from importlib import resources
from pathlib import Path
from typing import Any

import requests

from sora.config.core import Config, get_config, resolve_storage_path
from sora.config.download import (
    download_file,
    fetch_text,
    file_matches,
    filename_from_url,
    validate_http_url,
    write_json_atomic,
)

__all__ = ['PlanetaryKernelDB']


class PlanetaryKernelDB:
    _local_suffixes = frozenset({'.bsp'})

    def __init__(
            self,
            config: Config | None = None,
            *,
            session: requests.Session | None = None,
            refresh: bool = True,
    ):
        self.config = config or get_config()
        self._session = session
        self.metadata_url = self.config.ephem.planet_kernels_url
        self.update_age_days = self.config.ephem.update_age_days
        self.metadata_path = resolve_storage_path(
            self.config.data_path,
            self.config.ephem.planet_kernels,
        )
        self.download_dir = resolve_storage_path(
            self.config.data_path,
            self.config.ephem.data_dir,
        )
        self.state_path = self.metadata_path.with_name(
            f'{self.metadata_path.stem}_state.json'
        )
        self.metadata_path.parent.mkdir(parents=True, exist_ok=True)
        self.download_dir.mkdir(parents=True, exist_ok=True)
        self._refresh_state = self._load_refresh_state()
        self.kernels = self._load_metadata()
        if refresh and self.should_refresh_metadata():
            self.refresh_metadata()
        self.sync_local_files()

    @staticmethod
    def _validate_file(name: str, entry: Any) -> dict[str, Any]:
        if not isinstance(entry, dict):
            raise TypeError(f'Kernel file metadata for {name!r} must be an object')
        url = entry.get('url')
        if url is not None and not isinstance(url, str):
            raise TypeError(f'Kernel URL for {name!r} must be a string')
        if url is not None:
            validate_http_url(url)

        local_only = entry.get('local_only', url is None)
        if not isinstance(local_only, bool):
            raise TypeError(f'local_only for kernel {name!r} must be a boolean')
        if url is None and not local_only:
            raise ValueError(
                f'Kernel metadata for {name!r} requires a URL or local_only=true'
            )

        filename = entry.get('filename')
        if filename is None and url is not None:
            filename = filename_from_url(url)
        if not isinstance(filename, str):
            raise TypeError(f'Kernel filename for {name!r} must be a string')
        if Path(filename).name != filename or filename in {'.', '..'}:
            raise ValueError(f'Invalid filename for kernel {name!r}: {filename!r}')

        validated = {'filename': filename}
        if url is not None:
            validated['url'] = url
        if local_only:
            validated['local_only'] = True
        sha256 = entry.get('sha256')
        if sha256 is not None:
            if not isinstance(sha256, str):
                raise TypeError(f'SHA-256 for kernel {name!r} must be a string')
            normalized_hash = sha256.lower()
            if len(normalized_hash) != 64 or any(
                    char not in '0123456789abcdef' for char in normalized_hash
            ):
                raise ValueError(f'Invalid SHA-256 for kernel {name!r}')
            validated['sha256'] = normalized_hash
        size = entry.get('size')
        if size is not None:
            if not isinstance(size, int) or isinstance(size, bool) or size < 0:
                raise ValueError(f'Invalid size for kernel {name!r}')
            validated['size'] = size
        return validated

    @classmethod
    def _validate_entry(cls, name: str, entry: Any) -> dict[str, Any]:
        if not isinstance(name, str) or not name.strip():
            raise ValueError('Kernel names must be non-empty strings')
        if not isinstance(entry, dict):
            raise TypeError(f'Kernel metadata for {name!r} must be an object')

        files = entry.get('files')
        if files is None:
            return cls._validate_file(name, entry)
        if not isinstance(files, list) or not files:
            raise ValueError(
                f'Kernel metadata for {name!r} must contain a non-empty files list'
            )
        conflicting = {'url', 'filename', 'sha256', 'size', 'local_only'}.intersection(
            entry
        )
        if conflicting:
            raise ValueError(
                f'Kernel metadata for {name!r} cannot combine files with '
                f'{sorted(conflicting)}'
            )

        validated_files = [
            cls._validate_file(f'{name}[{index}]', file_entry)
            for index, file_entry in enumerate(files)
        ]
        filenames = [file_entry['filename'].casefold() for file_entry in validated_files]
        if len(filenames) != len(set(filenames)):
            raise ValueError(f'Kernel metadata for {name!r} has duplicate filenames')
        return {'files': validated_files}

    @staticmethod
    def _entry_files(entry: dict[str, Any]) -> list[dict[str, Any]]:
        files = entry.get('files')
        return files if isinstance(files, list) else [entry]

    @classmethod
    def _validate_document(
            cls,
            document: Any,
            *,
            source: str,
    ) -> dict[str, dict[str, Any]]:
        if not isinstance(document, dict):
            raise TypeError(f'{source} must contain a JSON object')

        kernels: dict[str, dict[str, Any]] = {}
        for name, entry in document.items():
            validated = cls._validate_entry(name, entry)
            normalized_name = name.casefold()
            if normalized_name in kernels:
                raise ValueError(
                    f'Duplicate planetary kernel name after normalization '
                    f'in {source}: {name!r}'
                )
            kernels[normalized_name] = validated
        return kernels

    def _load_metadata(self) -> dict[str, dict[str, Any]]:
        item = resources.files('sora.data').joinpath('planet_kernels.json')
        packaged_document = json.loads(item.read_text(encoding='utf-8'))
        packaged_kernels = self._validate_document(
            packaged_document,
            source='Packaged planetary kernel metadata',
        )
        metadata_exists = self.metadata_path.exists()
        if metadata_exists:
            try:
                with self.metadata_path.open('r', encoding='utf-8') as stream:
                    document = json.load(stream)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f'Invalid planetary kernel metadata in '
                    f'{self.metadata_path}: {exc}'
                ) from exc
        else:
            document = {}

        local_kernels = self._validate_document(
            document,
            source=f'Planetary kernel metadata in {self.metadata_path}',
        )
        kernels = {**packaged_kernels, **local_kernels}
        if not metadata_exists or kernels != local_kernels:
            write_json_atomic(self.metadata_path, kernels)
        return kernels

    def _load_refresh_state(self) -> dict[str, Any]:
        if not self.state_path.is_file():
            return {}
        try:
            with self.state_path.open('r', encoding='utf-8') as stream:
                state = json.load(stream)
        except json.JSONDecodeError:
            return {}
        return state if isinstance(state, dict) else {}

    def should_refresh_metadata(self) -> bool:
        """Return whether the remote catalogue should be checked."""
        if self._refresh_state.get('source') != self.metadata_url:
            return True
        last_update = self._refresh_state.get('last_update')
        if not isinstance(last_update, str):
            return True
        try:
            last_update_date = datetime.fromisoformat(last_update)
        except ValueError:
            return True
        return datetime.now() - last_update_date >= timedelta(
            days=self.update_age_days
        )

    def refresh_metadata(self) -> bool:
        """Refresh standard kernels from the remote catalogue."""
        try:
            document = json.loads(
                fetch_text(self.metadata_url, session=self._session)
            )
            remote_kernels = self._validate_document(
                document,
                source=f'Remote planetary kernel metadata at {self.metadata_url}',
            )
        except (
                json.JSONDecodeError,
                requests.RequestException,
                TypeError,
                ValueError,
        ):
            return False

        merged_kernels = {**self.kernels, **remote_kernels}
        changed = merged_kernels != self.kernels
        if changed:
            self.kernels = merged_kernels
            self._save_metadata()

        self._refresh_state = {
            'source': self.metadata_url,
            'last_update': datetime.now().isoformat(),
        }
        write_json_atomic(self.state_path, self._refresh_state)
        return changed

    def _save_metadata(self) -> None:
        write_json_atomic(self.metadata_path, self.kernels)

    def _local_files(self) -> dict[str, Path]:
        return {
            path.name.casefold(): path
            for path in sorted(
                self.download_dir.iterdir(),
                key=lambda candidate: candidate.name.casefold(),
            )
            if path.suffix.casefold() in self._local_suffixes
               and file_matches(path)
        }

    def sync_local_files(self) -> int:
        """Synchronize manually copied BSP files with the kernel metadata."""
        self.download_dir.mkdir(parents=True, exist_ok=True)
        local_files = self._local_files()
        changes = 0

        known_filenames: set[str] = set()
        for name, info in list(self.kernels.items()):
            retained_files = []
            for file_info in self._entry_files(info):
                normalized_filename = file_info['filename'].casefold()
                local_path = local_files.get(normalized_filename)
                if local_path is None and file_info.get('local_only'):
                    changes += 1
                    continue
                retained_files.append(file_info)
                if local_path is None:
                    continue
                known_filenames.add(normalized_filename)
                if file_info['filename'] != local_path.name:
                    file_info['filename'] = local_path.name
                    changes += 1

            if not retained_files:
                del self.kernels[name]
            elif 'files' in info and len(retained_files) != len(info['files']):
                info['files'] = retained_files

        for normalized_filename, path in local_files.items():
            if normalized_filename in known_filenames:
                continue
            name = path.stem.casefold()
            if name in self.kernels:
                continue
            self.kernels[name] = self._validate_entry(
                name,
                {
                    'filename': path.name,
                    'local_only': True,
                },
            )
            changes += 1

        if changes:
            self._save_metadata()
        return changes

    def add_kernel(
            self,
            name: str,
            url: str,
            overwrite: bool = False,
            *,
            sha256: str | None = None,
            size: int | None = None,
    ) -> bool:
        if not isinstance(name, str) or not name.strip():
            raise ValueError('Kernel names must be non-empty strings')
        normalized_name = name.casefold()
        if normalized_name in self.kernels and not overwrite:
            return False
        entry = {
            'url': url,
            'filename': filename_from_url(url),
            'sha256': sha256,
            'size': size,
        }
        self.kernels[normalized_name] = self._validate_entry(name, entry)
        self._save_metadata()
        return True

    def get_planetary_kernels(self, name: str, *, retries: int = 3) -> list[str]:
        if not isinstance(name, str) or not name.strip():
            raise ValueError('Kernel names must be non-empty strings')
        self.sync_local_files()
        normalized_name = name.casefold()
        if normalized_name not in self.kernels:
            raise ValueError(f"Kernel '{name}' not found in metadata.")

        paths = []
        for file_info in self._entry_files(self.kernels[normalized_name]):
            local_path = self.download_dir / file_info['filename']
            integrity = {
                'expected_size': file_info.get('size'),
                'expected_sha256': file_info.get('sha256'),
            }
            if not file_matches(local_path, **integrity):
                url = file_info.get('url')
                if url is None:
                    raise FileNotFoundError(
                        f"Local kernel '{name}' is no longer available at "
                        f'{local_path}'
                    )
                download_file(
                    url,
                    local_path,
                    session=self._session,
                    retries=retries,
                    description=local_path.name,
                    **integrity,
                )
            paths.append(str(local_path.resolve()))
        return paths

    def get_planetary_kernel(
            self,
            name: str,
            *,
            retries: int = 3,
    ) -> str | list[str]:
        """Return one path for legacy entries or all paths for multipart entries."""
        paths = self.get_planetary_kernels(name, retries=retries)
        return paths[0] if len(paths) == 1 else paths

    def list_kernels(
            self,
    ) -> list[tuple[str, str | list[str] | None]]:
        self.sync_local_files()
        kernels = []
        for name, data in sorted(self.kernels.items()):
            urls = [
                file_info['url']
                for file_info in self._entry_files(data)
                if 'url' in file_info
            ]
            source = urls[0] if len(urls) == 1 else urls or None
            kernels.append((name, source))
        return kernels

    def kernel_names(self) -> list[str]:
        """Return the names currently available in the metadata file."""
        self.sync_local_files()
        return sorted(self.kernels)

    def kernel_statuses(self) -> list[dict[str, Any]]:
        self.sync_local_files()
        statuses = []
        for name, info in sorted(self.kernels.items()):
            files = self._entry_files(info)
            paths = [self.download_dir / file_info['filename'] for file_info in files]
            urls = [
                file_info['url']
                for file_info in files
                if 'url' in file_info
            ]
            cached_files = [
                file_matches(
                    path,
                    expected_size=file_info.get('size'),
                    expected_sha256=file_info.get('sha256'),
                )
                for path, file_info in zip(paths, files)
            ]
            string_paths = [str(path) for path in paths]
            url_value = urls[0] if len(urls) == 1 else urls or None
            path_value = string_paths[0] if len(string_paths) == 1 else string_paths
            statuses.append(
                {
                    'name': name,
                    'urls': urls,
                    'url': url_value,
                    'local_only': all(
                        file_info.get('local_only', False)
                        for file_info in files
                    ),
                    'paths': string_paths,
                    'path': path_value,
                    'cached_files': cached_files,
                    'cached': all(cached_files),
                }
            )
        return statuses
