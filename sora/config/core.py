from __future__ import annotations

import os
import tempfile
from pathlib import Path
from threading import RLock
from types import MappingProxyType
from typing import Any, ClassVar, Mapping

import yaml
from platformdirs import PlatformDirs

from .meta import BaseConfigSection
from .sections import DEFAULT_SECTION_TYPES

__all__ = ['CONFIG_VERSION', 'Config', 'get_config', 'reload_config']

APP_NAME = 'sora'
CONFIG_NAME = 'config.yaml'
CONFIG_VERSION = 1
_VERSION_KEY = 'config_version'


class Config:
    """Load and persist registered SORA configuration sections.

    Concrete configuration sections are implemented and registered separately.
    This class only owns paths, YAML loading, version checks, section
    construction, and persistence of explicit user overrides.

    The distributed defaults must declare the current ``config_version``.
    A local override file may omit it and is then interpreted using the current
    schema; an explicitly declared incompatible version is rejected.
    """

    SECTION_TYPES: ClassVar[Mapping[str, type[BaseConfigSection]]] = (
        DEFAULT_SECTION_TYPES
    )

    def __init__(
        self,
        default_path: Path | str | None = None,
        *,
        local_path: Path | str | None = None,
        data_path: Path | str | None = None,
        section_types: Mapping[str, type[BaseConfigSection]] | None = None,
    ):
        platform_paths = PlatformDirs(APP_NAME, ensure_exists=False)
        self._default_path = (
            Path(default_path)
            if default_path is not None
            else Path(__file__).with_name('default_config.yaml')
        )
        self._local_path = (
            Path(local_path)
            if local_path is not None
            else platform_paths.user_config_path / CONFIG_NAME
        )
        self._data_path = (
            Path(data_path)
            if data_path is not None
            else platform_paths.user_data_path
        )
        self._section_types = self._prepare_section_types(section_types)

        default_document = self._load_yaml(self._default_path, required=True)
        local_document = self._load_yaml(self._local_path, required=False)
        self._validate_document(
            default_document,
            path=self._default_path,
            require_version=True,
            require_sections=True,
        )
        self._validate_document(
            local_document,
            path=self._local_path,
            require_version=False,
            require_sections=False,
        )

        self._sections: dict[str, BaseConfigSection] = {}
        for section_name, section_type in self._section_types.items():
            section = section_type(
                self,
                section_name,
                default_document[section_name],
                local_document.get(section_name, {}),
            )
            self._sections[section_name] = section
            setattr(self, section_name, section)

    @property
    def config_path(self) -> Path:
        """Path to the user override file."""
        return self._local_path

    @property
    def data_path(self) -> Path:
        """Platform-specific directory for persistent SORA data."""
        return self._data_path

    @property
    def default_path(self) -> Path:
        """Path to the distributed default configuration."""
        return self._default_path

    @property
    def sections(self) -> Mapping[str, BaseConfigSection]:
        """Read-only mapping of registered configuration sections."""
        return MappingProxyType(self._sections)

    def _prepare_section_types(
        self,
        section_types: Mapping[str, type[BaseConfigSection]] | None,
    ) -> dict[str, type[BaseConfigSection]]:
        registry = dict(
            self.SECTION_TYPES if section_types is None else section_types
        )
        reserved_names = {_VERSION_KEY, *dir(type(self)), *self.__dict__}
        for section_name, section_type in registry.items():
            if not isinstance(section_name, str) or not section_name:
                raise TypeError('Configuration section names must be non-empty strings.')
            if section_name in reserved_names:
                raise ValueError(
                    f"Configuration section name '{section_name}' is reserved."
                )
            if (
                not isinstance(section_type, type)
                or not issubclass(section_type, BaseConfigSection)
            ):
                raise TypeError(
                    f"Configuration section '{section_name}' must use a "
                    'BaseConfigSection subclass.'
                )
        return registry

    @staticmethod
    def _load_yaml(path: Path, *, required: bool) -> dict[str, Any]:
        if not path.exists():
            if required:
                raise FileNotFoundError(f'Configuration file not found: {path}')
            return {}

        try:
            with path.open('r', encoding='utf-8') as stream:
                document = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            raise ValueError(f'Invalid YAML in configuration file {path}: {exc}') from exc

        if document is None:
            return {}
        if not isinstance(document, Mapping):
            raise TypeError(
                f'Configuration file {path} must contain a mapping at its root.'
            )
        return dict(document)

    def _validate_document(
        self,
        document: Mapping[str, Any],
        *,
        path: Path,
        require_version: bool,
        require_sections: bool,
    ) -> None:
        if any(not isinstance(key, str) for key in document):
            raise TypeError(
                f'Configuration file {path} must use string keys at its root.'
            )

        has_version = _VERSION_KEY in document
        if require_version and not has_version:
            raise ValueError(
                f'Missing configuration version in {path}: '
                f'expected {CONFIG_VERSION}.'
            )
        if has_version:
            version = document[_VERSION_KEY]
            if version != CONFIG_VERSION:
                raise ValueError(
                    f'Unsupported configuration version in {path}: '
                    f'expected {CONFIG_VERSION}, found {version!r}.'
                )

        known_keys = {_VERSION_KEY, *self._section_types}
        unknown = set(document).difference(known_keys)
        if unknown:
            keys = ', '.join(sorted(unknown))
            raise KeyError(f'Unknown configuration sections in {path}: {keys}.')

        if require_sections:
            missing = set(self._section_types).difference(document)
            if missing:
                keys = ', '.join(sorted(missing))
                raise KeyError(
                    f'Missing default configuration sections in {path}: {keys}.'
                )

        for section_name in self._section_types:
            if section_name in document and not isinstance(
                document[section_name], Mapping
            ):
                raise TypeError(
                    f"Configuration section '{section_name}' in {path} "
                    'must be a mapping.'
                )

    def to_dict(self) -> dict[str, Any]:
        """Return the complete effective configuration."""
        return {
            _VERSION_KEY: CONFIG_VERSION,
            **{
                section_name: section.to_dict()
                for section_name, section in self._sections.items()
            },
        }

    def to_local_dict(self) -> dict[str, Any]:
        """Return only explicit user overrides and the format version."""
        document: dict[str, Any] = {_VERSION_KEY: CONFIG_VERSION}
        for section_name, section in self._sections.items():
            local_values = section.to_dict_local()
            if local_values:
                document[section_name] = local_values
        return document

    def save_local(self) -> None:
        """Atomically persist explicit user overrides."""
        self._local_path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f'.{self._local_path.name}.',
            suffix='.tmp',
            dir=self._local_path.parent,
            text=True,
        )
        temporary_path = Path(temporary_name)

        try:
            with os.fdopen(descriptor, 'w', encoding='utf-8') as stream:
                yaml.safe_dump(
                    self.to_local_dict(),
                    stream,
                    sort_keys=False,
                    default_flow_style=False,
                )
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary_path, self._local_path)
        except Exception:
            temporary_path.unlink(missing_ok=True)
            raise

    save = save_local

    def update(self, section: str, key: str, value: Any) -> None:
        """Set, validate, and persist one user-overridable value."""
        section_object = self._get_section(section)
        previous_overrides = section_object._snapshot_local_overrides()
        section_object.set_local(key, value)
        try:
            self.save_local()
        except Exception:
            section_object._apply_local_overrides(previous_overrides)
            raise

    def remove_override(self, section: str, key: str) -> None:
        """Remove and persist one user override."""
        section_object = self._get_section(section)
        previous_overrides = section_object._snapshot_local_overrides()
        section_object.remove_local(key)
        try:
            self.save_local()
        except Exception:
            section_object._apply_local_overrides(previous_overrides)
            raise

    def _get_section(self, section: str) -> BaseConfigSection:
        try:
            return self._sections[section]
        except KeyError as exc:
            raise KeyError(f"Unknown configuration section '{section}'.") from exc

    def __str__(self) -> str:
        return yaml.safe_dump(
            self.to_dict(),
            sort_keys=False,
            default_flow_style=False,
        )

    def __repr__(self) -> str:
        return (
            f'<Config(default={self._default_path!s}, '
            f'local={self._local_path!s})>'
        )


_active_config: Config | None = None
_active_config_lock = RLock()


def get_config() -> Config:
    """Return the shared process-wide configuration instance.

    The configuration is loaded lazily on the first call. Application code
    should use this function instead of constructing ``Config`` directly.
    """
    global _active_config

    with _active_config_lock:
        if _active_config is None:
            _active_config = Config()
        return _active_config


def reload_config() -> Config:
    """Reload and replace the shared configuration instance.

    Existing references still point to the previous instance; subsequent
    calls to ``get_config`` return the newly loaded configuration.
    """
    global _active_config

    with _active_config_lock:
        _active_config = Config()
        return _active_config
