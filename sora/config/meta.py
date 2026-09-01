from __future__ import annotations

from abc import ABC, abstractmethod
from copy import deepcopy
from pathlib import Path
from typing import Any, ClassVar, Mapping, TYPE_CHECKING

import yaml

if TYPE_CHECKING:
    from .core import Config

__all__ = [
    'BaseConfigSection',
    'deep_merge',
]


def deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    """Recursively merge mappings without mutating either input."""
    if not isinstance(base, Mapping):
        raise TypeError('base must be a mapping')
    if not isinstance(override, Mapping):
        raise TypeError('override must be a mapping')

    result = deepcopy(dict(base))
    for key, value in override.items():
        if (
            key in result
            and isinstance(result[key], Mapping)
            and isinstance(value, Mapping)
        ):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result


class BaseConfigSection(ABC):
    """Base class for a typed, explicitly declared configuration section.

    Subclasses must list every serializable attribute in ``FIELDS`` and the
    subset users may override in ``LOCAL_KEYS``. They are responsible for
    assigning all fields in ``_initialize`` and may validate the resulting
    values in ``_validate``. Assigning a declared field after initialization
    validates and persists it as a user override.
    """

    FIELDS: ClassVar[tuple[str, ...]] = ()
    LOCAL_KEYS: ClassVar[frozenset[str]] = frozenset()

    def __setattr__(self, name: str, value: Any) -> None:
        """Persist assignments to declared configuration fields."""
        if (
            name in type(self).FIELDS
            and '_parent' in self.__dict__
            and not self.__dict__.get('_applying_field_values', False)
        ):
            self._parent.update(self._section_name, name, value)
            return
        object.__setattr__(self, name, value)

    def __init__(
        self,
        parent: Config,
        section_name: str,
        default_data: Mapping[str, Any],
        local_data: Mapping[str, Any],
    ):
        self._parent = parent
        self._section_name = section_name
        self._validate_declaration()
        self._default_data = self._validate_default_data(default_data)
        self._local_overrides: dict[str, Any] = {}
        self._applying_field_values = False
        self._apply_local_overrides(local_data)

    def _validate_declaration(self) -> None:
        invalid_fields = [
            field
            for field in self.FIELDS
            if not isinstance(field, str) or not field
        ]
        if invalid_fields:
            raise TypeError(
                f"Configuration section '{self._section_name}' field names "
                'must be non-empty strings.'
            )
        if len(self.FIELDS) != len(set(self.FIELDS)):
            raise TypeError(
                f"Configuration section '{self._section_name}' declares duplicate fields."
            )

        local_keys = set(self.LOCAL_KEYS)
        if any(not isinstance(key, str) or not key for key in local_keys):
            raise TypeError(
                f"Configuration section '{self._section_name}' LOCAL_KEYS "
                'must contain only non-empty strings.'
            )

        invalid_local_keys = local_keys.difference(self.FIELDS)
        if invalid_local_keys:
            keys = ', '.join(sorted(invalid_local_keys))
            raise TypeError(
                f"Configuration section '{self._section_name}' declares unknown "
                f'LOCAL_KEYS: {keys}.'
            )

    def _validate_default_data(
        self, default_data: Mapping[str, Any]
    ) -> dict[str, Any]:
        if not isinstance(default_data, Mapping):
            raise TypeError(
                f"Default configuration section '{self._section_name}' must be a mapping."
            )
        if any(not isinstance(key, str) for key in default_data):
            raise TypeError(
                f"Default configuration section '{self._section_name}' "
                'must use string keys.'
            )

        unknown = set(default_data).difference(self.FIELDS)
        if unknown:
            keys = ', '.join(sorted(unknown))
            raise KeyError(
                f"Unknown default configuration keys in '{self._section_name}': {keys}."
            )
        return deepcopy(dict(default_data))

    def _validate_local_data(
        self, local_data: Mapping[str, Any]
    ) -> dict[str, Any]:
        if not isinstance(local_data, Mapping):
            raise TypeError(
                f"User configuration section '{self._section_name}' must be a mapping."
            )
        if any(not isinstance(key, str) for key in local_data):
            raise TypeError(
                f"User configuration section '{self._section_name}' "
                'must use string keys.'
            )

        local_keys = set(local_data)
        unknown = local_keys.difference(self.FIELDS)
        if unknown:
            keys = ', '.join(sorted(unknown))
            raise KeyError(
                f"Unknown user configuration keys in '{self._section_name}': {keys}."
            )

        forbidden = local_keys.difference(self.LOCAL_KEYS)
        if forbidden:
            keys = ', '.join(sorted(forbidden))
            raise KeyError(
                f"Configuration keys are not user-overridable in "
                f"'{self._section_name}': {keys}."
            )
        return deepcopy(dict(local_data))

    def _apply_local_overrides(self, local_data: Mapping[str, Any]) -> None:
        local_overrides = self._validate_local_data(local_data)
        effective_data = deep_merge(self._default_data, local_overrides)
        previous_values = {
            field: deepcopy(getattr(self, field))
            for field in self.FIELDS
            if hasattr(self, field)
        }
        previous_overrides = deepcopy(self._local_overrides)

        self._applying_field_values = True
        try:
            for field in self.FIELDS:
                if hasattr(self, field):
                    delattr(self, field)
            self._initialize(self._default_data, effective_data)

            missing = [
                field for field in self.FIELDS if not hasattr(self, field)
            ]
            if missing:
                keys = ', '.join(missing)
                raise TypeError(
                    f"Configuration section '{self._section_name}' did not "
                    f'initialize fields: {keys}.'
                )

            self._validate()
            self._local_overrides = local_overrides
        except Exception:
            for field in self.FIELDS:
                if hasattr(self, field):
                    delattr(self, field)
            for field, value in previous_values.items():
                setattr(self, field, value)
            self._local_overrides = previous_overrides
            raise
        finally:
            self._applying_field_values = False

    @abstractmethod
    def _initialize(
        self,
        default_data: Mapping[str, Any],
        effective_data: Mapping[str, Any],
    ) -> None:
        """Assign every field from the default and effective configuration."""
        raise NotImplementedError

    def _validate(self) -> None:
        """Validate initialized values. Subclasses may override this hook."""

    @staticmethod
    def _serialize(value: Any) -> Any:
        if isinstance(value, BaseConfigSection):
            return value.to_dict()
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, Mapping):
            return {
                key: BaseConfigSection._serialize(item)
                for key, item in value.items()
            }
        if isinstance(value, tuple):
            return [BaseConfigSection._serialize(item) for item in value]
        if isinstance(value, list):
            return [BaseConfigSection._serialize(item) for item in value]
        return deepcopy(value)

    def to_dict(self) -> dict[str, Any]:
        """Return all explicitly declared effective values."""
        return {
            field: self._serialize(getattr(self, field))
            for field in self.FIELDS
        }

    def to_dict_local(self) -> dict[str, Any]:
        """Return only values explicitly overridden by the user."""
        return {
            field: self._serialize(self._local_overrides[field])
            for field in self.FIELDS
            if field in self._local_overrides
        }

    def set_local(self, key: str, value: Any) -> None:
        """Set and validate one user-overridable value in memory."""
        if key not in self.FIELDS:
            raise KeyError(
                f"Unknown configuration key '{self._section_name}.{key}'."
            )
        if key not in self.LOCAL_KEYS:
            raise KeyError(
                f"Configuration key '{self._section_name}.{key}' "
                'is not user-overridable.'
            )

        local_overrides = deepcopy(self._local_overrides)
        local_overrides[key] = value
        self._apply_local_overrides(local_overrides)

    def _snapshot_local_overrides(self) -> dict[str, Any]:
        """Return an internal, non-serialized snapshot for rollback."""
        return deepcopy(self._local_overrides)

    def remove_local(self, key: str) -> None:
        """Remove a user override and restore the corresponding default."""
        if key not in self.FIELDS:
            raise KeyError(
                f"Unknown configuration key '{self._section_name}.{key}'."
            )
        if key not in self.LOCAL_KEYS:
            raise KeyError(
                f"Configuration key '{self._section_name}.{key}' "
                'is not user-overridable.'
            )

        local_overrides = deepcopy(self._local_overrides)
        local_overrides.pop(key, None)
        self._apply_local_overrides(local_overrides)

    def save(self) -> None:
        """Persist all current user overrides through the parent config."""
        self._parent.save_local()

    def __repr__(self) -> str:
        return f'<{self.__class__.__name__}({self.to_dict()})>'

    def __str__(self) -> str:
        return yaml.safe_dump(
            self.to_dict(),
            sort_keys=False,
            default_flow_style=False,
        )
