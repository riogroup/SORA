from __future__ import annotations

from math import isfinite
from typing import Any, Mapping
from urllib.parse import urlparse

from .meta import BaseConfigSection

__all__ = [
    'BodyPlotConfig',
    'OccMapConfig',
    'PredictionConfig',
    'ServicesConfig',
]


def _is_positive_number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and isfinite(value)
        and value > 0
    )


def _is_non_negative_number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and isfinite(value)
        and value >= 0
    )


class ServicesConfig(BaseConfigSection):
    """Connection preferences for external services used by SORA."""

    FIELDS = ('linea_tap_url', 'sbdb_cache')
    LOCAL_KEYS = frozenset(FIELDS)

    def _initialize(
        self,
        default_data: Mapping[str, Any],
        effective_data: Mapping[str, Any],
    ) -> None:
        self.linea_tap_url = effective_data['linea_tap_url']
        self.sbdb_cache = effective_data['sbdb_cache']

    def _validate(self) -> None:
        if not isinstance(self.linea_tap_url, str):
            raise TypeError('services.linea_tap_url must be a string')

        parsed_url = urlparse(self.linea_tap_url)
        if parsed_url.scheme not in {'http', 'https'} or not parsed_url.netloc:
            raise ValueError(
                'services.linea_tap_url must be an absolute HTTP or HTTPS URL'
            )

        if not isinstance(self.sbdb_cache, bool):
            raise TypeError('services.sbdb_cache must be a boolean')


class BodyPlotConfig(BaseConfigSection):
    """Visual defaults for plots of body shape models."""

    FIELDS = (
        'show_pole',
        'north_pole_color',
        'south_pole_color',
        'pole_length_scale',
        'default_surface_color',
        'surface_alpha',
    )
    LOCAL_KEYS = frozenset(FIELDS)

    def _initialize(
        self,
        default_data: Mapping[str, Any],
        effective_data: Mapping[str, Any],
    ) -> None:
        self.show_pole = effective_data['show_pole']
        self.north_pole_color = effective_data['north_pole_color']
        self.south_pole_color = effective_data['south_pole_color']
        self.pole_length_scale = effective_data['pole_length_scale']
        self.default_surface_color = effective_data['default_surface_color']
        self.surface_alpha = effective_data['surface_alpha']

    def _validate(self) -> None:
        if not isinstance(self.show_pole, bool):
            raise TypeError('body_plot.show_pole must be a boolean')

        for field in ('north_pole_color', 'south_pole_color'):
            value = getattr(self, field)
            if not isinstance(value, str) or not value.strip():
                raise TypeError(f'body_plot.{field} must be a non-empty string')

        if not _is_positive_number(self.pole_length_scale):
            raise ValueError('body_plot.pole_length_scale must be positive')

        if (
            not isinstance(self.default_surface_color, (list, tuple))
            or len(self.default_surface_color) != 3
            or not all(
                _is_non_negative_number(value) and value <= 1
                for value in self.default_surface_color
            )
        ):
            raise ValueError(
                'body_plot.default_surface_color must contain three numbers '
                'between 0 and 1'
            )

        if (
            not _is_non_negative_number(self.surface_alpha)
            or self.surface_alpha > 1
        ):
            raise ValueError('body_plot.surface_alpha must be between 0 and 1')


class PredictionConfig(BaseConfigSection):
    """Operational defaults for catalogue-based occultation predictions."""

    FIELDS = (
        'default_catalogue',
        'catalogue_row_limit',
        'catalogue_timeout',
        'fallback_catalogue',
        'fallback_on_timeout',
    )
    LOCAL_KEYS = frozenset(FIELDS)

    def _initialize(
        self,
        default_data: Mapping[str, Any],
        effective_data: Mapping[str, Any],
    ) -> None:
        self.default_catalogue = effective_data['default_catalogue']
        self.catalogue_row_limit = effective_data['catalogue_row_limit']
        self.catalogue_timeout = effective_data['catalogue_timeout']
        self.fallback_catalogue = effective_data['fallback_catalogue']
        self.fallback_on_timeout = effective_data['fallback_on_timeout']

    def _validate(self) -> None:
        for field in ('default_catalogue', 'fallback_catalogue'):
            value = getattr(self, field)
            if not isinstance(value, str) or not value.strip():
                raise TypeError(f'prediction.{field} must be a non-empty string')

        if (
            not isinstance(self.catalogue_row_limit, int)
            or isinstance(self.catalogue_row_limit, bool)
            or self.catalogue_row_limit == 0
            or self.catalogue_row_limit < -1
        ):
            raise ValueError(
                'prediction.catalogue_row_limit must be -1 or a positive integer'
            )

        if not _is_positive_number(self.catalogue_timeout):
            raise ValueError('prediction.catalogue_timeout must be positive')

        if not isinstance(self.fallback_on_timeout, bool):
            raise TypeError('prediction.fallback_on_timeout must be a boolean')


class OccMapConfig(BaseConfigSection):
    """User preferences for generated occultation maps."""

    FIELDS = (
        'style',
        'resolution',
        'format',
        'dpi',
        'size_cm',
        'meridian_interval',
        'parallel_interval',
        'center_point_interval',
        'site_name_scale',
        'country_name_scale',
        'site_marker_scale',
        'center_marker_scale',
        'positive_site_color',
        'negative_site_color',
        'site_marker',
        'site_name_offset_km',
        'states',
        'labels',
        'site_names',
        'arrow',
        'night_alpha',
        'site_box_alpha',
    )
    LOCAL_KEYS = frozenset(FIELDS)

    def _initialize(
        self,
        default_data: Mapping[str, Any],
        effective_data: Mapping[str, Any],
    ) -> None:
        self.style = effective_data['style']
        self.resolution = effective_data['resolution']
        self.format = effective_data['format']
        self.dpi = effective_data['dpi']
        self.size_cm = effective_data['size_cm']
        self.meridian_interval = effective_data['meridian_interval']
        self.parallel_interval = effective_data['parallel_interval']
        self.center_point_interval = effective_data['center_point_interval']
        self.site_name_scale = effective_data['site_name_scale']
        self.country_name_scale = effective_data['country_name_scale']
        self.site_marker_scale = effective_data['site_marker_scale']
        self.center_marker_scale = effective_data['center_marker_scale']
        self.positive_site_color = effective_data['positive_site_color']
        self.negative_site_color = effective_data['negative_site_color']
        self.site_marker = effective_data['site_marker']
        self.site_name_offset_km = effective_data['site_name_offset_km']
        self.states = effective_data['states']
        self.labels = effective_data['labels']
        self.site_names = effective_data['site_names']
        self.arrow = effective_data['arrow']
        self.night_alpha = effective_data['night_alpha']
        self.site_box_alpha = effective_data['site_box_alpha']

    def _validate(self) -> None:
        if not isinstance(self.style, int) or isinstance(self.style, bool):
            raise TypeError('occ_map.style must be an integer')
        if self.style not in {1, 2}:
            raise ValueError('occ_map.style must be 1 or 2')

        if not isinstance(self.resolution, int) or isinstance(
            self.resolution, bool
        ):
            raise TypeError('occ_map.resolution must be an integer')
        if self.resolution not in {1, 2, 3}:
            raise ValueError('occ_map.resolution must be 1, 2, or 3')

        if not isinstance(self.format, str) or not self.format.strip():
            raise TypeError('occ_map.format must be a non-empty string')

        if (
            not isinstance(self.dpi, int)
            or isinstance(self.dpi, bool)
            or self.dpi <= 0
        ):
            raise ValueError('occ_map.dpi must be a positive integer')

        if (
            not isinstance(self.size_cm, (list, tuple))
            or len(self.size_cm) != 2
            or not all(_is_positive_number(value) for value in self.size_cm)
        ):
            raise ValueError(
                'occ_map.size_cm must contain two positive numbers'
            )

        for field in ('meridian_interval', 'parallel_interval'):
            if not _is_positive_number(getattr(self, field)):
                raise ValueError(f'occ_map.{field} must be positive')

        if (
            not isinstance(self.center_point_interval, int)
            or isinstance(self.center_point_interval, bool)
            or self.center_point_interval <= 0
        ):
            raise ValueError(
                'occ_map.center_point_interval must be a positive integer'
            )

        for field in (
            'site_name_scale',
            'country_name_scale',
            'site_marker_scale',
            'center_marker_scale',
        ):
            if not _is_non_negative_number(getattr(self, field)):
                raise ValueError(f'occ_map.{field} must be non-negative')

        for field in ('positive_site_color', 'negative_site_color', 'site_marker'):
            value = getattr(self, field)
            if not isinstance(value, str) or not value.strip():
                raise TypeError(f'occ_map.{field} must be a non-empty string')

        if (
            not isinstance(self.site_name_offset_km, (list, tuple))
            or len(self.site_name_offset_km) != 2
            or not all(
                isinstance(value, (int, float))
                and not isinstance(value, bool)
                and isfinite(value)
                for value in self.site_name_offset_km
            )
        ):
            raise ValueError(
                'occ_map.site_name_offset_km must contain two finite numbers'
            )

        for field in ('states', 'labels', 'site_names', 'arrow'):
            if not isinstance(getattr(self, field), bool):
                raise TypeError(f'occ_map.{field} must be a boolean')

        for field in ('night_alpha', 'site_box_alpha'):
            value = getattr(self, field)
            if not _is_non_negative_number(value) or value > 1:
                raise ValueError(f'occ_map.{field} must be between 0 and 1')


DEFAULT_SECTION_TYPES = {
    'services': ServicesConfig,
    'prediction': PredictionConfig,
    'occ_map': OccMapConfig,
    'body_plot': BodyPlotConfig,
}
