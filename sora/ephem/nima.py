"""Local NIMA catalogue and BSP cache management."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import pandas as pd
import requests
from tqdm import tqdm

from sora.config.core import Config, get_config
from sora.config.database import BaseDatabase
from sora.config.download import download_file, filename_from_url, file_matches

__all__ = ['NimaDB']


class NimaDB(BaseDatabase):
    """Manage the local NIMA catalogue and its associated BSP files."""

    _update_description = 'Updating NIMA database.'

    def __init__(
            self,
            config: Config | None = None,
            *,
            session: requests.Session | None = None,
    ):
        """Initialize NIMA storage from the active SORA configuration."""
        self.config = config or get_config()
        self.database_url = self.config.nima.database_url
        self._urls = {'nima_table': self.database_url}
        super().__init__(
            self.config.nima,
            data_path=self.config.data_path,
            session=session,
        )
        self.sync_local_files()

    def show_table(
            self,
            limit: int | None = None,
            filters: Mapping[str, Any] | None = None,
    ) -> pd.DataFrame:
        """Return rows from the NIMA table using parameterized filters.

        Parameters
        ----------
        limit : `int`, optional
            Maximum number of rows to return.
        filters : mapping, optional
            Exact column-value matches combined with ``AND``.

        Returns
        -------
        `pandas.DataFrame`
            Rows selected from the locally cached NIMA catalogue.
        """
        if limit is not None and (
                not isinstance(limit, int) or isinstance(limit, bool) or limit <= 0
        ):
            raise ValueError('limit must be a positive integer')
        self.update_database()

        query = 'SELECT * FROM "nima_table"'
        parameters: list[Any] = []
        with self.connection() as connection:
            columns = {
                row[1]
                for row in connection.execute('PRAGMA table_info("nima_table")')
            }
            clauses: list[str] = []
            for column, value in (filters or {}).items():
                if column not in columns:
                    raise KeyError(f'Unknown NIMA table column: {column}')
                clauses.append(f'{self._quote_identifier(column)} = ?')
                parameters.append(value)
            if clauses:
                query += ' WHERE ' + ' AND '.join(clauses)
            if limit is not None:
                query += ' LIMIT ?'
                parameters.append(limit)
            return pd.read_sql_query(query, connection, params=parameters)

    @staticmethod
    def _normalized_name(value: Any) -> str:
        """Normalize an object identifier for NIMA catalogue lookup."""
        return ''.join(str(value).split()).casefold()

    @staticmethod
    def _local_candidate(data_dir: Path, bsp_url: str) -> Path:
        """Return the expected local path for a catalogue BSP URL."""
        return data_dir / filename_from_url(bsp_url)

    def _find_bsp_record(self, name: Any):
        """Find the BSP URL, SPK ID, and cached path for an object."""
        normalized_name = self._normalized_name(name)
        query = """
                SELECT bspfile, idSPK, bspfilepath
                FROM nima_table
                WHERE (
                    REPLACE(LOWER(CAST(name AS TEXT)), ' ', '') = ?
                        OR REPLACE(LOWER(CAST(designation AS TEXT)), ' ', '') = ?
                        OR REPLACE(LOWER(CAST(number AS TEXT)), ' ', '') = ?
                    )
                  AND bspfile IS NOT NULL \
                """
        with self.connection() as connection:
            return connection.execute(
                query,
                (normalized_name, normalized_name, normalized_name),
            ).fetchone()

    def _record_local_path(self, idspk, path: Path) -> None:
        """Persist the resolved local BSP path for one catalogue row."""
        with self.connection() as connection:
            connection.execute(
                'UPDATE nima_table SET bspfilepath = ? WHERE idSPK = ?',
                (str(path.resolve()), idspk),
            )

    def _path_in_data_dir(self, path: Path) -> bool:
        """Return whether a path belongs to the configured BSP directory."""
        try:
            path.resolve().relative_to(self.data_dir.resolve())
        except ValueError:
            return False
        return True

    def sync_local_files(self) -> int:
        """Associate manually copied BSP files with their NIMA catalogue rows."""
        if not self._database_has_required_tables():
            return 0

        self.data_dir.mkdir(parents=True, exist_ok=True)
        self._ensure_column('nima_table', 'bspfilepath')
        local_files = {
            path.name.casefold(): path
            for path in sorted(
                self.data_dir.iterdir(),
                key=lambda candidate: candidate.name.casefold(),
            )
            if file_matches(path)
        }
        with self.connection() as connection:
            rows = connection.execute(
                'SELECT bspfile, idSPK, bspfilepath '
                'FROM nima_table WHERE bspfile IS NOT NULL'
            ).fetchall()

            updates: list[tuple[str | None, Any]] = []
            for bsp_url, idspk, stored_path in rows:
                try:
                    filename = filename_from_url(bsp_url)
                except (TypeError, ValueError):
                    continue
                local_path = local_files.get(filename.casefold())
                resolved_path = (
                    str(local_path.resolve()) if local_path is not None else None
                )
                if stored_path != resolved_path:
                    updates.append((resolved_path, idspk))

            if updates:
                connection.executemany(
                    'UPDATE nima_table SET bspfilepath = ? WHERE idSPK = ?',
                    updates,
                )
        return len(updates)

    def update_database(self, force: bool = False) -> bool:
        """Refresh the catalogue and synchronize manually copied BSP files."""
        catalogue_updated = super().update_database(force=force)
        return bool(self.sync_local_files()) or catalogue_updated

    def status(self) -> dict[str, Any]:
        """Return catalogue, storage, and locally indexed BSP status."""
        synchronized = self.sync_local_files()
        status = super().status()
        indexed_files = 0
        if self._database_has_required_tables():
            with self.connection() as connection:
                indexed_files = connection.execute(
                    'SELECT COUNT(*) FROM nima_table '
                    'WHERE bspfilepath IS NOT NULL'
                ).fetchone()[0]
        status.update(
            {
                'indexed_bsp_files': indexed_files,
                'synchronized_paths': synchronized,
            }
        )
        return status

    def download_all_bspfiles(self, retries: int = 3) -> list[Path]:
        """Download every BSP referenced by the NIMA catalogue.

        Existing valid files are reused. If one or more downloads fail, all
        successful paths are recorded before a summary ``RuntimeError`` is
        raised.
        """
        self.update_database()
        self._ensure_column('nima_table', 'bspfilepath')
        self.data_dir.mkdir(parents=True, exist_ok=True)

        with self.connection() as connection:
            rows = connection.execute(
                'SELECT name, bspfile, idSPK, bspfilepath '
                'FROM nima_table WHERE bspfile IS NOT NULL'
            ).fetchall()

        downloaded: list[Path] = []
        failures: list[str] = []
        updates: list[tuple[str, Any]] = []
        for name, bsp_url, idspk, stored_path in tqdm(
                rows,
                desc='NIMA: Downloading BSP files',
        ):
            try:
                configured_path = self._local_candidate(self.data_dir, bsp_url)
                stored_candidate = Path(stored_path) if stored_path else None
                if (
                        stored_candidate is not None
                        and self._path_in_data_dir(stored_candidate)
                        and file_matches(stored_candidate)
                ):
                    local_path = stored_candidate
                elif file_matches(configured_path):
                    local_path = configured_path
                else:
                    local_path = download_file(
                        bsp_url,
                        configured_path,
                        session=self._session,
                        retries=retries,
                        description=str(name),
                    )
                    downloaded.append(local_path)
                updates.append((str(local_path.resolve()), idspk))
            except Exception as exc:
                failures.append(f'{name}: {exc}')

        if updates:
            with self.connection() as connection:
                connection.executemany(
                    'UPDATE nima_table SET bspfilepath = ? WHERE idSPK = ?',
                    updates,
                )
        if failures:
            raise RuntimeError(
                f'{len(failures)} NIMA BSP download(s) failed: '
                + '; '.join(failures)
            )
        return downloaded

    def get_bspfile(self, name: Any) -> tuple[str, Any]:
        """Return the local NIMA BSP path and SPK ID for an object.

        The catalogue is refreshed when stale, and a missing local BSP is
        downloaded into the currently configured data directory.
        """
        self.update_database()
        self._ensure_column('nima_table', 'bspfilepath')

        row = self._find_bsp_record(name)
        if row is None:
            normalized_name = self._normalized_name(name)
            raise ValueError(f'NIMA: Object {normalized_name} not found on database')

        bsp_url, idspk, stored_path = row
        self.data_dir.mkdir(parents=True, exist_ok=True)
        configured_path = self._local_candidate(self.data_dir, bsp_url)
        stored_candidate = Path(stored_path) if stored_path else None

        # Ignore a path recorded for an older destination after data_dir changes.
        if (
                stored_candidate is not None
                and self._path_in_data_dir(stored_candidate)
                and file_matches(stored_candidate)
        ):
            return str(stored_candidate.resolve()), idspk
        if file_matches(configured_path):
            self._record_local_path(idspk, configured_path)
            return str(configured_path.resolve()), idspk

        local_path = download_file(
            bsp_url,
            configured_path,
            session=self._session,
            description=configured_path.name,
        )
        self._record_local_path(idspk, local_path)
        return str(local_path.resolve()), idspk
