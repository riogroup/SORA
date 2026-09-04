"""Reusable helpers for remote tabular data cached in SQLite databases."""

from __future__ import annotations

import hashlib
import json
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timedelta
from io import StringIO
from pathlib import Path
from typing import Any, Iterator, Mapping

import pandas as pd
import numpy as np
import requests
from tqdm import tqdm

from sora.config.core import get_config, resolve_storage_path
from sora.config.download import fetch_text, write_json_atomic


def load_json(json_file: Path | str) -> dict[str, Any]:
    """Load a JSON metadata object, returning an empty mapping if absent."""
    path = Path(json_file)
    if not path.exists():
        return {}
    try:
        with path.open('r', encoding='utf-8') as stream:
            document = json.load(stream)
    except json.JSONDecodeError as exc:
        raise ValueError(f'Invalid database metadata in {path}: {exc}') from exc
    if not isinstance(document, dict):
        raise TypeError(f'Database metadata in {path} must be a JSON object')
    return document


def save_json(file_path: Path | str, data: Mapping[str, Any]) -> None:
    """Persist database metadata using an atomic file replacement."""
    write_json_atomic(file_path, data)


class BaseDatabase:
    """Mirror remote CSV tables into a transactional local SQLite database.

    Subclasses declare a mapping from table names to source URLs in ``_urls``
    and may optionally describe foreign keys in ``_fk_map``.
    """

    _update_description = 'Updating database.'
    _fk_map: Mapping[str, list[dict[str, str]]] = {}
    _urls: Mapping[str, str] = {}

    def __init__(
            self,
            config_section_data,
            *,
            data_path: Path | str | None = None,
            session: requests.Session | None = None,
    ):
        """Initialize paths, request session, and persisted refresh metadata."""
        self._section_data = config_section_data
        self._data_path = Path(
            get_config().data_path if data_path is None else data_path
        )
        self._session = session
        self.conn: sqlite3.Connection | None = None
        self._load_config()

    def _load_config(self) -> None:
        """Resolve configured paths and load persisted refresh metadata."""
        self.db_path = resolve_storage_path(
            self._data_path,
            self._section_data.database,
        )
        self._json_file = resolve_storage_path(
            self._data_path,
            self._section_data.json_data,
        )
        self.data_dir = resolve_storage_path(
            self._data_path,
            self._section_data.data_dir,
        )
        self.update_age_days = self._section_data.update_age_days
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._json_file.parent.mkdir(parents=True, exist_ok=True)
        self._json_data = load_json(self._json_file)

    @contextmanager
    def connection(self) -> Iterator[sqlite3.Connection]:
        """Yield a transactional SQLite connection and always close it."""
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        connection = sqlite3.connect(self.db_path)
        connection.execute('PRAGMA foreign_keys = ON')
        try:
            yield connection
            connection.commit()
        except Exception:
            connection.rollback()
            raise
        finally:
            connection.close()

    def open_connection(self) -> sqlite3.Connection:
        """Open a compatibility connection for external callers."""
        connection = self.conn
        if connection is None:
            self.db_path.parent.mkdir(parents=True, exist_ok=True)
            connection = sqlite3.connect(self.db_path)
            connection.execute('PRAGMA foreign_keys = ON')
            self.conn = connection
        return connection

    def close_connection(self, *, commit: bool = True) -> None:
        """Close a compatibility connection, committing by default."""
        if self.conn is None:
            return
        try:
            if commit:
                self.conn.commit()
            else:
                self.conn.rollback()
        finally:
            self.conn.close()
            self.conn = None

    def _database_has_required_tables(self) -> bool:
        """Return whether the database contains every configured source table."""
        if not self.db_path.is_file():
            return False
        required_tables = set(self._urls)
        if not required_tables:
            return True
        try:
            with sqlite3.connect(self.db_path) as connection:
                existing = {
                    row[0]
                    for row in connection.execute(
                        "SELECT name FROM sqlite_master WHERE type = 'table'"
                    )
                }
        except sqlite3.DatabaseError:
            return False
        return required_tables.issubset(existing)

    @staticmethod
    def _table_exists(connection: sqlite3.Connection, table: str) -> bool:
        """Return whether a table exists in an open SQLite connection."""
        return connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
            (table,),
        ).fetchone() is not None

    def should_update(self) -> bool:
        """Return whether data is missing, stale, or uses different sources."""
        if not self._database_has_required_tables():
            return True
        stored_sources = self._json_data.get('sources')
        if not isinstance(stored_sources, Mapping):
            return True
        if dict(stored_sources) != dict(self._urls):
            return True
        last_update = self._json_data.get('last_update')
        if not isinstance(last_update, str) or not last_update:
            return True
        try:
            last = datetime.fromisoformat(last_update)
        except (TypeError, ValueError):
            return True
        return datetime.now() - last > timedelta(days=self.update_age_days)

    @staticmethod
    def _quote_identifier(identifier: str) -> str:
        """Quote a validated SQLite identifier without treating it as SQL."""
        if not isinstance(identifier, str) or not identifier:
            raise ValueError('SQL identifiers must be non-empty strings')
        return '"' + identifier.replace('"', '""') + '"'

    @staticmethod
    def _sqlite_value(value):
        """Convert pandas and NumPy values into SQLite-compatible objects."""
        if pd.isna(value):
            return None
        return value.item() if isinstance(value, np.generic) else value

    def _create_table(
            self,
            connection: sqlite3.Connection,
            name: str,
            dataframe: pd.DataFrame,
            foreign_keys=None,
    ) -> dict[str, str]:
        """Replace a table from a dataframe within the caller's transaction."""
        if dataframe.empty:
            raise ValueError(f'Remote table {name!r} is empty')
        columns = [str(column) for column in dataframe.columns]
        if not columns or len({column.casefold() for column in columns}) != len(columns):
            raise ValueError(f'Remote table {name!r} has invalid or duplicate columns')

        schema: dict[str, str] = {}
        declarations: list[str] = []
        for column in columns:
            series = dataframe[column]
            if pd.api.types.is_bool_dtype(series) or pd.api.types.is_integer_dtype(series):
                column_type = 'INTEGER'
            elif pd.api.types.is_float_dtype(series):
                column_type = 'REAL'
            else:
                column_type = 'TEXT'
            primary_key = ' PRIMARY KEY' if column.casefold() == 'id' else ''
            declarations.append(
                f'{self._quote_identifier(column)} {column_type}{primary_key}'
            )
            schema[column] = column_type

        for foreign_key in foreign_keys or ():
            declarations.append(
                'FOREIGN KEY ({source}) REFERENCES {table}({target}) '
                'ON DELETE CASCADE'.format(
                    source=self._quote_identifier(foreign_key['from']),
                    table=self._quote_identifier(foreign_key['table']),
                    target=self._quote_identifier(foreign_key['to']),
                )
            )

        quoted_table = self._quote_identifier(name)
        connection.execute(f'DROP TABLE IF EXISTS {quoted_table}')
        connection.execute(
            f'CREATE TABLE {quoted_table} ({", ".join(declarations)})'
        )

        quoted_columns = ', '.join(self._quote_identifier(column) for column in columns)
        placeholders = ', '.join('?' for _ in columns)
        rows = (
            tuple(self._sqlite_value(value) for value in row)
            for row in dataframe.itertuples(index=False, name=None)
        )
        connection.executemany(
            f'INSERT INTO {quoted_table} ({quoted_columns}) VALUES ({placeholders})',
            rows,
        )
        return schema

    def _get_data(self, url: str) -> str:
        """Fetch one configured remote table as text."""
        return fetch_text(url, session=self._session)

    def update_database(self, force: bool = False) -> bool:
        """Refresh changed remote tables and publish their update metadata.

        Parameters
        ----------
        force : `bool`, optional
            Ignore the configured maximum age and check every remote source.

        Returns
        -------
        `bool`
            Whether an update check was performed.
        """
        if not force and not self.should_update():
            return False

        stored_hashes = self._json_data.get('hash', {})
        stored_schema = self._json_data.get('schema_cache', {})
        current_hashes = (
            dict(stored_hashes) if isinstance(stored_hashes, Mapping) else {}
        )
        schema_cache = (
            dict(stored_schema) if isinstance(stored_schema, Mapping) else {}
        )
        new_hashes = current_hashes.copy()

        # Replace all changed tables in one transaction. Refresh metadata is
        # written only after that transaction commits successfully.
        with self.connection() as connection:
            for key, url in tqdm(
                    self._urls.items(),
                    desc=self._update_description,
            ):
                text = self._get_data(url)
                hash_value = hashlib.sha256(text.encode()).hexdigest()
                if (
                        current_hashes.get(key) == hash_value
                        and self._table_exists(connection, key)
                ):
                    continue

                dataframe = pd.read_csv(StringIO(text), sep=',', quotechar='"')
                schema_cache[key] = self._create_table(
                    connection,
                    key,
                    dataframe,
                    self._fk_map.get(key),
                )
                new_hashes[key] = hash_value

        new_data = {
            'sources': dict(self._urls),
            'hash': new_hashes,
            'schema_cache': schema_cache,
            'last_update': datetime.now().isoformat(),
        }
        save_json(self._json_file, new_data)
        self._json_data = new_data
        return True

    def _ensure_column(self, table: str, colname: str, coltype: str = 'TEXT') -> None:
        """Add a column when absent while restricting its declared type."""
        normalized_type = coltype.upper()
        if normalized_type not in {'TEXT', 'INTEGER', 'REAL', 'BLOB', 'NUMERIC'}:
            raise ValueError(f'Unsupported SQLite column type: {coltype!r}')
        quoted_table = self._quote_identifier(table)
        quoted_column = self._quote_identifier(colname)
        with self.connection() as connection:
            columns = {
                row[1]
                for row in connection.execute(f'PRAGMA table_info({quoted_table})')
            }
            if not columns:
                raise sqlite3.OperationalError(f'no such table: {table}')
            if colname in columns:
                return
            connection.execute(
                f'ALTER TABLE {quoted_table} ADD COLUMN '
                f'{quoted_column} {normalized_type}'
            )

    def status(self) -> dict[str, Any]:
        """Return storage locations, freshness, and row counts by table."""
        tables: dict[str, int] = {}
        if self.db_path.is_file():
            try:
                with sqlite3.connect(self.db_path) as connection:
                    for table in self._urls:
                        quoted_table = self._quote_identifier(table)
                        if connection.execute(
                                "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
                                (table,),
                        ).fetchone():
                            tables[table] = connection.execute(
                                f'SELECT COUNT(*) FROM {quoted_table}'
                            ).fetchone()[0]
            except sqlite3.DatabaseError:
                tables = {}
        return {
            'database': str(self.db_path),
            'metadata': str(self._json_file),
            'data_directory': str(self.data_dir),
            'sources': dict(self._urls),
            'database_exists': self.db_path.is_file(),
            'last_update': self._json_data.get('last_update'),
            'needs_update': self.should_update(),
            'tables': tables,
        }
