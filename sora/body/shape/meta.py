from abc import ABC, abstractmethod


class BaseShape(ABC):
    """Abstract base class for shape models."""

    def __init__(self) -> None:
        super(BaseShape, self).__init__()

    def __repr__(self) -> str:
        """Return the unambiguous string representation of the shape."""
        return '<{}>'.format(self.__str__())

    def __str__(self) -> str:
        """Return the display name of the shape."""
        return '{}: {}'.format(self.__class__.__name__, getattr(self, 'name', ''))

    @abstractmethod
    def get_limb(self):
        """Return the projected limb of the shape."""
        pass
