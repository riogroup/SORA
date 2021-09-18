from abc import ABC, abstractmethod


class BaseShape(ABC):
    def __init__(self) -> None:
        super(BaseShape, self).__init__()

    def __repr__(self) -> str:
        return '<{}: {}>'.format(self.__class__.__name__, getattr(self, 'name', ''))

    def __str__(self) -> str:
        return super(BaseShape, self).__str__()

    @abstractmethod
    def get_limb(self):
        pass
