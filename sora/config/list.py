from collections import OrderedDict

__all__ = ['List']


class List(OrderedDict):
    """Base class for SORA ordered list containers.

    This class inherits from `collections.OrderedDict`, which keeps the order
    in which objects were added. Items can be accessed by name or by insertion
    order.

    """
    _allowed_types = None  # _allowed_types attribute must be defined for each class with the list of allowed types.
    _set_func = "_add_item"  # _set_func attribute must be defined for each class with the name of the function to add item.

    def __setitem__(self, name, item):
        """Forbid direct assignment with ``obj[name] = item``.

        Parameters
        ----------
        name : `str`
            Name that would be assigned to the item.
        item
            Item that would be stored.

        Raises
        ------
        ValueError
            Always raised to force the use of the configured add method.
        """
        raise ValueError("{} cannot be set directly. Please use {}".format(self.__class__.__name__, self._set_func))

    # __add_item must be defined inside an add method of Child Class.
    def _add_item(self, name, item):
        """Add an item to the list with the given name.

        Parameters
        ----------
        name : `str`
            Name of the item. It must be unique in the list.

        item : `_allowed_types`
            The item to be added to the list. Its type must be one of the types
            in the attribute ``self._allowed_types``.

        Raises
        ------
        TypeError
            If `name` is not a string or `item` is not an allowed type.
        ValueError
            If `name` is empty or already exists in the list.
        """
        if not isinstance(name, str):
            raise TypeError("name must be a string")
        if name == '':
            raise ValueError("name can not be an empty string")
        if name in self.keys():
            raise ValueError("{} already exists. It must be removed before redefined".format(name))
        if not isinstance(item, self._allowed_types):
            raise TypeError("item must be an allowed type: {}".format([k.__name__ for k in self._allowed_types]))
        super().__setitem__(name, item)

    def __getitem__(self, key):
        """Return an item by name or insertion order.

        Parameters
        ----------
        key : `str`, `int`
            Key used to find an item in the list. It can be the name given in
            ``_add_item`` or the index corresponding to the order in which the
            item is stored in the list.

        Returns
        -------
        item
            Stored item matching `key`.

        Raises
        ------
        IndexError
            If an integer key is outside the list bounds.
        TypeError
            If `key` is neither a string nor an integer.
        """
        if isinstance(key, int):
            try:
                key = list(self.keys())[key]
            except IndexError:
                raise IndexError("Cannot get item {} for list with size {}".format(key, len(self)))
            return self[key]
        elif isinstance(key, str):
            return super().__getitem__(key)
        else:
            raise TypeError("{} can only be indexed with the "
                            "named keys and integers.".format(self.__class__.__name__))

    def __delitem__(self, key):
        """Delete one item from the list by name or insertion order.

        Parameters
        ----------
        key : `str`, `int`
            Key used to find an item in the list. It can be the name given in
            ``_add_item`` or the index corresponding to the order in which the
            item is stored in the list.

        Raises
        ------
        IndexError
            If an integer key is outside the list bounds.
        TypeError
            If `key` is neither a string nor an integer.
        """
        if isinstance(key, int):
            try:
                key = list(self.keys())[key]
            except IndexError:
                raise IndexError("Cannot get item {} for list with size {}".format(key, len(self)))
        if not isinstance(key, str):
            raise TypeError("{} can only be indexed with the "
                            "named keys and integers.".format(self.__class__.__name__))
        super().__delitem__(key)

    def __str__(self):
        """Return a printable representation of all items in the list."""
        strings = []
        for key in self.keys():
            strings.append(str(self[key]))
        return "\n".join(strings)

    def __repr__(self):
        """Return the developer representation of the list.

        <ClassListName:
            0: Type of item0 (Name of item0)
            1: Type of item1 (Name of item1)
            ...>
        """
        data = ['\n    {}: {}({})'.format(i, self[i].__class__.__name__, self[i].name) for i in range(len(self))]
        return '<{}:{}>'.format(self.__class__.__name__, ''.join(data))
