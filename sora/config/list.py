from collections import OrderedDict

__all__ = ['List']


class List(OrderedDict):
    """Abstract Class to handle SORA List Classes.

    This Class inherits from OrderedDict, which is a dictionary that keeps the
    order that the objects were added. With this, the user can access a item by
    its order as well.

    """
    _allowed_types = None  # _allowed_types attribute must be defined for each class with the list of allowed types.
    _set_func = "_add_item"  # _set_func attribute must be defined for each class with the name of the function to add item.

    def __setitem__(self, name, item):
        """Overwrites the default __setitem__ from OrderedDict.

        The reason is to forbid the user to define an item as ``obj[name] = item``
        and force the use of the function 'add_item'.
        """
        raise ValueError("{} cannot be set directly. Please use {}".format(self.__class__.__name__, self._set_func))

    # __add_item must be defined inside an add method of Child Class.
    def _add_item(self, name, item):
        """Adds an item to the list with the name given.

        Parameters
        ----------
        name : `str`
            The name of the item. It must be unique for a item in the list.

        item : `_allowed_types`
            The item to be added to the list. Its type must be one of the types
            in the attribute ``self._allowed_type``.
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
        """"Defines behaviour for ``obj[name]``.

        Parameters
        ----------
        key : `str`, `int`
            The key used to find item in the list. It can the name given in
            `_add_item` or the number corresponding to the order the item is
            stored in the list.
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
        """Defines behaviour of ``del(obj[item])`` to delete only one item in the list.

        Parameters
        ----------
        key : `str`, `int`
            The key used to find item in the list. It can the name given in
            ``_add_item`` or the number corresponding to the order the item is
            stored in the list.
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
        """Defines behaviour for str(obj) or print(obj)
        """
        strings = []
        for key in self.keys():
            strings.append(str(self[key]))
        return "\n".join(strings)

    def __repr__(self):
        """Defines the object's string representation.

        <ClassListName:
            0: Type of item0 (Name of item0)
            1: Type of item1 (Name of item1)
            ...>
        """
        data = ['\n    {}: {}({})'.format(i, self[i].__class__.__name__, self[i].name) for i in range(len(self))]
        return '<{}:{}>'.format(self.__class__.__name__, ''.join(data))
