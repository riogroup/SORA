# remove this block for v1.0
class _PositionDict(dict):
    """This is a modified Dictionary object to allow switching on/off of data points.
    It also avoids user to change data.
    """

    def __setitem__(self, key, value):
        """ Redefines how to set a value to a key in the dictionary.
        It only sets a value if the key starts with '_occ_'. Otherwise, it only allows for the user to provide
        'on' or 'off' which is passed only to change the 'on' keyword.
        """
        status = {'on': True, 'off': False}
        n = 0
        if key.startswith('_occ_'):
            super().__setitem__(key[5:], value)
            n = 1
        elif key in self.keys():
            n = 1
            if value not in status.keys():
                raise ValueError("Value must be 'on' or 'off' only.")
            if type(self[key]) == _PositionDict:
                for k in self[key].keys():
                    self[key][k] = value
            elif key == 'on':
                super().__setitem__('on', status[value])
                if status[value]:
                    self['enable']()
                else:
                    self['disable']()
        else:
            if value not in status.keys():
                raise ValueError("Value must be 'on' or 'off' only.")
            for key1 in self.keys():
                if type(self[key1]) == _PositionDict:
                    if key in self[key1].keys():
                        n = 1
                        self[key1][key] = value
        if n == 0:
            raise KeyError('Key "{}" does not exist'.format(key))

    def __str__(self):
        out = []
        for key in self.keys():
            if key not in ['enable', 'disable']:
                out.append('{}: {}'.format(key, self[key]))
        out = '\n' + '\n'.join(out)
        return out.replace('\n', '\n  ')

    def __repr__(self):
        out = []
        for key in self.keys():
            if key not in ['enable', 'disable']:
                out.append('\'{}\': {}'.format(key, self[key].__repr__()))
        out = ',\n'.join(out)
        return '{' + out.replace('\n', '\n  ') + '}'
# end of block removal
