# Useful classes
class colors():
    '''
    Docstring
    '''
    __cores = {
    'positive': 'blue',
    'negative': 'green',
    'error': 'red'
    }
    def __init__(self):
        pass

    @property
    def positive_color(self):
        return self.__cores['positive']

    @positive_color.setter
    def positive_color(self, value):
        self.__cores['positive']=value

    @property
    def negative_color(self):
        return self.__cores['negative']

    @negative_color.setter
    def negative_color(self, value):
        self.__cores['negative']=value

    @property
    def error_bar(self):
        return self.__cores['error']

    @error_bar.setter
    def error_bar(self, value):
        self.__cores['error']=value


# useful functions
def test_attr(attr, typ, name):
    """
    This function tests if the attribute "attr" belongs to the type "typ",
    if not it gives an error message informing the name of the variable
    """
    try:
        return typ(attr)
    except:
        raise TypeError('"{}" keyword must be a {} object'.format(name, typ))

