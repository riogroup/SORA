### Useful classes
class colors():
    '''
    Docstring
    '''
    positive_color = 'blue'
    negative_color = 'green'
    error_bar = 'red'
    def __init__(self):
        return
    
### useful functions
def test_attr(attr, typ, name):
    """
    This function tests if the attr belongs to the type typ,
    if not it gives an error message.
    """
    try:
        return typ(attr)
    except:
        raise ValueError('{} keyword must be a {} object'.format(name, typ))