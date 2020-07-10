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

