import sys


def progressbar(it, prefix="", size=40, file=sys.stdout):
    """Yield items while displaying a text progress bar.

    Parameters
    ----------
    it : iterable
        Iterable object with a defined length.
    prefix : `str`, optional
        Text printed before the progress bar. Default is an empty string.
    size : `int`, optional
        Width of the progress bar, in characters. Default is 40.
    file : file-like, optional
        Output stream where the progress bar is written. Default is
        `sys.stdout`.

    Yields
    ------
    item
        Items from `it`, in the original order.
    """
    count = len(it)

    def show(j):
        x = int(size*j/count)
        file.write("%s |%s%s|  - %i%% \r" % (prefix, "█"*x, "."*(size-x), 100*j/count))
        file.flush()

    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()


def progressbar_show(j, count, prefix='', size=40, file=sys.stdout):
    """Display a text progress bar for a given progress value.

    Parameters
    ----------
    j : `int`
        Current progress value.
    count : `int`
        Total number of steps.
    prefix : `str`, optional
        Text printed before the progress bar. Default is an empty string.
    size : `int`, optional
        Width of the progress bar, in characters. Default is 40.
    file : file-like, optional
        Output stream where the progress bar is written. Default is
        `sys.stdout`.
    """
    x = int(size*j/count)
    file.write("%s |%s%s|  - %i%% \r" % (prefix, "█"*x, "."*(size-x), 100*j/count))
    file.flush()
    if x >= size:
        file.write("\n")
        file.flush()
