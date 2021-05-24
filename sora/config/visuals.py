import sys


def progressbar(it, prefix="", size=40, file=sys.stdout):
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
    x = int(size*j/count)
    file.write("%s |%s%s|  - %i%% \r" % (prefix, "█"*x, "."*(size-x), 100*j/count))
    file.flush()
    if x >= size:
        file.write("\n")
        file.flush()
