# -*- coding: utf-8 -*-

__version__ = '1.0.dev'
__date__ = '2016/01/13'

try:
    from ._githash import __githash__, __dev_value__
    if '.dev' in __version__:
        __version__ += __dev_value__
except Exception:
    pass
