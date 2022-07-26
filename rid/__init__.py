__author__ = 'Yanze Wang; Jiahao Fan'
__author_email__ = 'yanze039@mit.edu; jiahaofan@pku.edu.cn'
try:
    from ._version import version as __version__
except ImportError:
    from .__about__ import __version__

SRC_ROOT = __path__[0]
