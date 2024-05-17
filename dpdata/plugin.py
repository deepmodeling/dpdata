"""Base of plugin systems."""

from __future__ import annotations


class Plugin:
    """A class to register plugins.

    Examples
    --------
    >>> example_plugin = Plugin()
    >>> @example_plugin.register("xx")
        def xxx():
            pass
    >>> print(example_plugin.plugins['xx'])
    """

    def __init__(self):
        self.plugins = {}

    def register(self, key):
        """Register a plugin.

        Parameters
        ----------
        key : str
            Key of the plugin.
        """

        def decorator(object):
            self.plugins[key] = object
            return object

        return decorator

    def get_plugin(self, key):
        return self.plugins[key]

    def __add__(self, other):
        self.plugins.update(other.plugins)
        return self
