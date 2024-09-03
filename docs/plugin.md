# Plugins

One can follow a simple example under `plugin_example/` directory to add their own format by creating and installing plugins.
It's critical to add the {class}`Format` class to `entry_points['dpdata.plugins']` in `pyproject.toml`:

```toml
[project.entry-points.'dpdata.plugins']
random = "dpdata_random:RandomFormat"
```
