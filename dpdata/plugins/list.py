from __future__ import annotations

from dpdata.format import Format


@Format.register("list")
class ListFormat(Format):
    """In-memory list of one-frame System objects.

    This write-only convenience format splits a multi-frame System or
    LabeledSystem into independent one-frame objects; it does not create an
    on-disk file.
    """

    def to_system(self, data, **kwargs):
        """Split system data into a list of one-frame systems.

        Parameters
        ----------
        data : dict
            System or labeled-system data.
        **kwargs : dict
            Additional format arguments accepted for API compatibility.

        Returns
        -------
        list[System] or list[LabeledSystem]
            Empty for zero frames; otherwise one object per frame.
        """
        from dpdata import LabeledSystem, System

        if "forces" in data:
            system = LabeledSystem(data=data)
        else:
            system = System(data=data)
        if len(system) == 0:
            return []
        if len(system) == 1:
            return [system]
        else:
            systems = []
            for ii in range(len(system)):
                systems.append(system.sub_system([ii]))
            return systems
