from dpdata.format import Format


@Format.register("list")
class ListFormat(Format):
    def to_system(self, data, **kwargs):
        """Convert system to list, usefull for data collection."""
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
