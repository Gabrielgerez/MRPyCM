class StepFunction:
    def __init__(self, cavity, inside=1.0, outside=2.0):
        self.C = cavity
        self.inside = inside  # permittivity of free space
        self.outside = outside  # permittivity of solvent

    def __call__(self):
        raise NotImplementedError()
