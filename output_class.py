class CastepOutput:
    def __init__(self, path) -> None:
        with open(path, 'r') as f:
            self.energy = 1