class Stoichiometry(dict):
    """
    Extends dict to return 0 when getting an item that is not a key in the dict and to support subtraction
    """
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            return 0

    def __sub__(self, other: "Stoichiometry"):
        keys = set(self.keys()).union(other.keys())
        result = {}
        for k in keys:
            result[k] = self[k] - other[k]
        return Stoichiometry(result)