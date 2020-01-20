import importlib


def find(search_string):
    """
    Performs a find operation on available functions.

    Parameters
    ----------
    search_string : str
        Name of the function to be searched.

    Returns
    -------
    list
        A list of available functions related to ``search_string``.

    """
    search_string = search_string.lower()
    current_module = "einsteinpy.symbolic.predefined"
    predefined = importlib.import_module(current_module)
    return [
        fname
        for fname in dir(predefined)
        if ((search_string in fname.lower()) and (not fname.islower()))
    ]
