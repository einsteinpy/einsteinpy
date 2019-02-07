import alabaster

import einsteinpy

project = "EinsteinPy"
copyright = "2019, EinsteinPy Development Team"

version = einsteinpy.__version__

release = version
highlight_language = "python"
pygments_style = "sphinx"

add_function_parentheses = True

add_module_names = True

needs_sphinx = "1.3"
extensions = ["alabaster"]

templates_path = ["_templates"]

source_suffix = ".rst"

master_doc = "index"

html_theme = "alabaster"

html_theme_path = [alabaster.get_path()]

html_title = "EinsteinPy"

html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        "relations.html",
        "searchbox.html",
        "donate.html",
    ]
}
